#include "ForceDiagramStorage.h"
#include "LinkNodesOnDevice.h"
#include "WLCSolveOnDevice.h"
#include "TorsionSolveOnDevice.h"
#include "PltForceDevice.h"
#include "PltForceFieldDevice.h"
#include "PltonPltForceFieldDevice.h"
#include "PltTndrlDevice.h"
#include "PltonPltTndrlDevice.h"
#include "PltVlmPushDevice.h"
#include "AdvancePositionOnDevice.h"
#include "BucketSchemeOnDevice.h"
#include "NodeSystemDevice.h"



void NodeSystemDevice::setBucketScheme() {

	initDimensionBucketScheme(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);

	buildBucketScheme(nodeInfoVecs, pltInfoVecs, domainParams,
		auxVecs, generalParams);

	extendBucketScheme(nodeInfoVecs, pltInfoVecs, domainParams, auxVecs);
};

void NodeSystemDevice::solveForcesOnDevice() {

	//RESET FORCE TO ZERO AT BEGINNING/////////////////////////////////////////////////
	thrust::fill(nodeInfoVecs.nodeForceX.begin(),nodeInfoVecs.nodeForceX.end(),0);
	thrust::fill(nodeInfoVecs.nodeForceY.begin(),nodeInfoVecs.nodeForceY.end(),0);
	thrust::fill(nodeInfoVecs.nodeForceZ.begin(),nodeInfoVecs.nodeForceZ.end(),0);

	if (generalParams.linking == true) {
			LinkNodesOnDevice(
					nodeInfoVecs,
					wlcInfoVecs,
					auxVecs,
					generalParams);
	}
	TorsionSolveOnDevice(nodeInfoVecs, torsionInfoVecs, generalParams);

	//std::cout<<"prewlc"<<std::endl;
	WLCSolveOnDevice(nodeInfoVecs, wlcInfoVecs, generalParams);

	//platetelet-node forces
	//RESETS PLATELET FORCES
	if (generalParams.pltfrcfld == true) {// note: this force-field includes both pulling and pushing
		PltForceFieldOnDevice(//platelet on node force field
			nodeInfoVecs,
			wlcInfoVecs,
			generalParams,
			pltInfoVecs,
			auxVecs);
		if (generalParams.pltonplt == true) {
			PltInteractionPltOnDevice(//platelet on platelet interaction through force field
				generalParams,
				pltInfoVecs,
				auxVecs);
		}

	}
	else if (generalParams.plttndrl == true) { //note for now force-field type has priority over tndrl-type
		//initialize Trnl-Node Id list
	  if (generalParams.currentTime==0.0){
	    thrust::fill(pltInfoVecs.tndrlNodeId.begin(),pltInfoVecs.tndrlNodeId.end(), generalParams.maxIdCount);
			thrust::fill(pltInfoVecs.tndrlNodeType.begin(),pltInfoVecs.tndrlNodeType.end(), 0);
	    }

		// Tndrl-node pulling
		PltTndrlOnDevice(
		  nodeInfoVecs,
		  wlcInfoVecs,
		  generalParams,
		  pltInfoVecs,
		  auxVecs);

		//Tndrl-Plt pulling
		if (generalParams.pltonplt == true) {
			PltonPltTndrlOnDevice(//platelet on platelet interaction through tndrl
				generalParams,
				pltInfoVecs,
				auxVecs);
		}

		PltVlmPushOnDevice(//push for volume exclusion
			nodeInfoVecs,
			wlcInfoVecs,
			generalParams,
			pltInfoVecs,
			auxVecs);

	}




};


void NodeSystemDevice::solveSystemDevice() {


	//set initial epsilon
	generalParams.epsilon = (1.0) *
		sqrt(6.0*generalParams.kB * generalParams.temperature * generalParams.dtTemp / generalParams.viscousDamp);

	while (generalParams.runSim == true) {
		generalParams.iterationCounter++;
		generalParams.currentTime += generalParams.dtTemp;



		AdvancePositionOnDevice(
			nodeInfoVecs,
			pltInfoVecs,
		 	generalParams);

		setBucketScheme();

		solveForcesOnDevice(); //resets and solves forces for next time step


		if (generalParams.iterationCounter % 50 == 0) {
			storage->print_VTK_File();
			//store sum of all forces on each node. Used in stress calculations
			//store before upadting storage class.
			thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())) + generalParams.maxNodeCount,
				nodeInfoVecs.sumForcesOnNode.begin(),//save vector
				NormFunctor());
				//platelets
				thrust::transform(
					thrust::make_zip_iterator(
						thrust::make_tuple(
							pltInfoVecs.pltForceX.begin(),
							pltInfoVecs.pltForceY.begin(),
							pltInfoVecs.pltForceZ.begin())),
					thrust::make_zip_iterator(
						thrust::make_tuple(
							pltInfoVecs.pltForceX.begin(),
							pltInfoVecs.pltForceY.begin(),
							pltInfoVecs.pltForceZ.begin())) + generalParams.maxPltCount,
					pltInfoVecs.sumForcesOnPlt.begin(),//save vector
					NormFunctor());


			generalParams.epsilon = (1.0) *
				sqrt(6.0 * generalParams.kB * generalParams.temperature * generalParams.dtTemp / generalParams.viscousDamp);
		}

	}

};

NodeSystemDevice::NodeSystemDevice()  {};

void NodeSystemDevice::assignForceDiagramStorage(std::shared_ptr<ForceDiagramStorage> _storage) {
	storage = _storage;
}

//__host__ __device__
void NodeSystemDevice::initializeSystem(
	thrust::host_vector<bool>& hostIsNodeFixed,
	thrust::host_vector<double>& hostPosX,
	thrust::host_vector<double>& hostPosY,
	thrust::host_vector<double>& hostPosZ,
	thrust::host_vector<unsigned>& hostWLCEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCEdgeRight,
	thrust::host_vector<double>& hostWLCLenZero,

	thrust::host_vector<unsigned>& hostWLCSubEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCSubEdgeRight,
	thrust::host_vector<double>& hostWLCSubLenZero,
	thrust::host_vector<unsigned>& hostTorsionIndexLeft,
	thrust::host_vector<unsigned>& hostTorsionIndexCenter,
	thrust::host_vector<unsigned>& hostTorsionIndexRight,
	thrust::host_vector<double>& hostTorsionAngleZero,
	//platelets
	thrust::host_vector<bool>& hostIsPltFixed,
	thrust::host_vector<double>& hostPltPosX,
	thrust::host_vector<double>& hostPltPosY,
	thrust::host_vector<double>& hostPltPosZ) {

	std::cout<< "total Edge Count: "<< generalParams.originEdgeCount << std::endl;
	std::cout << "max num nodes: " << generalParams.maxNodeCount << std::endl;
	//platelets

	std::cout << "max num platelets in device: " << generalParams.maxPltCount << std::endl;



	setPltVecs(
		hostIsPltFixed,
		hostPltPosX,
		hostPltPosY,
		hostPltPosZ);

	setNodeVecs(//calls initDimensionBucketScheme
		hostIsNodeFixed,
		hostPosX,
		hostPosY,
		hostPosZ);


	setTorsionVecs(
		hostTorsionIndexLeft,
		hostTorsionIndexCenter,
		hostTorsionIndexRight,
		hostTorsionAngleZero);

	setWLCVecs(
		hostWLCEdgeLeft,
		hostWLCEdgeRight,
		hostWLCLenZero );

};

void NodeSystemDevice::setNodeVecs(
	thrust::host_vector<bool>& hostIsNodeFixed,
	thrust::host_vector<double>& hostPosX,
	thrust::host_vector<double>& hostPosY,
	thrust::host_vector<double>& hostPosZ) {


	nodeInfoVecs.sumForcesOnNode.resize(generalParams.maxNodeCount);

	nodeInfoVecs.nodeVelocity.resize(generalParams.maxNodeCount);

	nodeInfoVecs.nodeLocX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeLocY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeLocZ.resize(generalParams.maxNodeCount);

	nodeInfoVecs.nodeForceX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeForceY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeForceZ.resize(generalParams.maxNodeCount);

	nodeInfoVecs.discretizedEdgeStrain.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	nodeInfoVecs.discretizedEdgeAlignment.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);

	//sized larger for input later
	nodeInfoVecs.deviceEdgeLeft.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	nodeInfoVecs.deviceEdgeRight.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);


	thrust::fill(nodeInfoVecs.discretizedEdgeStrain.begin(), nodeInfoVecs.discretizedEdgeStrain.end(),0.0);
	thrust::fill(nodeInfoVecs.deviceEdgeRight.begin(), nodeInfoVecs.deviceEdgeRight.end(), 0);	//fill force and velocity with zeros for computation.
	thrust::fill(nodeInfoVecs.deviceEdgeLeft.begin(), nodeInfoVecs.deviceEdgeLeft.end(), 0);	//fill force and velocity with zeros for computation.

	thrust::fill(nodeInfoVecs.sumForcesOnNode.begin(), nodeInfoVecs.sumForcesOnNode.end(), 0);


	thrust::copy(hostPosX.begin(), hostPosX.end(), nodeInfoVecs.nodeLocX.begin());
	thrust::copy(hostPosY.begin(), hostPosY.end(), nodeInfoVecs.nodeLocY.begin());
	thrust::copy(hostPosZ.begin(), hostPosZ.end(), nodeInfoVecs.nodeLocZ.begin());


	//copy fixed positions
	nodeInfoVecs.isNodeFixed.resize(generalParams.maxNodeCount);
	thrust::copy(hostIsNodeFixed.begin(), hostIsNodeFixed.end(), nodeInfoVecs.isNodeFixed.begin());

	nodeInfoVecs.linksThreadMade.resize(generalParams.maxNodeCount);
	nodeInfoVecs.delinksThreadMade.resize(generalParams.maxNodeCount);
	nodeInfoVecs.idMadeTempLeft.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);
	nodeInfoVecs.idMadeTempRight.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);

	//at this point all nodes are filled, so we can generate domainParams
	initDimensionBucketScheme(
		nodeInfoVecs,
		pltInfoVecs,
		domainParams,
		auxVecs,
		generalParams);


	domainParams.originMinX = domainParams.minX;
	domainParams.originMaxX = domainParams.maxX;
	domainParams.originMinY = domainParams.minY;
	domainParams.originMaxY = domainParams.maxY;
	domainParams.originMinZ = domainParams.minZ;
	domainParams.originMaxZ = domainParams.maxZ;

	std::cout<< "node count : " <<nodeInfoVecs.nodeLocY.size()<< std::endl;


	auxVecs.id_bucket.resize(generalParams.maxNodeCount);
	auxVecs.id_value.resize(generalParams.maxNodeCount);
	auxVecs.id_bucket_expanded.resize(27 * (generalParams.maxNodeCount));
	auxVecs.id_value_expanded.resize(27 *( generalParams.maxNodeCount ));

};

//platelet
void NodeSystemDevice::setPltVecs(
	thrust::host_vector<bool>& hostIsPltFixed,
	thrust::host_vector<double>& hostPltPosX,
	thrust::host_vector<double>& hostPltPosY,
	thrust::host_vector<double>& hostPltPosZ) {


	pltInfoVecs.sumForcesOnPlt.resize(generalParams.maxPltCount);

	pltInfoVecs.pltVelocity.resize(generalParams.maxPltCount);

	pltInfoVecs.pltLocX.resize(generalParams.maxPltCount);
	pltInfoVecs.pltLocY.resize(generalParams.maxPltCount);
	pltInfoVecs.pltLocZ.resize(generalParams.maxPltCount);

	pltInfoVecs.pltForceX.resize(generalParams.maxPltCount);
	pltInfoVecs.pltForceY.resize(generalParams.maxPltCount);
	pltInfoVecs.pltForceZ.resize(generalParams.maxPltCount);

	pltInfoVecs.pltImagingConnection.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
	pltInfoVecs.nodeImagingConnection.resize(generalParams.maxPltCount * generalParams.pltMaxConn);

	pltInfoVecs.nodeUnreducedId.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
	pltInfoVecs.nodeUnreducedForceX.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
	pltInfoVecs.nodeUnreducedForceY.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
	pltInfoVecs.nodeUnreducedForceZ.resize(generalParams.maxPltCount * generalParams.pltMaxConn);

	pltInfoVecs.nodeReducedId.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
	pltInfoVecs.nodeReducedForceX.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
	pltInfoVecs.nodeReducedForceY.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
	pltInfoVecs.nodeReducedForceZ.resize(generalParams.maxPltCount * generalParams.pltMaxConn);

	thrust::fill(pltInfoVecs.sumForcesOnPlt.begin(), pltInfoVecs.sumForcesOnPlt.end(), 0);


	thrust::copy(hostPltPosX.begin(), hostPltPosX.end(), pltInfoVecs.pltLocX.begin());
	thrust::copy(hostPltPosY.begin(), hostPltPosY.end(), pltInfoVecs.pltLocY.begin());
	thrust::copy(hostPltPosZ.begin(), hostPltPosZ.end(), pltInfoVecs.pltLocZ.begin());


	std::cout<<"num platelets: "<< pltInfoVecs.pltLocX.size() << std::endl;
	std::cout<<"num platelets var: "<< generalParams.maxPltCount << std::endl;
	//copy fixed positions
	pltInfoVecs.isPltFixed.resize(generalParams.maxPltCount);
	thrust::fill(pltInfoVecs.isPltFixed.begin(), pltInfoVecs.isPltFixed.end(), false);
	//thrust::copy(hostIsPltFixed.begin(), hostIsPltFixed.end(), pltInfoVecs.isPltFixed.begin());


	auxVecs.idPlt_bucket.resize(generalParams.maxPltCount);
	auxVecs.idPlt_value.resize(generalParams.maxPltCount);
	auxVecs.idPlt_bucket_expanded.resize(27 *( generalParams.maxPltCount ));
	auxVecs.idPlt_value_expanded.resize(27 * (generalParams.maxPltCount));

	if (generalParams.currentTime==0.0){
    pltInfoVecs.tndrlNodeId.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
		pltInfoVecs.tndrlNodeType.resize(generalParams.maxPltCount * generalParams.pltMaxConn);
  }

};

void NodeSystemDevice::setTorsionVecs(
	thrust::host_vector<unsigned>& hostTorsionIndexLeft,
	thrust::host_vector<unsigned>& hostTorsionIndexCenter,
	thrust::host_vector<unsigned>& hostTorsionIndexRight,
	thrust::host_vector<double>& hostTorsionAngleZero) {


	torsionInfoVecs.leftIndex.resize(generalParams.totalTorsionCount);
	torsionInfoVecs.centerIndex.resize(generalParams.totalTorsionCount);
	torsionInfoVecs.rightIndex.resize(generalParams.totalTorsionCount);
	torsionInfoVecs.angleZero.resize(generalParams.totalTorsionCount);

	thrust::copy(hostTorsionIndexLeft.begin(), hostTorsionIndexLeft.end(), torsionInfoVecs.leftIndex.begin());
	thrust::copy(hostTorsionIndexCenter.begin(), hostTorsionIndexCenter.end(), torsionInfoVecs.centerIndex.begin());
	thrust::copy(hostTorsionIndexRight.begin(), hostTorsionIndexRight.end(), torsionInfoVecs.rightIndex.begin());

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				torsionInfoVecs.leftIndex.begin(),
				torsionInfoVecs.centerIndex.begin(),
				torsionInfoVecs.rightIndex.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				torsionInfoVecs.leftIndex.begin(),
				torsionInfoVecs.centerIndex.begin(),
				torsionInfoVecs.rightIndex.begin())) + generalParams.totalTorsionCount,
			torsionInfoVecs.angleZero.begin(),//save vector
		TorsionAngleFunctor(
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data())));

	torsionInfoVecs.forceX.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.forceY.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.forceZ.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceX.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceY.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceZ.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);

	thrust::fill(torsionInfoVecs.forceX.begin(), torsionInfoVecs.forceX.end(), 0.0);
	thrust::fill(torsionInfoVecs.forceY.begin(), torsionInfoVecs.forceY.end(), 0.0);
	thrust::fill(torsionInfoVecs.forceZ.begin(), torsionInfoVecs.forceZ.end(), 0.0);

	torsionInfoVecs.tempTorIndices.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
	torsionInfoVecs.reducedIds.resize(torsionInfoVecs.factor * generalParams.totalTorsionCount);
};

void NodeSystemDevice::setWLCVecs(
	thrust::host_vector<unsigned>& hostWLCSubEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCSubEdgeRight,
	thrust::host_vector<double>& hostWLCSubLenZero ) {

	wlcInfoVecs.globalNeighbors.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	wlcInfoVecs.currentNodeEdgeCountVector.resize(generalParams.maxNodeCount);

	wlcInfoVecs.lengthZero.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	wlcInfoVecs.numOriginalNeighborsNodeVector.resize(generalParams.maxNodeCount);

	//default value is maxNodeCount
	thrust::fill(wlcInfoVecs.globalNeighbors.begin(), wlcInfoVecs.globalNeighbors.end(), generalParams.maxNodeCount);
	thrust::fill(wlcInfoVecs.currentNodeEdgeCountVector.begin(), wlcInfoVecs.currentNodeEdgeCountVector.end(),0);
	thrust::fill(wlcInfoVecs.lengthZero.begin(), wlcInfoVecs.lengthZero.end(), 0.0);



	nodeInfoVecs.deviceEdgeLeft = hostWLCSubEdgeLeft;
	nodeInfoVecs.deviceEdgeRight = hostWLCSubEdgeRight;

	//scan through hostAdj and put in device.
	for (unsigned id = 0; id < hostWLCSubLenZero.size(); id++) {

		unsigned idL = hostWLCSubEdgeLeft[id];
		unsigned idR = hostWLCSubEdgeRight[id];

		double edgeLen = hostWLCSubLenZero[id];
		//we use the lengthZero vector to identify edges as well.
		//node id is row, column node is connected to row node.

		//add edge for left node
		unsigned edgeNumL = wlcInfoVecs.currentNodeEdgeCountVector[idL]; //number of edges on (nodeId = row)	is that entry in cECV
		unsigned indexL = idL*generalParams.maxNeighborCount + edgeNumL;
		wlcInfoVecs.lengthZero[indexL] = edgeLen;
		wlcInfoVecs.globalNeighbors[indexL] = idR;
		(wlcInfoVecs.currentNodeEdgeCountVector[idL])++; //right connects to left

		//add edge for right node
		unsigned edgeNumR = wlcInfoVecs.currentNodeEdgeCountVector[idR]; //number of edges on (nodeId = row)	is that entry in cECV
		unsigned indexR = idR*generalParams.maxNeighborCount + edgeNumR;
		wlcInfoVecs.lengthZero[indexR] = edgeLen;
		wlcInfoVecs.globalNeighbors[indexR] = idL;
		(wlcInfoVecs.currentNodeEdgeCountVector[idR])++; //left connects to right

		generalParams.currentEdgeCount++;

	}

	//at this point currentNodeEdgeCountVector holds the number of edges, copy this to
	thrust::copy(wlcInfoVecs.currentNodeEdgeCountVector.begin(), wlcInfoVecs.currentNodeEdgeCountVector.end(), wlcInfoVecs.numOriginalNeighborsNodeVector.begin());
};
