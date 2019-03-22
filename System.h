#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "SystemStructures.h"


//Data Structure for node location. velocity and force
struct NodeInfoVecs {
	thrust::device_vector<unsigned> deviceEdgeLeft;
	thrust::device_vector<unsigned> deviceEdgeRight;


	thrust::device_vector<unsigned> idMadeTempLeft;//left id linked
	thrust::device_vector<unsigned> idMadeTempRight;//right id linked
	thrust::device_vector<unsigned> linksThreadMade;//max links each thread has made each timestep
	thrust::device_vector<unsigned> delinksThreadMade;//max de-links each thread has made each timestep

	//holds sum of forces for node on given time_step
	thrust::device_vector<double> sumForcesOnNode;

	thrust::device_vector<double> discretizedEdgeStrain; //counts strain of edge
	thrust::device_vector<double> addedEdgeStrain;

	thrust::device_vector<double> discretizedEdgeAlignment;

	thrust::device_vector<bool> isNodeFixed;
	thrust::device_vector<bool> isNodeInPltVol;


	// X,Y,Z, location, velocity and force of all
	thrust::device_vector<double> nodeLocX;
	thrust::device_vector<double> nodeLocY;
	thrust::device_vector<double> nodeLocZ;

	thrust::device_vector<double> nodeVelocity;

//holds forces to advance position and velocity
	thrust::device_vector<double> nodeForceX;
	thrust::device_vector<double> nodeForceY;
	thrust::device_vector<double> nodeForceZ;


};

//Data Structure for platelet location. velocity and force
struct PltInfoVecs {
	//holds sum of forces for platelet on given time_step
	// thrust::device_vector<unsigned> deviceEdgeLeft;
	// thrust::device_vector<unsigned> deviceEdgeRight;
	//
	// thrust::host_vector<unsigned> idEdgesMadeHost;
	//
	// thrust::device_vector<unsigned> idEdgesMadeTemp;

	thrust::device_vector<double> sumForcesOnPlt;

	// thrust::device_vector<double> discretizedEdgeStrain; //counts strain of edge
	// thrust::device_vector<double> addedEdgeStrain;
	//
	// thrust::device_vector<double> discretizedEdgeAlignment;

	thrust::device_vector<bool> isPltFixed;


	// X,Y,Z, location, velocity and force of all
	thrust::device_vector<double> pltLocX;
	thrust::device_vector<double> pltLocY;
	thrust::device_vector<double> pltLocZ;

	thrust::device_vector<double> pltVelocity;

//holds forces to advance position and velocity
	thrust::device_vector<double> pltForceX;
	thrust::device_vector<double> pltForceY;
	thrust::device_vector<double> pltForceZ;

	unsigned numConnections=0;
	thrust::device_vector<unsigned> pltImagingConnection;//used for imaging
	thrust::device_vector<unsigned> nodeImagingConnection;//used for imaging

	thrust::device_vector<unsigned> nodeUnreducedId;
	thrust::device_vector<double> nodeUnreducedForceX;
	thrust::device_vector<double> nodeUnreducedForceY;
	thrust::device_vector<double> nodeUnreducedForceZ;

	thrust::device_vector<unsigned> nodeReducedId;
	thrust::device_vector<double> nodeReducedForceX;
	thrust::device_vector<double> nodeReducedForceY;
	thrust::device_vector<double> nodeReducedForceZ;

	thrust::device_vector<unsigned> tndrlNodeId;
	thrust::device_vector<unsigned> tndrlNodeType;

//


};

//struct used for linking of nodes in network
struct AuxVecs {

	//bucket set to interact with network-network
	thrust::device_vector<unsigned> id_bucket_net_intc; //bucket id
	thrust::device_vector<unsigned> id_value_net_intc; //node id
	thrust::device_vector<unsigned> id_bucket_expanded_net_intc;
	thrust::device_vector<unsigned> id_value_expanded_net_intc;

	thrust::device_vector<unsigned> keyBegin_net_intc;
	thrust::device_vector<unsigned> keyEnd_net_intc;

	unsigned endIndexBucketKeys_net_intc;
	

	//////////////////////////////
	//Same Grid used for both network->plt and plt->network. 
	//////////////////////////////
	//bucket set to interact with network->plt
	thrust::device_vector<unsigned> id_bucket_plt_intc; //bucket id
	thrust::device_vector<unsigned> id_value_plt_intc; //node id
	
	//Do I need the following???
	thrust::device_vector<unsigned> id_bucket_expanded_plt_intc;
	thrust::device_vector<unsigned> id_value_expanded_plt_intc;

	thrust::device_vector<unsigned> keyBegin_plt_intc;
	thrust::device_vector<unsigned> keyEnd_plt_intc;

	unsigned endIndexBucketKeys_plt_intc;

  	//platelets
	  
	//bucket set to interact with plt->network
    thrust::device_vector<unsigned> idPlt_bucket; //bucket id
    thrust::device_vector<unsigned> idPlt_value; //node id
  	thrust::device_vector<unsigned> idPlt_bucket_expanded;
  	thrust::device_vector<unsigned> idPlt_value_expanded;

  	thrust::device_vector<unsigned> keyPltBegin;
  	thrust::device_vector<unsigned> keyPltEnd;

  	unsigned endIndexBucketPltKeys_plt_intc;

};



struct DomainParams {
	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;
  	double pltminX;
	double pltmaxX;
	double pltminY;
	double pltmaxY;
	double pltminZ;
	double pltmaxZ;
	double originMinX;
	double originMaxX;
	double originMinY;
	double originMaxY;
	double originMinZ;
	double originMaxZ;

	double gridSpacing_plt_intc;
	double gridSpacing_net_intc;

	unsigned XBucketCount_net_intc;
	unsigned YBucketCount_net_intc;
	unsigned ZBucketCount_net_intc;
	unsigned totalBucketCount_net_intc = 0;

	unsigned XBucketCount_plt_intc;
	unsigned YBucketCount_plt_intc;
	unsigned ZBucketCount_plt_intc;
	unsigned totalBucketCount_plt_intc = 0;//initialized as zero to set system.
};


//Data for edge node id's
struct WLCInfoVecs {
	double percentOriginalEdgesUnderStrain1 = 0.0;
	double percentOriginalEdgesExtended = 0.0;
	double percentOriginalEdgesCompressed = 0.0;
	double percentAddedEdgesUnderStrain1 = 0.0;
	double percentAddedEdgesExtended = 0.0;
	double percentAddedEdgesCompressed = 0.0;
	double averageStrainAddedEdges = 0.0;
	double averageStrainOriginalEdges = 0.0;

    thrust::device_vector<unsigned> numOriginalNeighborsNodeVector;//holds how many original edges a node was connected to.
	thrust::device_vector<unsigned> globalNeighbors;
	thrust::device_vector<unsigned> currentNodeEdgeCountVector;
	thrust::device_vector<unsigned> currentBindCountPerOriginalEdgeVector;

	thrust::device_vector<double> lengthZero;

};

struct TorsionInfoVecs {
	unsigned factor = 3;
	unsigned currentSpringCount = 0;

	thrust::device_vector<unsigned> leftIndex;
	thrust::device_vector<unsigned> centerIndex;
	thrust::device_vector<unsigned> rightIndex;
	thrust::device_vector<double> angleZero;
	thrust::device_vector<double> forceX;//triple leftIndexLength
	thrust::device_vector<double> forceY;
	thrust::device_vector<double> forceZ;
	thrust::device_vector<double> tempForceX;//triple leftIndexLength
	thrust::device_vector<double> tempForceY;
	thrust::device_vector<double> tempForceZ;

	thrust::device_vector<unsigned> tempTorIndices;
	thrust::device_vector<unsigned> reducedIds;

};


struct GeneralParams{
	//general computation
	bool runSim = true; //default true to begin sim. Sim ends when runSim == false


	double lagTime = 1.0;//set in main.cpp usin

	unsigned maxNeighborCount = 100;
	unsigned maxNodeCount;//after discretize
	unsigned originNodeCount;//pre discretize
  	
	  //platelets
  	unsigned maxPltCount;//after discretize
	unsigned originPltCount;//pre discretize

	unsigned originLinkCount;//constant unsubdivided count of edges
	unsigned originEdgeCount = 0; //total links set at beginning. Constants
	unsigned currentEdgeCount = 0;//total non constant if links are made

	unsigned totalTorsionCount;//total bending springs
	unsigned subNodeCount = 0;//maximal subnode division for longest edge

	double kB, CLM, temperature, torsionStiffness, viscousDamp_Fibrin, viscousDamp_Plt;
	double nodeMass = 1;
	double persistenceLengthMon;

	double fiberDiameter = 0.1;

  //platelet parameters
  //look in builder for default parameter settings.
	double pltForce;
	unsigned plt_other_intrct;//used for unreduced force size. Set large so able to push out many nodes
	unsigned plt_tndrl_intrct;//used for tendril arm counts
	double pltR;
	double pltRForce;
	double pltRAdhesion;
	double pltMass;
	double pltDensity;
	bool pltfrcfld = false;
	bool plttndrl = true; 
	bool pltrelease = true;
	bool plthandhand = true;
	bool pltonplt = true;
	bool agg_release = true;
	unsigned maxIdCountFlag;//flag variable

	//parameters for advancing timestep and determining equilibrium
	double df, dtMax, dtTemp, epsilon, maxForce;
	double currentTime = 0.0;

	//total equilibrium iters and linking determiner
	unsigned iterationCounter = 0;
	bool linking = true;
	unsigned maxLinksPerIteration = 20;


};


class Storage;

class System {
public:
	DomainParams domainParams;
	NodeInfoVecs nodeInfoVecs;
  	PltInfoVecs pltInfoVecs;
	AuxVecs auxVecs;
	WLCInfoVecs wlcInfoVecs;
	TorsionInfoVecs torsionInfoVecs;
	GeneralParams generalParams;

	std::shared_ptr<Storage> storage;

public:

	System();

	void assignStorage(std::shared_ptr<Storage> );

	void initializeSystem(
		thrust::host_vector<bool>& _hostIsNodeFixed,
		thrust::host_vector<double>& _hostPosX,
		thrust::host_vector<double>& _hostPosY,
		thrust::host_vector<double>& _hostPosZ,
		thrust::host_vector<unsigned>& hostWLCEdgeLeft,
		thrust::host_vector<unsigned>& hostWLCEdgeRight,
		thrust::host_vector<double>& hostWLCLenZero,
		thrust::host_vector<unsigned>& hostWLCEdgeSubLeft,
		thrust::host_vector<unsigned>& hostWLCEdgeSubRight,
		thrust::host_vector<double>& hostWLCLenSubZero,
		thrust::host_vector<unsigned>& _hostTorsionIndexLeft,
		thrust::host_vector<unsigned>& _hostTorsionIndexCenter,
		thrust::host_vector<unsigned>& _hostTorsionIndexRight,
		thrust::host_vector<double>& _hostTorsionAngleZero,
    //platelets
    thrust::host_vector<bool>& _hostIsPltFixed,
		thrust::host_vector<double>& _hostPltPosX,
		thrust::host_vector<double>& _hostPltPosY,
		thrust::host_vector<double>& _hostPltPosZ);


	//use from cpu side to begin solving system.
	void solveSystem();

	void solveForces();

	void setBucketScheme();

	void determineBounds();//used for strain.

	void setNodeVecs(
		thrust::host_vector<bool>& hostIsNodeFixed,
		thrust::host_vector<double>& hostPosX,
		thrust::host_vector<double>& hostPosY,
		thrust::host_vector<double>& hostPosZ);

  	void setPltVecs(
  	  thrust::host_vector<bool>& hostIsPltFixed,
  	  thrust::host_vector<double>& hostPltPosX,
  	  thrust::host_vector<double>& hostPltPosY,
  	  thrust::host_vector<double>& hostPltPosZ);

	void setTorsionVecs(
		thrust::host_vector<unsigned>& hostTorsionIndexLeft,
		thrust::host_vector<unsigned>& hostTorsionIndexCenter,
		thrust::host_vector<unsigned>& hostTorsionIndexRight,
		thrust::host_vector<double>& hostTorsionAngleZero);

	void setWLCVecs(
		thrust::host_vector<unsigned>& hostWLCSubEdgeLeft,
		thrust::host_vector<unsigned>& hostWLCSubEdgeRight,
		thrust::host_vector<double>& hostWLCSubLenZero );

};


#endif
