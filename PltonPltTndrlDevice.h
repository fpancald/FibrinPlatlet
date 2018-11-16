#ifndef PLTONPLTTNDRLDEVICE_H_
#define PLTONPLTTNDRLDEVICE_H_

#include <vector>
#include "SystemStructures.h"

void PltonPltTndrlOnDevice(
	GeneralParams& generalParams,
	PltInfoVecs& pltInfoVecs,
	AuxVecs& auxVecs);

	struct PltonPltTndrlForceFunctor : public thrust::unary_function<U2CVec6, CVec3>  {
  unsigned pltmaxConn;
  double pltRForce;
  double pltForce;
  double pltR;
  unsigned maxPltCount;
  unsigned maxIdCount;

  double* pltLocXAddr;
	double* pltLocYAddr;
	double* pltLocZAddr;

  double* pltForceXAddr;
	double* pltForceYAddr;
	double* pltForceZAddr;

	unsigned* bucketNbrsExp;
	unsigned* keyBegin;
	unsigned* keyEnd;
	
	unsigned* tndrlNodeId;
	unsigned* tndrlNodeType;


   __host__ __device__
   //
	   PltonPltTndrlForceFunctor(
			unsigned& _pltmaxConn,
			double& _pltRForce,
			double& _pltForce,
			double& _pltR,
			unsigned& _maxPltCount,
			unsigned& _maxIdCount,

			double* _pltLocXAddr,
			double* _pltLocYAddr,
			double* _pltLocZAddr,

			double* _pltForceXAddr,
			double* _pltForceYAddr,
			double* _pltForceZAddr,

				unsigned* _bucketNbrsExp,
				unsigned* _keyBegin,
				unsigned* _keyEnd,
				
		unsigned* _tndrlNodeId,
		unsigned* _tndrlNodeType) :

	pltmaxConn(_pltmaxConn),
	pltRForce(_pltRForce),
	pltForce(_pltForce),
	pltR(_pltR),
	maxPltCount(_maxPltCount),
	maxIdCount(_maxIdCount),

	pltLocXAddr(_pltLocXAddr),
		pltLocYAddr(_pltLocYAddr),
		pltLocZAddr(_pltLocZAddr),

	pltForceXAddr(_pltForceXAddr),
	pltForceYAddr(_pltForceYAddr),
	pltForceZAddr(_pltForceZAddr),

		bucketNbrsExp(_bucketNbrsExp),
		keyBegin(_keyBegin),
		keyEnd(_keyEnd),

			tndrlNodeId(_tndrlNodeId),
tndrlNodeType(_tndrlNodeType){}


   __device__
		void operator()(const U2CVec3 &u2d3) {

		unsigned bucketId = thrust::get<0>(u2d3);
		unsigned pltId = thrust::get<1>(u2d3);

		//beginning and end of attempted interaction network nodes.
			__attribute__ ((unused)) unsigned beginIndex = keyBegin[bucketId];
			__attribute__ ((unused)) unsigned endIndex = keyEnd[bucketId];
			
		unsigned storageLocation = pltId * pltmaxConn;

		double pltLocX = thrust::get<2>(u2d3);
		double pltLocY = thrust::get<3>(u2d3);
		double pltLocZ = thrust::get<4>(u2d3);


		double sumPltForceX = 0.0;
		double sumPltForceY = 0.0;
		double sumPltForceZ = 0.0;

					double vecN_PX = 0.0;
				double vecN_PY = 0.0;
				double vecN_PZ = 0.0;
				double dist = 0.0;

			//pulling
			//Loop through the number of available tendrils
			for(unsigned interactionCounter = 0; interactionCounter < pltmaxConn; interactionCounter++) {

			  //check if tendril pulls a plt
			   if (tndrlNodeId[storageLocation + interactionCounter]!=maxIdCount && tndrlNodeType[storageLocation + interactionCounter]==1){//note this happens only if plt-plt interaction is on

				//Calculate distance from plt to node.
				unsigned pullPlt_id = tndrlNodeId[storageLocation + interactionCounter];//bucketNbrsExp[i];
				//Get position of plt
				double vecN_PX = pltLocX - pltLocXAddr[pullPlt_id];
				double vecN_PY = pltLocY - pltLocYAddr[pullPlt_id];
				double vecN_PZ = pltLocZ - pltLocZAddr[pullPlt_id];
				double dist = sqrt(
					(vecN_PX) * (vecN_PX) +
					(vecN_PY) * (vecN_PY) +
					(vecN_PZ) * (vecN_PZ));

				//check if the plt is not pulled  anymore
				if ((dist >= 2*pltRForce) || (dist <= 2*pltR) ){
				  //empty tendril
				  tndrlNodeId[storageLocation + interactionCounter]=maxIdCount;
				}
			  }

			  // check if tendril still has no node or platelet to pull
			if (tndrlNodeId[storageLocation + interactionCounter]==maxIdCount){
			  //try to find a platelet to pull
			  for (unsigned j=0; j<maxPltCount; j++){
				unsigned newpullPlt_id=j;

				   double vecN_PX = pltLocX - pltLocXAddr[newpullPlt_id];
				   double vecN_PY = pltLocY - pltLocYAddr[newpullPlt_id];
				   double vecN_PZ = pltLocZ - pltLocZAddr[newpullPlt_id];
				  //Calculate distance from plt to node.
				   double dist = sqrt(
					  (vecN_PX) * (vecN_PX) +
					  (vecN_PY) * (vecN_PY) +
					  (vecN_PZ) * (vecN_PZ));

				  //check if new node is in interaction range and fill tenril with new node than break neighbors loop
				  if ((dist < 2*pltRForce) && (dist > 2*pltR) ) {//pull this node
					  tndrlNodeId[storageLocation + interactionCounter]=newpullPlt_id;//bucketNbrsExp[i];
											tndrlNodeType[storageLocation + interactionCounter]=1;//assign type
					  break;
				  }

			  }
			}

			  //check if tendril has been filled with plt and apply pulling forces. Note if filled direction and distence of forces are already calculated
			  if (tndrlNodeId[storageLocation + interactionCounter]!=maxIdCount && tndrlNodeType[storageLocation + interactionCounter]==1){
				//Determine direction of force based on positions and multiply magnitude force
				double forceNodeX = (vecN_PX / dist) * (pltForce);
				double forceNodeY = (vecN_PY / dist) * (pltForce);
				double forceNodeZ = (vecN_PZ / dist) * (pltForce);

				//count force for plt.
				sumPltForceX += (-1.0) * forceNodeX;
				sumPltForceY += (-1.0) * forceNodeY;
				sumPltForceZ += (-1.0) * forceNodeZ;

			  }

			}

	}

};

#endif