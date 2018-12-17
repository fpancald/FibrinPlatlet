#ifndef FUNCTOR_PLT_ARM_PLT_H_
#define FUNCTOR_PLT_ARM_PLT_H_

#include "SystemStructures.h"

struct functor_plt_arm_plt : public thrust::unary_function<U2CVec6, CVec3>  {
  	unsigned plt_tndrl_intrct;
  	double pltRForce;
  	double pltForce;
  	double pltR;

  	unsigned maxPltCount;
  	unsigned maxIdCountFlag;
	bool pltrelease;

	double* pltLocXAddr;
	double* pltLocYAddr;
	double* pltLocZAddr;

	unsigned* idPlt_value_expanded;
	unsigned* keyPltBegin;
	unsigned* keyPltEnd;
	
	unsigned* tndrlNodeId;
	unsigned* tndrlNodeType;


   __host__ __device__
   //
	functor_plt_arm_plt(
		unsigned& _plt_tndrl_intrct,
		double& _pltRForce,
		double& _pltForce,
		double& _pltR,

		unsigned& _maxPltCount,
		unsigned& _maxIdCountFlag,
		bool& _pltrelease,

		double* _pltLocXAddr,
		double* _pltLocYAddr,
		double* _pltLocZAddr,

		unsigned* _idPlt_value_expanded,
		unsigned* _keyPltBegin,
		unsigned* _keyPltEnd,
				
		unsigned* _tndrlNodeId,
		unsigned* _tndrlNodeType) :

	plt_tndrl_intrct(_plt_tndrl_intrct),
	pltRForce(_pltRForce),
	pltForce(_pltForce),
	pltR(_pltR),

	maxPltCount(_maxPltCount),
	maxIdCountFlag(_maxIdCountFlag),
	pltrelease(_pltrelease),

	pltLocXAddr(_pltLocXAddr),
	pltLocYAddr(_pltLocYAddr),
	pltLocZAddr(_pltLocZAddr),

	idPlt_value_expanded(_idPlt_value_expanded),
	keyPltBegin(_keyPltBegin),
	keyPltEnd(_keyPltEnd),

	tndrlNodeId(_tndrlNodeId),
	tndrlNodeType(_tndrlNodeType) {}


   __device__
		CVec3 operator()(const U2CVec6 &u2d6) {

		unsigned pltId = thrust::get<0>(u2d6);
		unsigned bucketId = thrust::get<1>(u2d6);

		//beginning and end of attempted interaction network nodes.
		unsigned beginIndex = keyPltBegin[bucketId];
		unsigned endIndex = keyPltEnd[bucketId];
			
		unsigned storageLocation = pltId * plt_tndrl_intrct;

		double pltLocX = thrust::get<2>(u2d6);
		double pltLocY = thrust::get<3>(u2d6);
		double pltLocZ = thrust::get<4>(u2d6);

        //use for return. 
        double pltCurrentForceX = thrust::get<5>(u2d6);
        double pltCurrentForceY = thrust::get<6>(u2d6);
        double pltCurrentForceZ = thrust::get<7>(u2d6);

        double sumPltForceX = pltCurrentForceX;
        double sumPltForceY = pltCurrentForceY;
        double sumPltForceZ = pltCurrentForceZ;

		double vecN_PX = 0.0;
		double vecN_PY = 0.0;
		double vecN_PZ = 0.0;
		double dist = 0.0;

		//pulling
		//Loop through the number of available tendrils
		for( unsigned interactionCounter = 0; interactionCounter < plt_tndrl_intrct; interactionCounter++ ) {

			//check if tendril pulls a plt
			if ( ((tndrlNodeId[storageLocation + interactionCounter]) != maxIdCountFlag) && 
				((tndrlNodeType[storageLocation + interactionCounter]) == 1)) {
				//note this happens only if plt-plt interaction is on

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
				//check if the plt is not pulled  anymore and pltrelease is true
				if ( ((dist >= 2.0 * pltRForce) || (dist <= 2.0 * pltR)) && pltrelease==true ){
				  	//empty tendril
				  	tndrlNodeId[storageLocation + interactionCounter] = maxIdCountFlag;
				}
			}

			  // check if tendril still has no node or platelet to pull
			if (tndrlNodeId[storageLocation + interactionCounter] == maxIdCountFlag){
			  	//then there was nothing to pull 
				//so try to find a platelet to pu
			  	for (unsigned iter = beginIndex; iter < endIndex; iter++){
					unsigned newpullPlt_id = idPlt_value_expanded[iter];

				   	double vecN_PX = pltLocX - pltLocXAddr[newpullPlt_id];
				   	double vecN_PY = pltLocY - pltLocYAddr[newpullPlt_id];
				   	double vecN_PZ = pltLocZ - pltLocZAddr[newpullPlt_id];
				  	//Calculate distance from plt to node.
				   	double dist = sqrt(
						  (vecN_PX) * (vecN_PX) +
						  (vecN_PY) * (vecN_PY) +
						  (vecN_PZ) * (vecN_PZ));

					//check if new node is in interaction range and fill tenril with new node than break neighbors loop
					if ((dist < 2.0 * pltRForce) && (dist > 2.0 * pltR) ) {//pull this node
						tndrlNodeId[storageLocation + interactionCounter] = newpullPlt_id;//bucketNbrsExp[i];
						tndrlNodeType[storageLocation + interactionCounter] = 1;//assign type
						break;
					}
			  	}
			}

			//check if tendril has been filled with plt and apply pulling forces. Note if filled direction and distence of forces are already calculated
			if ( ((tndrlNodeId[storageLocation + interactionCounter]) != maxIdCountFlag) && 
			  	( (tndrlNodeType[storageLocation + interactionCounter]) == 1) && (dist < 2.0 * pltRForce) && (dist > 2.0 * pltR)){
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
		//after counting force from pulling other platelets, apply to self. 
    	//return platelet forces
    	return thrust::make_tuple(sumPltForceX, sumPltForceY, sumPltForceZ);
	}

};
#endif
