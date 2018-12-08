#ifndef PLTTNDRLDEVICE_H_
#define PLTTNDRLDEVICE_H_

#include <vector>
#include "SystemStructures.h"

    //Tndrl-like force
void PltTndrlOnDevice(
  NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
  PltInfoVecs& pltInfoVecs,
  AuxVecs& auxVecs);

struct PltTndrlonNodeForceFunctor : public thrust::unary_function< U2CVec6, CVec3>  {
  	unsigned plt_tndrl_intrct;
  	double pltRForce;
  	double pltForce;
  	double pltR;

  	unsigned maxPltCount;
  	double fiberDiameter;
  	unsigned maxNodeCount;
  	unsigned maxIdCountFlag;
  	unsigned maxNeighborCount;
	bool pltrelease;
	bool plthandhand;

  	double* nodeLocXAddr;
	double* nodeLocYAddr;
	double* nodeLocZAddr;
	double* nodeUForceXAddr;
	double* nodeUForceYAddr;
	double* nodeUForceZAddr;

  	unsigned* nodeUId;
  	unsigned* pltUId;

	unsigned* id_value_expanded;
	unsigned* keyBegin;
	unsigned* keyEnd;

  	unsigned* tndrlNodeId;
  	unsigned* tndrlNodeType;
  	unsigned* glblNghbrsId;

  	double* pltLocXAddr;
	double* pltLocYAddr;
	double* pltLocZAddr;


   __host__ __device__
   //
       PltTndrlonNodeForceFunctor(
            unsigned& _plt_tndrl_intrct,
            double& _pltRForce,
            double& _pltForce,
            double& _pltR,

            unsigned& _maxPltCount,
            double& _fiberDiameter,
            unsigned& _maxNodeCount,
            unsigned& _maxIdCountFlag,
            unsigned& _maxNeighborCount,
			bool& _pltrelease,
			bool& _plthandhand,

            double* _nodeLocXAddr,
            double* _nodeLocYAddr,
            double* _nodeLocZAddr,
            double* _nodeUForceXAddr,
            double* _nodeUForceYAddr,
            double* _nodeUForceZAddr,

      		unsigned* _nodeUId,
            unsigned* _pltUId,

      		unsigned* _id_value_expanded,
      		unsigned* _keyBegin,
      		unsigned* _keyEnd,

            unsigned* _tndrlNodeId,
            unsigned* _tndrlNodeType,
            unsigned* _glblNghbrsId,

            double* _pltLocXAddr,
            double* _pltLocYAddr,
            double* _pltLocZAddr) :

    plt_tndrl_intrct(_plt_tndrl_intrct),
    pltRForce(_pltRForce),
    pltForce(_pltForce),
    pltR(_pltR),

    maxPltCount(_maxPltCount),
    fiberDiameter(_fiberDiameter),
    maxNodeCount(_maxNodeCount),
    maxIdCountFlag(_maxIdCountFlag),
    maxNeighborCount(_maxNeighborCount),
	pltrelease(_pltrelease),
	plthandhand(_plthandhand),

    nodeLocXAddr(_nodeLocXAddr),
	nodeLocYAddr(_nodeLocYAddr),
	nodeLocZAddr(_nodeLocZAddr),
	nodeUForceXAddr(_nodeUForceXAddr),
	nodeUForceYAddr(_nodeUForceYAddr),
	nodeUForceZAddr(_nodeUForceZAddr),

    nodeUId(_nodeUId),
    pltUId(_pltUId),

	id_value_expanded(_id_value_expanded),
	keyBegin(_keyBegin),
	keyEnd(_keyEnd),

    tndrlNodeId(_tndrlNodeId),
    tndrlNodeType(_tndrlNodeType),
    glblNghbrsId(_glblNghbrsId),

    pltLocXAddr(_pltLocXAddr),
	pltLocYAddr(_pltLocYAddr),
	pltLocZAddr(_pltLocZAddr){}


   __device__
 		CVec3 operator()(const U2CVec6 &u2d6) {

        unsigned pltId = thrust::get<0>(u2d6);
        unsigned bucketId = thrust::get<1>(u2d6);

        //beginning and end of attempted interaction network nodes.
		unsigned beginIndex = keyBegin[bucketId];
		unsigned endIndex = keyEnd[bucketId];


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

    	//double vecN_PX = 0.0;
    	//double vecN_PY = 0.0;
    	//double vecN_PZ = 0.0;
    	//double dist = 0.0;

        //pulling
        //Loop through the number of available tendrils
        for(unsigned interactionCounter = 0; interactionCounter < plt_tndrl_intrct; interactionCounter++) {

            unsigned pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
            unsigned pullNode_type = tndrlNodeType[storageLocation + interactionCounter];

            if (( pullNode_id != maxIdCountFlag) &&
               ( pullNode_type == 0)) {
			   //then we have a node connected to the plt.
                //TYPE 0 is network

                //Calculate distance from plt to node.
               // unsigned pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
                //Get position of node
                double vecN_PX = pltLocX - nodeLocXAddr[pullNode_id];
                double vecN_PY = pltLocY - nodeLocYAddr[pullNode_id];
                double vecN_PZ = pltLocZ - nodeLocZAddr[pullNode_id];
                double dist = sqrt(
                    (vecN_PX) * (vecN_PX) +
                    (vecN_PY) * (vecN_PY) +
                    (vecN_PZ) * (vecN_PZ));

                //check if the node is not in pulling range anymore.
                if ((dist >= pltRForce) || (dist <= (pltR + fiberDiameter / 2.0 ) ) ){

                  	//then node is out of range, so we empty tendril if pltrelease is true

					if (pltrelease==true){
						tndrlNodeId[storageLocation + interactionCounter] = maxIdCountFlag;//reset
					}
                  	//try to find a new node to pull within connections of previous node if plthandhand is true
					if (plthandhand==true){
						unsigned startIndex = maxNeighborCount * pullNode_id;
						unsigned endIndex = startIndex + maxNeighborCount;

						for (unsigned nbr_loc = startIndex; nbr_loc < endIndex; nbr_loc++){

								unsigned newpullNode_id = glblNghbrsId[ nbr_loc ];
								//check tentative node is not already connected
								for (unsigned checkId = 0; checkId < plt_tndrl_intrct; checkId++) {
									if (newpullNode_id != tndrlNodeId[storageLocation + checkId]) {
									  break;
									}
								}
							//check neighbor not empty
							if (newpullNode_id < maxNodeCount){//maxNodeCount is default neighbor value.
								 vecN_PX = pltLocX - nodeLocXAddr[newpullNode_id];
								 vecN_PY = pltLocY - nodeLocYAddr[newpullNode_id];
								 vecN_PZ = pltLocZ - nodeLocZAddr[newpullNode_id];
								//Calculate distance from plt to node.
								 dist = sqrt(
									(vecN_PX) * (vecN_PX) +
									(vecN_PY) * (vecN_PY) +
									(vecN_PZ) * (vecN_PZ));

								//check if new node is in interaction range and fill tendril with new node than break neighbors loop
								if ((dist < pltRForce) && (dist > (pltR + fiberDiameter / 2.0) ) ) {//pull this node
									tndrlNodeId[storageLocation + interactionCounter] = newpullNode_id;//bucketNbrsExp[i];
									tndrlNodeType[storageLocation + interactionCounter] = 0;//assign type
									break;
								}
							}
						}
					}
                }
            }

        	//check if tendril instead still pulls a plt
        	else if ( (pullNode_id != maxIdCountFlag) &&
        	    ( pullNode_type == 1) ) {//note this happens only if plt-plt interaction is on
        	  	//then we have a plt connected to the plt.
				//TYPE 1 is plt

        	    //Calc range
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
        	    if (((dist >= 2.0 * pltRForce) || (dist <= 2.0 * pltR)) && pltrelease==true ){
        	        //then plt is out of range so we disconnect


					tndrlNodeId[storageLocation + interactionCounter] = maxIdCountFlag;
        	    }
        	}

			//after check, re generate choice of node and type.
            pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
            pullNode_type = tndrlNodeType[storageLocation + interactionCounter];

        	// check if tendril still has no node or platelet to pull
        	if (pullNode_id == maxIdCountFlag){
				//then we have nothing to pull.
        	    //try to find a node to pull by searching.

				//ISSUE HERE: we need a random permutation of nodes.
        	    for (unsigned newpull_index = beginIndex; newpull_index < endIndex; newpull_index++){
        	        unsigned newpullNode_id = id_value_expanded[ newpull_index ];
					bool node_is_new = true;

					//check tentative node is not already connected
				    for (unsigned checkId = 0; checkId < plt_tndrl_intrct; checkId++){
        	        	if (newpullNode_id == tndrlNodeId[storageLocation + checkId]){
							node_is_new = false;
        	        		break;
        	        	}
        	      	}
					if ((node_is_new)){

						double vecN_PX = pltLocX - nodeLocXAddr[newpullNode_id];
						double vecN_PY = pltLocY - nodeLocYAddr[newpullNode_id];
						double vecN_PZ = pltLocZ - nodeLocZAddr[newpullNode_id];
						//Calculate distance from plt to node.
						double dist = sqrt(
						(vecN_PX) * (vecN_PX) +
						(vecN_PY) * (vecN_PY) +
						(vecN_PZ) * (vecN_PZ));

						//check if new node is in interaction range and fill tenril with new unique node
						//then break searching loop
						if ((dist < pltRForce) &&
							(dist > (pltR + fiberDiameter / 2.0) ) ) {
							//pull this node

							tndrlNodeId[storageLocation + interactionCounter] = newpullNode_id;//bucketNbrsExp[i];
							tndrlNodeType[storageLocation + interactionCounter] = 0;//assign type
							break;
						}
					}
        	    }

        	}


        	//WARNING: ONLY PULLING NODES
			//UP TO HERE WE HAVE CHECKED FOR NODES
        	//check if tendril has been filled with a node and apply pulling forces. Note if filled direction and distence of forces are already calculated

			//after last check, re generate choice of node and type.
            pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
            pullNode_type = tndrlNodeType[storageLocation + interactionCounter];

			if ( (pullNode_id != maxIdCountFlag)
        	   && ( pullNode_type == 0) ) {
        	    //then we have a post-search node we can pull.
				//Add force to it and the current platelet.

				unsigned pullNode_id = tndrlNodeId[storageLocation + interactionCounter];
				double vecN_PX = pltLocX - nodeLocXAddr[pullNode_id];
        	    double vecN_PY = pltLocY - nodeLocYAddr[pullNode_id];
        	    double vecN_PZ = pltLocZ - nodeLocZAddr[pullNode_id];
        	    //Calculate distance from plt to node.
        	    double dist = sqrt(
        	       (vecN_PX) * (vecN_PX) +
        	       (vecN_PY) * (vecN_PY) +
        	       (vecN_PZ) * (vecN_PZ));
				if ((dist < pltRForce) && (dist > (pltR + fiberDiameter / 2.0))){
					//Determine direction of force based on positions and multiply magnitude force
					double forceNodeX = (vecN_PX / dist) * (pltForce);
					double forceNodeY = (vecN_PY / dist) * (pltForce);
					double forceNodeZ = (vecN_PZ / dist) * (pltForce);

					//count force for self plt.
					sumPltForceX += (-1.0) * forceNodeX;
					sumPltForceY += (-1.0) * forceNodeY;
					sumPltForceZ += (-1.0) * forceNodeZ;

					//store force in temporary vector if a node is pulled. Call reduction later.
					nodeUForceXAddr[storageLocation + interactionCounter] = forceNodeX;
					nodeUForceYAddr[storageLocation + interactionCounter] = forceNodeY;
					nodeUForceZAddr[storageLocation + interactionCounter] = forceNodeZ;
					nodeUId[storageLocation + interactionCounter] = pullNode_id;
					pltUId[storageLocation + interactionCounter] = pltId;
				}

        	}

        }
    //return platelet forces
    return thrust::make_tuple(sumPltForceX, sumPltForceY, sumPltForceZ);

   }
};

#endif