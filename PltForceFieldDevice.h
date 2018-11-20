#ifndef PLTFORCEFIELDDEVICE_H_
#define PLTFORCEFIELDDEVICE_H_

#include <vector>
#include "SystemStructures.h"

//Force field-like mode
void PltForceFieldOnDevice(
  NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
  PltInfoVecs& pltInfoVecs,
  AuxVecs& auxVecs);

  struct PltonNodeForceFunctor : public thrust::unary_function<U2CVec6, CVec3>  {
    unsigned plt_other_intrct;
    double pltRForce;
    double pltForce;
    double pltR;

    unsigned maxPltCount;
    double fiberDiameter;
    unsigned maxNodeCount;

    double* nodeLocXAddr;
  	double* nodeLocYAddr;
  	double* nodeLocZAddr;
  	double* nodeUForceXAddr;
  	double* nodeUForceYAddr;
  	double* nodeUForceZAddr;

    unsigned* nodeUId;

  	unsigned* id_value_expanded;
  	unsigned* keyBegin;
  	unsigned* keyEnd;


     __host__ __device__
     //
         PltonNodeForceFunctor(
              unsigned& _plt_other_intrct,
              double& _pltRForce,
              double& _pltForce,
              double& _pltR,

              unsigned& _maxPltCount,
              double& _fiberDiameter,
              unsigned& _maxNodeCount,

              double* _nodeLocXAddr,
              double* _nodeLocYAddr,
              double* _nodeLocZAddr,
              double* _nodeUForceXAddr,
              double* _nodeUForceYAddr,
              double* _nodeUForceZAddr,

        			unsigned* _nodeUId,

        			unsigned* _id_value_expanded,
        			unsigned* _keyBegin,
        			unsigned* _keyEnd) :

      plt_other_intrct(_plt_other_intrct),
      pltRForce(_pltRForce),
      pltForce(_pltForce),
      pltR(_pltR),

      maxPltCount(_maxPltCount),
      fiberDiameter(_fiberDiameter),
      maxNodeCount(_maxNodeCount),

      nodeLocXAddr(_nodeLocXAddr),
  		nodeLocYAddr(_nodeLocYAddr),
  		nodeLocZAddr(_nodeLocZAddr),
  		nodeUForceXAddr(_nodeUForceXAddr),
  		nodeUForceYAddr(_nodeUForceYAddr),
  		nodeUForceZAddr(_nodeUForceZAddr),

      nodeUId(_nodeUId),

  		id_value_expanded(_id_value_expanded),
  		keyBegin(_keyBegin),
  		keyEnd(_keyEnd) {}


     __device__
   		CVec3 operator()(const U2CVec6 &u2d6) {

          unsigned pltId = thrust::get<0>(u2d6);
          unsigned bucketId = thrust::get<1>(u2d6);

          //beginning and end of attempted interaction network nodes.
  		    unsigned beginIndex = keyBegin[bucketId];
  		    unsigned endIndex = keyEnd[bucketId];


          unsigned storageLocation = pltId * plt_other_intrct;

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

          //Loop through the number of available neighbors for each plt.
          unsigned interactionCounter = 0;

          //change with bucket counter
          for(unsigned i = beginIndex; i < endIndex; i++) {

            //Choose a neighbor.
            unsigned pullNode_id = id_value_expanded[i];
            //
            //Get position of node
            double vecN_PX = pltLocX - nodeLocXAddr[pullNode_id];
            double vecN_PY = pltLocY - nodeLocYAddr[pullNode_id];
            double vecN_PZ = pltLocZ - nodeLocZAddr[pullNode_id];
            //Calculate distance from plt to node.
            double dist = sqrt(
                (vecN_PX) * (vecN_PX) +
                (vecN_PY) * (vecN_PY) +
                (vecN_PZ) * (vecN_PZ));


            //only pull as many as are arms.
            if (interactionCounter < plt_other_intrct) {
                //attraction if platelet and fiber are within interaction distance but not overlapping
                if ((dist < pltRForce) && (dist > (pltR + fiberDiameter / 2.0 ) ) ) {
                    //node only affects plt position if it is pulled.
                    //Determine direction of force based on positions and multiply magnitude force
                    double forceNodeX = (vecN_PX / dist) * (pltForce);
                    double forceNodeY = (vecN_PY / dist) * (pltForce);
                    double forceNodeZ = (vecN_PZ / dist) * (pltForce);

                    //count force for plt.
                    sumPltForceX += (-1.0) * forceNodeX;
                    sumPltForceY += (-1.0) * forceNodeY;
                    sumPltForceZ += (-1.0) * forceNodeZ;

                    //store force in temporary vector. Call reduction later.

                    nodeUForceXAddr[storageLocation + interactionCounter] = forceNodeX;
                    nodeUForceYAddr[storageLocation + interactionCounter] = forceNodeY;
                    nodeUForceZAddr[storageLocation + interactionCounter] = forceNodeZ;
                    nodeUId[storageLocation + interactionCounter] = pullNode_id;
                    
                    interactionCounter++;

                }
                //repulsion if fiber and platelet overlap
                else if (dist < (pltR + fiberDiameter / 2 ) ) {
                    //node only affects plt position if it is pulled.
                    //Determine direction of force based on positions and multiply magnitude force
                    double forceNodeX = -(vecN_PX / dist) * (pltForce);
                    double forceNodeY = -(vecN_PY / dist) * (pltForce);
                    double forceNodeZ = -(vecN_PZ / dist) * (pltForce);

                    //count force for plt.
                    sumPltForceX += (-1.0) * forceNodeX;
                    sumPltForceY += (-1.0) * forceNodeY;
                    sumPltForceZ += (-1.0) * forceNodeZ;

                    //store force in temporary vector. Call reduction later.

                    nodeUForceXAddr[storageLocation + interactionCounter] = forceNodeX;
                    nodeUForceYAddr[storageLocation + interactionCounter] = forceNodeY;
                    nodeUForceZAddr[storageLocation + interactionCounter] = forceNodeZ;
                    nodeUId[storageLocation + interactionCounter] = pullNode_id;
                    
                    interactionCounter++;
                }
            }


      }
      return thrust::make_tuple(sumPltForceX, sumPltForceY, sumPltForceZ);

     }
  };

  #endif