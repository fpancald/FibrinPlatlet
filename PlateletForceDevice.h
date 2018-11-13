
#ifndef PLATELETFORCEDEVICE_H_
#define PLATELETFORCEDEVICE_H_

#include <vector>
#include "SystemStructures.h"



void PltForceOnDevice(
  NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
  PltInfoVecs& pltInfoVecs,
  AuxVecs& auxVecs);


void PltInteractionOnDevice(
  	GeneralParams& generalParams,
  	PltInfoVecs& pltInfoVecs,
  	AuxVecs& auxVecs);


struct AddPltonNodeForceFunctor {//same as torsion
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	__host__ __device__
	//
		AddPltonNodeForceFunctor(
				double* _forceXAddr,
				double* _forceYAddr,
				double* _forceZAddr) :
			forceXAddr(_forceXAddr),
			forceYAddr(_forceYAddr),
			forceZAddr(_forceZAddr) {}

	__device__
	void operator() (const Tuddd& u1d3) {
			unsigned idToAssign = thrust::get<0>(u1d3);
			if (!isnan(thrust::get<1>(u1d3)) && !isnan(thrust::get<2>(u1d3)) && !isnan(thrust::get<3>(u1d3))) {

			forceXAddr[idToAssign] += thrust::get<1>(u1d3);
			forceYAddr[idToAssign] += thrust::get<2>(u1d3);
			forceZAddr[idToAssign] += thrust::get<3>(u1d3);
			}

	}

};

//go through and add appropriate entries for input
struct PltonNodeForceFunctor : public thrust::unary_function<U2CVec3, CVec3>  {
  unsigned pltmaxConn;
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
  unsigned* pltUId;

	unsigned* bucketNbrsExp;
	unsigned* keyBegin;
	unsigned* keyEnd;


   __host__ __device__
   //
       PltonNodeForceFunctor(
            unsigned& _pltmaxConn,
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
            unsigned* _pltUId,

      			unsigned* _bucketNbrsExp,
      			unsigned* _keyBegin,
      			unsigned* _keyEnd) :

    pltmaxConn(_pltmaxConn),
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
    pltUId(_pltUId),

		bucketNbrsExp(_bucketNbrsExp),
		keyBegin(_keyBegin),
		keyEnd(_keyEnd) {}


   __device__
 		CVec3 operator()(const U2CVec3 &u2d3) {

        unsigned bucketId = thrust::get<0>(u2d3);
        unsigned pltId = thrust::get<1>(u2d3);

        //beginning and end of attempted interaction network nodes.
		    unsigned beginIndex = keyBegin[bucketId];
		    unsigned endIndex = keyEnd[bucketId];


        unsigned storageLocation = pltId * pltmaxConn;

        double pltLocX = thrust::get<2>(u2d3);
        double pltLocY = thrust::get<3>(u2d3);
        double pltLocZ = thrust::get<4>(u2d3);


        double sumPltForceX = 0.0;
        double sumPltForceY = 0.0;
        double sumPltForceZ = 0.0;

        //Loop through the number of available neighbors for each plt.
        unsigned interactionCounter=0;

        for(unsigned i = 0; i < maxNodeCount; i++) {

          //Choose a neighbor.
          unsigned pullNode_id = i;//bucketNbrsExp[i];
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
          if (interactionCounter < pltmaxConn) {
              //attraction if platelet and fiber are within interaction distance but not overlapping
              if ((dist < pltRForce) && (dist > pltR + fiberDiameter/2) ) {
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
                  pltUId[storageLocation + interactionCounter] = pltId;

                  interactionCounter++;

              }
              //repulsion if fiber and platelet overlap
              else if (dist < pltR + fiberDiameter/2)  {
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
                  pltUId[storageLocation + interactionCounter] = pltId;
                  
                  interactionCounter++;
              }
          }


    }
    return thrust::make_tuple(sumPltForceX, sumPltForceY, sumPltForceZ);

   }
};

//go through and add appropriate entries for input
struct PltonPltForceFunctor : public thrust::unary_function<U2CVec6, CVec3>  {
  unsigned pltmaxConn;
  double pltRForce;
  double pltForce;
  double pltR;
  unsigned maxPltCount;

  double* pltLocXAddr;
	double* pltLocYAddr;
	double* pltLocZAddr;
  
  double* pltForceXAddr;
	double* pltForceYAddr;
	double* pltForceZAddr;

	unsigned* bucketNbrsExp;
	unsigned* keyBegin;
	unsigned* keyEnd;


   __host__ __device__
   //
       PltonPltForceFunctor(
            unsigned& _pltmaxConn,
            double& _pltRForce,
            double& _pltForce,
            double& _pltR,
            unsigned& _maxPltCount,

            double* _pltLocXAddr,
            double* _pltLocYAddr,
            double* _pltLocZAddr,

            double* _pltForceXAddr,
            double* _pltForceYAddr,
            double* _pltForceZAddr,

      			unsigned* _bucketNbrsExp,
      			unsigned* _keyBegin,
      			unsigned* _keyEnd) :

    pltmaxConn(_pltmaxConn),
    pltRForce(_pltRForce),
    pltForce(_pltForce),
    pltR(_pltR),
    maxPltCount(_maxPltCount),

    pltLocXAddr(_pltLocXAddr),
		pltLocYAddr(_pltLocYAddr),
		pltLocZAddr(_pltLocZAddr),

    pltForceXAddr(_pltForceXAddr),
    pltForceYAddr(_pltForceYAddr),
    pltForceZAddr(_pltForceZAddr),

		bucketNbrsExp(_bucketNbrsExp),
		keyBegin(_keyBegin),
		keyEnd(_keyEnd) {}


   __device__
 		void operator()(const U2CVec3 &u2d3) {

        unsigned bucketId = thrust::get<0>(u2d3);
        unsigned pltId = thrust::get<1>(u2d3);

        //beginning and end of attempted interaction network nodes.
		    unsigned beginIndex = keyBegin[bucketId];
		    unsigned endIndex = keyEnd[bucketId];

        double pltLocX = thrust::get<2>(u2d3);
        double pltLocY = thrust::get<3>(u2d3);
        double pltLocZ = thrust::get<4>(u2d3);


        double sumPltForceX = 0.0;
        double sumPltForceY = 0.0; 
        double sumPltForceZ = 0.0;

        //Loop through the number of available neighbors for each plt.
        unsigned interactionCounter=0;

//for now no bucket scheme
        for(unsigned i = 0; i < maxPltCount; i++) {

          //Choose a neighbor.
          unsigned pullPlt_id = i;//bucketNbrsExp[i];
          //
          if ( (pullPlt_id != pltId) && (pullPlt_id < maxPltCount) ) {
            //Get position of plt
            double vecN_PX = pltLocX - pltLocXAddr[pullPlt_id];
            double vecN_PY = pltLocY - pltLocYAddr[pullPlt_id];
            double vecN_PZ = pltLocZ - pltLocZAddr[pullPlt_id];
            //Calculate distance from plt to plt.
            double dist = sqrt(
                (vecN_PX) * (vecN_PX) +
                (vecN_PY) * (vecN_PY) +
                (vecN_PZ) * (vecN_PZ));


            //only pull as many as are arms.
            if (interactionCounter < pltmaxConn) {
                //attraction if platelet and fiber are within interaction distance but not overlapping

                  if ((dist < 2.0 * pltRForce) && (dist > 2.0 * pltR) ) {
                      //plt only affects plt position if it is pulled.
                      //Determine direction of force based on positions and multiply magnitude force
                      double forcePltX = (vecN_PX / dist) * (pltForce);
                      double forcePltY = (vecN_PY / dist) * (pltForce);
                      double forcePltZ = (vecN_PZ / dist) * (pltForce);

                      //count force for plt.
                      sumPltForceX += 2.0 * (-1.0) * forcePltX;
                      sumPltForceY += 2.0 * (-1.0) * forcePltY;
                      sumPltForceZ += 2.0 * (-1.0) * forcePltZ;

                  }
                  //repulsion if fiber and platelet overlap
                  else if (dist < 2 * pltR )  {
                      //plt only affects plt position if it is pulled.
                      //Determine direction of force based on positions and multiply magnitude force
                      double forcePltX = - (vecN_PX / dist) * (pltForce);
                      double forcePltY = - (vecN_PY / dist) * (pltForce);
                      double forcePltZ = - (vecN_PZ / dist) * (pltForce);

                      //count force for plt.
                      sumPltForceX += 2.0 * (-1.0) * forcePltX;
                      sumPltForceY += 2.0 * (-1.0) * forcePltY;
                      sumPltForceZ += 2.0 * (-1.0) * forcePltZ;

                  }
            }
        }

      }
      pltForceXAddr[pltId] += sumPltForceX;
      pltForceYAddr[pltId] += sumPltForceY;
      pltForceZAddr[pltId] += sumPltForceZ;
    }
};
#endif