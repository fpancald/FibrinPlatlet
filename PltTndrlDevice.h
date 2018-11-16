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

struct PltTndrlonNodeForceFunctor : public thrust::unary_function<U2CVec3, CVec3>  {
  unsigned pltmaxConn;
  double pltRForce;
  double pltForce;
  double pltR;
  unsigned maxPltCount;
  double fiberDiameter;
  unsigned maxNodeCount;
  unsigned maxIdCount;
  unsigned maxNeighborCount;

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

  unsigned* tndrlNodeId;
  unsigned* tndrlNodeType;
  unsigned* glblNghbrsId;
  double* pltLocXAddr;
	double* pltLocYAddr;
	double* pltLocZAddr;


   __host__ __device__
   //
       PltTndrlonNodeForceFunctor(
            unsigned& _pltmaxConn,
            double& _pltRForce,
            double& _pltForce,
            double& _pltR,
            unsigned& _maxPltCount,
            double& _fiberDiameter,
            unsigned& _maxNodeCount,
            unsigned& _maxIdCount,
            unsigned& _maxNeighborCount,

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
      			unsigned* _keyEnd,

            unsigned* _tndrlNodeId,
            unsigned* _tndrlNodeType,
            unsigned* _glblNghbrsId,
            double* _pltLocXAddr,
            double* _pltLocYAddr,
            double* _pltLocZAddr) :

    pltmaxConn(_pltmaxConn),
    pltRForce(_pltRForce),
    pltForce(_pltForce),
    pltR(_pltR),
    maxPltCount(_maxPltCount),
    fiberDiameter(_fiberDiameter),
    maxNodeCount(_maxNodeCount),
    maxIdCount(_maxIdCount),
    maxNeighborCount(_maxNeighborCount),

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
		keyEnd(_keyEnd),
    tndrlNodeId(_tndrlNodeId),
    tndrlNodeType(_tndrlNodeType),
    glblNghbrsId(_glblNghbrsId),
    pltLocXAddr(_pltLocXAddr),
		pltLocYAddr(_pltLocYAddr),
		pltLocZAddr(_pltLocZAddr){}


   __device__
 		CVec3 operator()(const U2CVec3 &u2d3) {

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

          //check if current tendril is still connected to a node (i.e. <maxNodecount)
          if (tndrlNodeId[storageLocation + interactionCounter]!=maxIdCount && tndrlNodeType[storageLocation + interactionCounter]==0){

            //Calculate distance from plt to node.
            unsigned pullNode_id = tndrlNodeId[storageLocation + interactionCounter];//bucketNbrsExp[i];
            //Get position of node
            double vecN_PX = pltLocX - nodeLocXAddr[pullNode_id];
            double vecN_PY = pltLocY - nodeLocYAddr[pullNode_id];
            double vecN_PZ = pltLocZ - nodeLocZAddr[pullNode_id];
            double dist = sqrt(
                (vecN_PX) * (vecN_PX) +
                (vecN_PY) * (vecN_PY) +
                (vecN_PZ) * (vecN_PZ));

            //check if the node is not pulled  anymore
            if ((dist >= pltRForce) || (dist <= pltR + fiberDiameter/2) ){
              //empty tendril
              tndrlNodeId[storageLocation + interactionCounter]=maxIdCount;
              //try to find a new node to pull within connections of previous node
              for (unsigned j=0; j<maxNodeCount; j++){
                unsigned newpullNode_id=glblNghbrsId[pullNode_id*maxNodeCount+j];
                //check tentative node is not already connected
                for (unsigned checkId=0; checkId<pltmaxConn; checkId++){
                  if (newpullNode_id!=tndrlNodeId[storageLocation + checkId]){
                    break;
                  }
                }
                //check neighbor not empty
                if (newpullNode_id!=maxNodeCount){
                   vecN_PX = pltLocX - nodeLocXAddr[newpullNode_id];
                   vecN_PY = pltLocY - nodeLocYAddr[newpullNode_id];
                   vecN_PZ = pltLocZ - nodeLocZAddr[newpullNode_id];
                  //Calculate distance from plt to node.
                   dist = sqrt(
                      (vecN_PX) * (vecN_PX) +
                      (vecN_PY) * (vecN_PY) +
                      (vecN_PZ) * (vecN_PZ));

                  //check if new node is in interaction range and fill tendril with new node than break neighbors loop
                  if ((dist < pltRForce) && (dist > pltR + fiberDiameter/2) ) {//pull this node
                      tndrlNodeId[storageLocation + interactionCounter]=newpullNode_id;//bucketNbrsExp[i];
                      tndrlNodeType[storageLocation + interactionCounter]=0;//assign type
                      break;
                  }
                }
              }
            }
          }

          //check if tendril instead still pulls a plt
          else if (tndrlNodeId[storageLocation + interactionCounter]!=maxIdCount && tndrlNodeType[storageLocation + interactionCounter]==1){//note this happens only if plt-plt interaction is on

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

            //try to find a node to pull
            for (unsigned newpullNode_id=0; newpullNode_id<maxNodeCount; newpullNode_id++){
              //check tentative node is not already connected
              for (unsigned checkId=0; checkId<pltmaxConn; checkId++){
                if (newpullNode_id!=tndrlNodeId[storageLocation + checkId]){
                  break;
                }
              }
                 double vecN_PX = pltLocX - nodeLocXAddr[newpullNode_id];
                 double vecN_PY = pltLocY - nodeLocYAddr[newpullNode_id];
                 double vecN_PZ = pltLocZ - nodeLocZAddr[newpullNode_id];
                //Calculate distance from plt to node.
                 double dist = sqrt(
                    (vecN_PX) * (vecN_PX) +
                    (vecN_PY) * (vecN_PY) +
                    (vecN_PZ) * (vecN_PZ));

                //check if new node is in interaction range and fill tenril with new node than break neighbors loop
                if ((dist < pltRForce) && (dist > pltR + fiberDiameter/2) ) {//pull this node
                    tndrlNodeId[storageLocation + interactionCounter]=newpullNode_id;//bucketNbrsExp[i];
                    tndrlNodeType[storageLocation + interactionCounter]=0;//assign type
                    break;
                }
            }

            /*// if still empty go to platelets
            if (tndrlNodeId[storageLocation + interactionCounter]==maxNodeCount+maxPltCount){
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
                      break;
                  }

              }
            }*/
          }

          //check if tendril has been filled with a node and apply pulling forces. Note if filled direction and distence of forces are already calculated
          if (tndrlNodeId[storageLocation + interactionCounter]!=maxIdCount && tndrlNodeType[storageLocation + interactionCounter]==0){
            //node only affects plt position if it is pulled.
			      unsigned pullNode_id=tndrlNodeId[storageLocation + interactionCounter];
            //Determine direction of force based on positions and multiply magnitude force
            double forceNodeX = (vecN_PX / dist) * (pltForce);
            double forceNodeY = (vecN_PY / dist) * (pltForce);
            double forceNodeZ = (vecN_PZ / dist) * (pltForce);

            //count force for plt.
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
/*
        //pushing
        unsigned pushCounter=0;
        //go through all nodes that might be pushed
        for( unsigned i = 0; i < maxNodeCount; i++){
          unsigned pushNode_id=i;
          //
          //Get position of node
          double vecN_PX = pltLocX - nodeLocXAddr[pushNode_id];
          double vecN_PY = pltLocY - nodeLocYAddr[pushNode_id];
          double vecN_PZ = pltLocZ - nodeLocZAddr[pushNode_id];
          //Calculate distance from plt to node.
          double dist = sqrt(
              (vecN_PX) * (vecN_PX) +
              (vecN_PY) * (vecN_PY) +
              (vecN_PZ) * (vecN_PZ));


          //check pushcounter
          if (pushCounter < maxNeighborCount) {
              //repulsion if fiber and platelet overlap
               if (dist < pltR + fiberDiameter/2)  {
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

                  nodeUForceXAddr[storageLocation + pushCounter] = forceNodeX;
                  nodeUForceYAddr[storageLocation + pushCounter] = forceNodeY;
                  nodeUForceZAddr[storageLocation + pushCounter] = forceNodeZ;
                  nodeUId[storageLocation + pushCounter] = pushNode_id;
                  pltUId[storageLocation + pushCounter] = pltId;

                  pushCounter++;
              }
          }
        }

        //go through all plts that might be pushed
        for( unsigned i = 0; i < maxPltCount; i++){
          unsigned pushPlt_id=i;
          //
          //Get position of node
          double vecN_PX = pltLocX - pltLocXAddr[pushPlt_id];
          double vecN_PY = pltLocY - pltLocYAddr[pushPlt_id];
          double vecN_PZ = pltLocZ - pltLocZAddr[pushPlt_id];
          //Calculate distance from plt to node.
          double dist = sqrt(
              (vecN_PX) * (vecN_PX) +
              (vecN_PY) * (vecN_PY) +
              (vecN_PZ) * (vecN_PZ));


          //check pushcounter
          if (pushCounter < maxNeighborCount) {
              //repulsion if fiber and platelet overlap
               if (dist < 2*pltR )  {
                  //node only affects plt position if it is pulled.
                  //Determine direction of force based on positions and multiply magnitude force
                  double forcePltX = -(vecN_PX / dist) * (pltForce);
                  double forcePltY = -(vecN_PY / dist) * (pltForce);
                  double forcePltZ = -(vecN_PZ / dist) * (pltForce);

                  //count force for plt.
                  sumPltForceX += (-1.0) * forcePltX;
                  sumPltForceY += (-1.0) * forcePltY;
                  sumPltForceZ += (-1.0) * forcePltZ;

                  pushCounter++;
              }
          }
        }*/
    //return platelet forces
    return thrust::make_tuple(sumPltForceX, sumPltForceY, sumPltForceZ);

   }
};

#endif