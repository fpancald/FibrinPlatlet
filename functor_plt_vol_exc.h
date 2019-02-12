#ifndef FUNCTOR_PLT_VOL_EXC_H_
#define FUNCTOR_PLT_VOL_EXC_H_

#include <vector>
#include "SystemStructures.h"

//This functor performs a volume exclusion for each platelet against 
//fiber nodes and other platelets
struct functor_plt_vol_exc : public thrust::unary_function< U2CVec6, CVec3 >  {
    unsigned plt_other_intrct;
    double pltRForce;
	double pltRAdhesion;
    double pltForce;
    double pltR;

    unsigned maxPltCount;
    double fiberDiameter;
    unsigned maxNodeCount;
    unsigned maxNeighborCount;

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
    
  	unsigned* idPlt_value_expanded;
  	unsigned* keyPltBegin;
  	unsigned* keyPltEnd;

    double* pltLocXAddr;
  	double* pltLocYAddr;
  	double* pltLocZAddr;


     __host__ __device__
     //
         functor_plt_vol_exc(
              unsigned& _plt_other_intrct,
              double& _pltRForce,
			  double& _pltRAdhesion,
              double& _pltForce,
              double& _pltR,

              unsigned& _maxPltCount,
              double& _fiberDiameter,
              unsigned& _maxNodeCount,
              unsigned& _maxNeighborCount,

              double* _nodeLocXAddr,
              double* _nodeLocYAddr,
              double* _nodeLocZAddr,
              double* _nodeUForceXAddr,
              double* _nodeUForceYAddr,
              double* _nodeUForceZAddr,
              
        			unsigned* _nodeUId,

        			unsigned* _id_value_expanded,
        			unsigned* _keyBegin,
        			unsigned* _keyEnd,
              
        			unsigned* _idPlt_value_expanded,
        			unsigned* _keyPltBegin,
        			unsigned* _keyPltEnd,

              double* _pltLocXAddr,
              double* _pltLocYAddr,
              double* _pltLocZAddr) :

      plt_other_intrct(_plt_other_intrct),
      pltRForce(_pltRForce),
	  pltRAdhesion(_pltRAdhesion),
      pltForce(_pltForce),
      pltR(_pltR),

      maxPltCount(_maxPltCount),
      fiberDiameter(_fiberDiameter),
      maxNodeCount(_maxNodeCount),
      maxNeighborCount(_maxNeighborCount),

      nodeLocXAddr(_nodeLocXAddr),
  		nodeLocYAddr(_nodeLocYAddr),
  		nodeLocZAddr(_nodeLocZAddr),
  		nodeUForceXAddr(_nodeUForceXAddr),
  		nodeUForceYAddr(_nodeUForceYAddr),
  		nodeUForceZAddr(_nodeUForceZAddr),

      nodeUId(_nodeUId),

  		id_value_expanded(_id_value_expanded),//these are for network
  		keyBegin(_keyBegin),
  		keyEnd(_keyEnd),

      idPlt_value_expanded(_idPlt_value_expanded),//these are for plt
  		keyPltBegin(_keyPltBegin),
  		keyPltEnd(_keyPltEnd),
      
      pltLocXAddr(_pltLocXAddr),
  		pltLocYAddr(_pltLocYAddr),
  		pltLocZAddr(_pltLocZAddr) {}


     __device__
   		CVec3 operator()(const U2CVec6 &u2d6) {

        unsigned pltId = thrust::get<0>(u2d6);
        unsigned bucketId = thrust::get<1>(u2d6);

        //beginning and end of attempted interaction network nodes.
  		unsigned beginIndexNode = keyBegin[bucketId];
  		unsigned endIndexNode = keyEnd[bucketId];

        unsigned beginIndexPlt = keyPltBegin[bucketId];
  		unsigned endIndexPlt = keyPltEnd[bucketId];


        unsigned storageLocation = pltId * plt_other_intrct;

        double pltLocX = thrust::get<2>(u2d6);
        double pltLocY = thrust::get<3>(u2d6);
        double pltLocZ = thrust::get<4>(u2d6);
          
        //use for return. 
        double sumPltForceX = thrust::get<5>(u2d6);;
        double sumPltForceY = thrust::get<6>(u2d6);;
        double sumPltForceZ = thrust::get<7>(u2d6);;
        
        //pushing
        unsigned pushCounter = 0;
		
		// lennard-jones parameters. set so that when platelets start touching each other they first attract (for a short range and then repell)
		double sigma=2.0 * pltR*pltRAdhesion/pow(2,1/6);
		double eps=pltForce;

        //go through all nodes that might be pushed
        for( unsigned id_count = beginIndexNode; id_count < endIndexNode; id_count++){
            unsigned pushNode_id = id_value_expanded[id_count];
              
            if (pushCounter < plt_other_intrct) {//must be same as other counter.
                  
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
                //repulsion if fiber and platelet overlap more than adhesion fractional radius
                if (dist < (pltRAdhesion*pltR + fiberDiameter / 2.0) )  {
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
                    pushCounter++;
                }
				//adhesion if between fractional adhesion radius and radius
				else if ((dist > (pltRAdhesion*pltR + fiberDiameter / 2.0)) && (dist < (pltR + fiberDiameter / 2.0)) )  {
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
                    nodeUForceXAddr[storageLocation + pushCounter] = forceNodeX;
                    nodeUForceYAddr[storageLocation + pushCounter] = forceNodeY;
                    nodeUForceZAddr[storageLocation + pushCounter] = forceNodeZ;
                    nodeUId[storageLocation + pushCounter] = pushNode_id;
                    pushCounter++;
                }
            }
        }

          //IN THIS SECTION THERE IS NO WRITING TO VECTORS. NO NEED FOR INCREMENT
          //FORCE IS ONLY APPLIED TO SELF. 
          //go through all plts that might be pushed
          //only apply force to self. 
          for( unsigned i = beginIndexPlt; i < endIndexPlt; i++){
              unsigned pushPlt_id = idPlt_value_expanded[i];
              
              
              if (pushPlt_id != pltId) {
                  //then pushPlt id can be pushed since it is not the self
                  double vecN_PX = pltLocX - pltLocXAddr[pushPlt_id];
                  double vecN_PY = pltLocY - pltLocYAddr[pushPlt_id];
                  double vecN_PZ = pltLocZ - pltLocZAddr[pushPlt_id];
                  //Calculate distance from plt to node.
                  double dist = sqrt(
                      (vecN_PX) * (vecN_PX) +
                      (vecN_PY) * (vecN_PY) +
                      (vecN_PZ) * (vecN_PZ));
      
      
                  //check pushcounter
                  //repulsion if fiber and platelet overlap
                  if (dist < (2.0 * pltR ) )  {
                      //node only affects plt position if it is pulled.
                      //Determine direction of force based on positions and multiply magnitude force
					  
					  //Lennard-Jones force coeff
					  double LJF=48*eps*(pow(sigma,12)/pow(dist,13)-pow(sigma,6)/pow(dist,7));
					  
                      double forcePltX = -(vecN_PX / dist) * (pltForce) *LJF;// (dist * dist + dist - pltR) / (dist * dist);//(1.0 + 2.0 * pltR - dist);
                      double forcePltY = -(vecN_PY / dist) * (pltForce) *LJF;// (dist * dist + dist - pltR) / (dist * dist);//(1.0 + 2.0 * pltR - dist);
                      double forcePltZ = -(vecN_PZ / dist) * (pltForce) *LJF;// (dist * dist + dist - pltR) / (dist * dist);//(1.0 + 2.0 * pltR - dist);
                      //count force for plt.
                      sumPltForceX += (-1.0) * forcePltX;
                      sumPltForceY += (-1.0) * forcePltY;
                      sumPltForceZ += (-1.0) * forcePltZ;
                  }
              }
          }
      //return platelet forces
      return thrust::make_tuple(sumPltForceX, sumPltForceY, sumPltForceZ);

     }
};
#endif