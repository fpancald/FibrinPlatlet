#ifndef PLTONPLTFORCEFIELDDEVICE_H_
#define PLTONPLTFORCEFIELDDEVICE_H_

#include <vector>
#include "SystemStructures.h"

void PltInteractionPltOnDevice(
  	GeneralParams& generalParams,
  	PltInfoVecs& pltInfoVecs,
  	AuxVecs& auxVecs);

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
    		    __attribute__ ((unused)) unsigned beginIndex = keyBegin[bucketId];
    		    __attribute__ ((unused)) unsigned endIndex = keyEnd[bucketId];

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