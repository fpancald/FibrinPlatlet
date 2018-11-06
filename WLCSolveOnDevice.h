#ifndef WLCSOLVEONDEVICE_H_
#define WLCSOLVEONDEVICE_H_

#include "SystemStructures.h"


/*
the structure of lengthZero_index is 
0  1  2  3 
4  5  6  7 
8  9  10 11
12 13 14 15 for a 4 node system. 
index/4 = row,
index%4 = col. If you apply force to column node always or row node always then 
each thread will apply opposing forces to springs. 
if you decide to apply force to column instead of rows, you'll need sign change
LengthZero_value is symmetric, so values line up correctly.
*/


void WLCSolveOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams);
void GetStrainParameters(NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,  
	GeneralParams& generalParams,
	DomainParams& domainParams);

	
struct WLCfunctor {
	double* locXAddr;
	double* locYAddr;
	double* locZAddr;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	double Kb; //convert to nN and microns
	double PLengthMon;
	double CLM;
	double Temp;
	unsigned maxNeighborCount;
	unsigned maxNodeCount;

	double* lenZero;
	unsigned* edgeCountVec;
	unsigned* globalNeighbors;
	unsigned* numOriginalNeighborsVec;

	__host__ __device__

		WLCfunctor(
			double* _locXAddr, 
			double* _locYAddr, 
			double* _locZAddr,
			double* _forceXAddr, 
			double* _forceYAddr, 
			double* _forceZAddr, 

			double& _Kb, 
			double& _PLengthMon, 
			double& _CLM, 
			double& _Temp,
			unsigned& _maxNeighborCount,
			unsigned& _maxNodeCount,

			double* _lenZero,
			unsigned* _globalNeighbors,
			unsigned* _edgeCountVec,
			unsigned* _numOriginalNeighborsVec) :

		locXAddr(_locXAddr),
		locYAddr(_locYAddr),
		locZAddr(_locZAddr),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),

		Kb(_Kb), 
		PLengthMon(_PLengthMon), 
		CLM(_CLM), 
		Temp(_Temp),
		maxNeighborCount(_maxNeighborCount),
		maxNodeCount(_maxNodeCount),

		lenZero(_lenZero),
		globalNeighbors(_globalNeighbors),
		edgeCountVec(_edgeCountVec),
		numOriginalNeighborsVec(_numOriginalNeighborsVec) {}
 
	__device__
	void operator()(const Tub& u1b1) {
		//idA represents row.
		unsigned idA = thrust::get<0>(u1b1);
		//unsigned numOriginalNeighbors = numOriginalNeighborsVec[idA];

		bool isFixed = thrust::get<1>(u1b1);
		double sumForceX = 0;
		double sumForceY = 0;
		double sumForceZ = 0;

		if (!isFixed) {
			//only apply force if not fixed. 
			
			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			

			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX
				if (idB < maxNodeCount){
				
					double lengthZero = lenZero[i];
					if (lengthZero > 0) {
						
						double posXA_XB = locXAddr[idB] - locXAddr[idA];
						double posYA_YB = locYAddr[idB] - locYAddr[idA];
						double posZA_ZB = locZAddr[idB] - locZAddr[idA];
		
						double currentLength = sqrt(
							(posXA_XB) * (posXA_XB)+
							(posYA_YB) * (posYA_YB)+
							(posZA_ZB) * (posZA_ZB));
	
						double strain = ((currentLength - lengthZero) / lengthZero);

						double dL_norm = strain / ( CLM);//CLM is unitless since it was already normalized. 
						double magForce = (1100.0*(Kb*Temp) / PLengthMon) * ( 0.25 * pow(1.0 - dL_norm, -2.0) - 0.25 + dL_norm);
				
						double magForceX = (posXA_XB / currentLength) * magForce;
						double magForceY = (posYA_YB / currentLength) * magForce;
						double magForceZ = (posZA_ZB / currentLength) * magForce;
						
						sumForceX += magForceX;
						sumForceY += magForceY;
						sumForceZ += magForceZ;
			
					}
				}
			}
			
			
			if (isfinite(sumForceX))
				forceXAddr[idA] += sumForceX;
			
			if (isfinite(sumForceY))
				forceYAddr[idA] += sumForceY;
			
			if (isfinite(sumForceY))
				forceZAddr[idA] += sumForceZ;	
		}

	}
};


#endif /*WLCSOLVEONDEVICE*/
