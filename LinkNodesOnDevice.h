#ifndef LINKNODESONDEVICE_H_
#define LINKNODESONDEVICE_H_

#include "SystemStructures.h"

void LinkNodesOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);


struct LinkNodesFunctor {
	double* nodeLocXAddr;
	double* nodeLocYAddr;
	double* nodeLocZAddr;
	unsigned* globalNeighbors;
	unsigned* currentEdgeCountVec;
	double*	lengthZero;

	unsigned* bucketKeys;
	unsigned* bucketValues;
	unsigned* bucketNbrsExp;
	unsigned* keyBegin;
	unsigned* keyEnd;

	double fiberDiameter;
	unsigned maxNbrCount;
	unsigned maxNodeCount;

	unsigned maxLinksPerIteration;
	unsigned* idEdgesMadeTempVec;

	__host__ __device__
		//
		LinkNodesFunctor(
			double* _nodeLocXAddr,
			double* _nodeLocYAddr,
			double* _nodeLocZAddr,
			unsigned* _globalNeighbors,
			unsigned* _currentEdgeCountVec,
			double* _lengthZero,

			unsigned* _bucketKeys,
			unsigned* _bucketValues,
			unsigned* _bucketNbrsExp,
			unsigned* _keyBegin,
			unsigned* _keyEnd,

			double& _fiberDiameter,
			unsigned& _maxNbrCount,
			unsigned& _maxNodeCount,

			unsigned& _maxLinksPerIteration,
			unsigned* _idEdgesMadeTempVec) :

		nodeLocXAddr(_nodeLocXAddr),
		nodeLocYAddr(_nodeLocYAddr),
		nodeLocZAddr(_nodeLocZAddr),
		globalNeighbors(_globalNeighbors),
		currentEdgeCountVec(_currentEdgeCountVec),
		lengthZero(_lengthZero),

		bucketKeys(_bucketKeys),
		bucketValues(_bucketValues),
		bucketNbrsExp(_bucketNbrsExp),
		keyBegin(_keyBegin),
		keyEnd(_keyEnd),

		fiberDiameter(_fiberDiameter),
		maxNbrCount(_maxNbrCount),
		maxNodeCount(_maxNodeCount),

		maxLinksPerIteration(_maxLinksPerIteration),
		idEdgesMadeTempVec(_idEdgesMadeTempVec) {}

	__device__
		//
	void operator() (const Tuu& u2) {
		unsigned bucketId = thrust::get<0>(u2);//bucket containing nodeId
		unsigned nodeId = thrust::get<1>(u2);//node to attempt link from.
		unsigned finalEdgeIdPlaced = 0;

		unsigned possibleEdgeBegin = maxNbrCount * nodeId;
		unsigned possibleEdgeEnd = maxNbrCount * nodeId + maxNbrCount;//test all possible for now
		unsigned numPlacedLinks = 0;

		if (( nodeId < maxNodeCount ) && ( bucketId != ULONG_MAX))  {
			//then nodeId is not a DPD particle.

			//beginning and end of attempted linking end id's in bucketNbrsExp
			unsigned beginIndex = keyBegin[bucketId];
			unsigned endIndex = keyEnd[bucketId];

			for (unsigned iter = beginIndex; iter < endIndex; iter++ ) {
				if (numPlacedLinks >= maxLinksPerIteration) {
					break;
				}

				unsigned candidateId = bucketNbrsExp[iter];//test id
				if (candidateId < maxNodeCount) {
					//then candidateId is not a dpd particle.
					bool candidateIsNew = true;

					//we now attempt to link nodeId and candidateEdge.
					if ((candidateId != nodeId) && (currentEdgeCountVec[nodeId] < maxNbrCount)) {
						double dist = sqrt(
							((nodeLocXAddr[nodeId] - nodeLocXAddr[candidateId]) * (nodeLocXAddr[nodeId] - nodeLocXAddr[candidateId])) +
							((nodeLocYAddr[nodeId] - nodeLocYAddr[candidateId]) * (nodeLocYAddr[nodeId] - nodeLocYAddr[candidateId])) +
							((nodeLocZAddr[nodeId] - nodeLocZAddr[candidateId]) * (nodeLocZAddr[nodeId] - nodeLocZAddr[candidateId])));

						if ((dist < fiberDiameter)) {
							//then we have a possible link if no previous link was placed.

							//unsigned targetPossibleEdgeBegin = maxNbrCount * candidateId;
							//unsigned targetPossibleEdgeEnd = maxNbrCount * candidateId + maxNbrCount;//test all possile

							//make sure placement is possible for nodeID
							//other thread will take care of candidateID
							for (unsigned k = possibleEdgeBegin; k < possibleEdgeEnd; k++) {
								if (candidateId == (globalNeighbors[k])) {
									candidateIsNew = false;
									break;
								}
							}

							if (candidateIsNew == true) {

								//must regenerate for each new candidate
								unsigned currentEdgeCount = currentEdgeCountVec[nodeId];

								if (currentEdgeCount < maxNbrCount) {
									//then the new edge can be placed
									unsigned location;
									//decide placement of nbr only in row of nodeId,
									//another thread will place for candidateId
									for (unsigned place = possibleEdgeBegin; place < possibleEdgeEnd; place++ ) {

										if (place < (maxNbrCount * maxNodeCount)) {
											if (maxNodeCount < globalNeighbors[place]) {
												//then we have found location where we can place a node scine global neighbors default is ULONG_MAX

												location = place;
												break;
											}
										}
									}
									if (location < (maxNbrCount * maxNodeCount) ) {
										globalNeighbors[location] = candidateId; //1
										lengthZero[location ] = dist;//1

										(currentEdgeCountVec[nodeId])+=1;

										//only if a node is placed do we record it.
										//place in matrix format upper diagonal
										if (nodeId < candidateId) {
											finalEdgeIdPlaced = candidateId + maxNodeCount * nodeId;
										}
										else {
											finalEdgeIdPlaced = maxNodeCount * candidateId + nodeId;
										};
										//place edge and increment counter
										unsigned index = maxLinksPerIteration * nodeId + numPlacedLinks;
										idEdgesMadeTempVec[index] = finalEdgeIdPlaced;
										numPlacedLinks += 1;
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

//unused
struct DeLinkCopiesFunctor {

	unsigned* globalNeighbors;
	unsigned* currentEdgeCountVec;
	unsigned maxNbrCount;
	double* lengthZero;
	unsigned maxNodeCount;

	__host__ __device__
	DeLinkCopiesFunctor(

		unsigned* _globalNeighbors,
		double* _lengthZero,
		unsigned* _currentEdgeCountVec,
		unsigned _maxNbrCount,
		unsigned _maxNodeCount) :
		globalNeighbors(_globalNeighbors),
		lengthZero(_lengthZero),
		currentEdgeCountVec(_currentEdgeCountVec),
		maxNbrCount(_maxNbrCount),
		maxNodeCount(_maxNodeCount) {}

	__device__
	void operator() (const unsigned& nodeId) {
		unsigned possibleEdgeBegin = maxNbrCount * nodeId;
		unsigned possibleEdgeEnd = maxNbrCount * nodeId + maxNbrCount;//test all possible for now

		for (unsigned index = possibleEdgeBegin; index < possibleEdgeEnd; index++) {
		//index is the location of nbrId

			unsigned nbrId = globalNeighbors[index];//loop through all possible neighbors
			if (nbrId < maxNodeCount) {
				//if ((currentEdgeCountVec[nbrId] == maxNbrCount)) {
					bool isNodeNbrOfNodeId = false;//assume not nbr, and check if true.

					unsigned targetPossibleEdgeBegin = maxNbrCount * nbrId;
					unsigned targetPossibleEdgeEnd = maxNbrCount * nbrId + maxNbrCount;
					for (unsigned indexNbr = targetPossibleEdgeBegin; indexNbr < targetPossibleEdgeEnd; indexNbr++ ) {
						if ((globalNeighbors[indexNbr]) == nodeId) {
							isNodeNbrOfNodeId = true;//nbrId is connected to nodeId
							break;
						}
					}
					if (!isNodeNbrOfNodeId) {
						//then nbrId is not connected to nodeId, but nodeId is connected to nbrId.
						//we thus remove nbrId from nodeId.
						//these two mess everything up
						globalNeighbors[index] = ULONG_MAX;
						//lengthZero[index] = 0.0;
						(currentEdgeCountVec[nodeId])-=1;//subtract from edgecount if nbr removed.

					}

				//}
			}
		}
	}

};

#endif /*LINKNODESONDEVICE_H_*/
