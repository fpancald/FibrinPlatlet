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

	unsigned* id_value_expanded;
	unsigned* keyBegin;
	unsigned* keyEnd;

	double fiberDiameter;
	unsigned maxNbrCount;
	unsigned maxNodeCount;

	unsigned maxLinksPerIteration;
	unsigned* idTempVecLeft;
	unsigned* idTempVecRight;

	__host__ __device__
		//
		LinkNodesFunctor(
			double* _nodeLocXAddr,
			double* _nodeLocYAddr,
			double* _nodeLocZAddr,
			unsigned* _globalNeighbors,
			unsigned* _currentEdgeCountVec,
			double* _lengthZero,

			unsigned* _id_value_expanded,
			unsigned* _keyBegin,
			unsigned* _keyEnd,

			double& _fiberDiameter,
			unsigned& _maxNbrCount,
			unsigned& _maxNodeCount,

			unsigned& _maxLinksPerIteration,
			unsigned* _idTempVecLeft,
			unsigned* _idTempVecRight) :

		nodeLocXAddr(_nodeLocXAddr),
		nodeLocYAddr(_nodeLocYAddr),
		nodeLocZAddr(_nodeLocZAddr),
		globalNeighbors(_globalNeighbors),
		currentEdgeCountVec(_currentEdgeCountVec),
		lengthZero(_lengthZero),

		id_value_expanded(_id_value_expanded),
		keyBegin(_keyBegin),
		keyEnd(_keyEnd),

		fiberDiameter(_fiberDiameter),
		maxNbrCount(_maxNbrCount),
		maxNodeCount(_maxNodeCount),

		maxLinksPerIteration(_maxLinksPerIteration),
		idTempVecLeft(_idTempVecLeft),
		idTempVecRight(_idTempVecRight) {}

	__device__
		
	unsigned operator() (const Tuu& u2) {
		unsigned nodeId = thrust::get<0>(u2);//node to attempt link from.
		unsigned bucketId = thrust::get<1>(u2);//bucket containing nodeId
		
		unsigned final_id_left = 0;
		unsigned final_id_right = 0;

		unsigned possibleEdgeBegin = maxNbrCount * nodeId;
		unsigned possibleEdgeEnd = maxNbrCount * nodeId + maxNbrCount;//test all possible for now
		unsigned numPlacedLinks = 0;

		unsigned save_index = maxLinksPerIteration * nodeId;//changes after adding an edge
		unsigned save_index_end = maxLinksPerIteration * nodeId + maxLinksPerIteration;


		if ( nodeId < maxNodeCount )  {
			//tsafety

			//beginning and end of attempted linking end id's in id_value_expanded
			__attribute__ ((unused)) unsigned beginIndex = keyBegin[bucketId];
			__attribute__ ((unused)) unsigned endIndex = keyEnd[bucketId];

			for (unsigned iter = 0; iter < maxNodeCount; iter++ ) {

				unsigned candidateId = iter;//id_value_expanded[iter];//test id

				if ( (candidateId < maxNodeCount) && (save_index < save_index_end) ) {
					//then candidateId is not a dpd particle.
					bool candidateIsNew = true;

					//we now attempt to link nodeId and candidateEdge.
					if ((candidateId != nodeId) ){//&& (currentEdgeCountVec[nodeId] < maxNbrCount)) {
						double dist = sqrt(
							((nodeLocXAddr[nodeId] - nodeLocXAddr[candidateId]) * (nodeLocXAddr[nodeId] - nodeLocXAddr[candidateId])) +
							((nodeLocYAddr[nodeId] - nodeLocYAddr[candidateId]) * (nodeLocYAddr[nodeId] - nodeLocYAddr[candidateId])) +
							((nodeLocZAddr[nodeId] - nodeLocZAddr[candidateId]) * (nodeLocZAddr[nodeId] - nodeLocZAddr[candidateId])));

						if ((dist < fiberDiameter)) {
							//then we have a possible link if no previous link was placed.

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
								//unsigned currentEdgeCount = currentEdgeCountVec[nodeId];

								//if (currentEdgeCount < maxNbrCount) {
									//then the new edge can be placed
									unsigned location;
									//decide placement of nbr only in row of nodeId,
									//another thread will place for candidateId
									for (unsigned place = possibleEdgeBegin; place < possibleEdgeEnd; place++ ) {

										if (maxNodeCount == globalNeighbors[place]) {
											//maxNodeCount is the default value of the matrix
											//then we have found location where we can place a node scine global neighbors default is maxNodeCount
											location = place;
											break;
										}

									}
									//if (location < (maxNbrCount * maxNodeCount) ) {
										globalNeighbors[location] = candidateId;
										lengthZero[location ] = dist;

										//only if a node is placed do we record it.
										//place in matrix format upper diagonal
										/*if (nodeId < candidateId) {
											//on row node_id (mat notation)
											finalEdgeIdPlaced = candidateId + maxNodeCount * nodeId;
										}
										else {
											//on row cand_id (mat notation)
											finalEdgeIdPlaced = nodeId + maxNodeCount * candidateId;
										};*/
										final_id_left = min(nodeId, candidateId);
										final_id_right = max(nodeId, candidateId);

										idTempVecLeft[save_index] = final_id_left;
										idTempVecRight[save_index] = final_id_right;
										(currentEdgeCountVec[nodeId])+=1;

										save_index += 1;
										numPlacedLinks += 1;
									//}
								//}
							}
						}
					}
				}
			}
		}
		return (numPlacedLinks);
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
	unsigned operator() (const unsigned& nodeId) {
		unsigned possibleEdgeBegin = maxNbrCount * nodeId;
		unsigned possibleEdgeEnd = maxNbrCount * nodeId + maxNbrCount;//test all possible for now

		unsigned numDelinked=0;

		for (unsigned index = possibleEdgeBegin; index < possibleEdgeEnd; index++) {
		//index is the location of nbrId

			unsigned nbrId = globalNeighbors[index];//loop through all possible neighbors
			if (nbrId < maxNodeCount) {
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
					globalNeighbors[index] = maxNodeCount;
					lengthZero[index] = 0.0;
					(currentEdgeCountVec[nodeId])-=1;//subtract from edgecount if nbr removed.
					numDelinked +=1;
				}

			}
		}
		return numDelinked;
	}

};

#endif /*LINKNODESONDEVICE_H_*/
