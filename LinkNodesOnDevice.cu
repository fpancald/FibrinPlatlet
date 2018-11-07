
#include "LinkNodesOnDevice.h"
#include "NodeSystemDevice.h"


void LinkNodesOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {
 
 
		thrust::for_each(  
				thrust::make_zip_iterator(
					thrust::make_tuple(
						auxVecs.bucketKeys.begin(), 
						auxVecs.bucketValues.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						auxVecs.bucketKeys.begin(),
						auxVecs.bucketValues.begin())) + generalParams.maxNodeCount,
				
				LinkNodesFunctor(
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),
					thrust::raw_pointer_cast(auxVecs.bucketKeys.data()),
					thrust::raw_pointer_cast(auxVecs.bucketValues.data()),
					thrust::raw_pointer_cast(auxVecs.bucketValuesIncludingNeighbor.data()),
					thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
					thrust::raw_pointer_cast(auxVecs.keyEnd.data()), 
					generalParams.fiberDiameter,
					generalParams.maxNeighborCount,
					generalParams.maxNodeCount, 
					generalParams.maxLinksPerIteration, 
					thrust::raw_pointer_cast(nodeInfoVecs.idEdgesMadeTemp.data()) ) );
		

		unsigned numAdded = thrust::count_if( nodeInfoVecs.idEdgesMadeTemp.begin(), nodeInfoVecs.idEdgesMadeTemp.end(), isNotEqualToZero() );
		if (numAdded != 0) {
			//be aware that if an edge was made, it will appear twice, but if it was made once and removed, it will appear once.
			//sort in increasing order. Then when we hit zero, we are done.  
			thrust::sort(nodeInfoVecs.idEdgesMadeTemp.begin(), nodeInfoVecs.idEdgesMadeTemp.end(),thrust::greater<unsigned>() ); 
			
			unsigned idLast = nodeInfoVecs.idEdgesMadeTemp[0];
			//std::cout<<"numadded "<< numAdded<<std::endl;
			//std::cout<<"idlast "<< idLast<<std::endl;
			//std::cout<<"idnext "<< nodeInfoVecs.idEdgesMadeTemp[1] <<std::endl;
			
			unsigned count = 0;
			for (unsigned i = 1; i<nodeInfoVecs.idEdgesMadeTemp.size(); i++) {
				//add extra edges and preferred lengths. Notice the lower and upper must be added since each imparts force to one single node and 
				//not the neighboring node to the edge. This is b/c edges are solved per node and not per edge
				unsigned id = nodeInfoVecs.idEdgesMadeTemp[i];
				
				if (id == 0){
					break; 
				}
	
				if (id == idLast) { 
					//then id has shown up twice and can be added
					count +=1; 
				}  
				else { 
					count = 0;
				}
	
				if ((id != 0) && (count > 0) ) {
					//count edges 
					//std::cout<<"placing id: "<< id<<std::endl;
				
					//idEdgesMadeHost contains id's in matrix format
					nodeInfoVecs.idEdgesMadeHost.push_back(id);
					unsigned idLeft = id % generalParams.maxNodeCount;
					unsigned idRight = id / generalParams.maxNodeCount;
					nodeInfoVecs.deviceEdgeLeft[generalParams.currentEdgeCount] = (idLeft);
					nodeInfoVecs.deviceEdgeRight[generalParams.currentEdgeCount] = (idRight);
					generalParams.currentEdgeCount += 1;
				} 
				
				idLast = id;//set last id to current
			} 
			
		}																					
		

		thrust::fill(nodeInfoVecs.idEdgesMadeTemp.begin(), 
				nodeInfoVecs.idEdgesMadeTemp.end(), 0);
		
};
