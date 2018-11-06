#include "NodeSystemDevice.h"
#include "WLCSolveOnDevice.h" 

void WLCSolveOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,  
	GeneralParams& generalParams) {
 
 
	thrust::counting_iterator<unsigned> startEdgeIter(0);
			  
	//
	thrust::for_each( 
		thrust::make_zip_iterator( 
			thrust::make_tuple(startEdgeIter,
								nodeInfoVecs.isNodeFixed.begin() )),
		thrust::make_zip_iterator(
			thrust::make_tuple(startEdgeIter,
								nodeInfoVecs.isNodeFixed.begin() )) + generalParams.maxNodeCount,
		WLCfunctor(
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data()),
 
			generalParams.kB,
			generalParams.persistenceLengthMon,
			generalParams.CLM,
			generalParams.temperature,
			generalParams.maxNeighborCount,
			generalParams.maxNodeCount,

			thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.numOriginalNeighborsNodeVector.data()) ) );
};

