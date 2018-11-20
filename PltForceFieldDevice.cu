#include "PltForceFieldDevice.h"
#include "PltForceDevice.h"
#include "NodeSystemDevice.h"

//for a given platelet, search for network nodes, pull and push them
//Force field-like mode

void PltForceFieldOnDevice(
  	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
	PltInfoVecs& pltInfoVecs,
	AuxVecs& auxVecs) {


		thrust::fill(pltInfoVecs.nodeUnreducedForceX.begin(), pltInfoVecs.nodeUnreducedForceX.end(), 0.0);
		thrust::fill(pltInfoVecs.nodeUnreducedForceY.begin(), pltInfoVecs.nodeUnreducedForceY.end(), 0.0);
		thrust::fill(pltInfoVecs.nodeUnreducedForceZ.begin(), pltInfoVecs.nodeUnreducedForceZ.end(), 0.0);

		thrust::fill(pltInfoVecs.nodeReducedForceX.begin(), pltInfoVecs.nodeReducedForceX.end(), 0.0);
		thrust::fill(pltInfoVecs.nodeReducedForceY.begin(), pltInfoVecs.nodeReducedForceY.end(), 0.0);
		thrust::fill(pltInfoVecs.nodeReducedForceZ.begin(), pltInfoVecs.nodeReducedForceZ.end(), 0.0);



        //Call the plt force on nodes functor
		thrust::counting_iterator<unsigned> counter(0);

        thrust::transform(
        	thrust::make_zip_iterator(
        		thrust::make_tuple(
					counter,
   					auxVecs.idPlt_bucket.begin(),
        			pltInfoVecs.pltLocX.begin(),
        			pltInfoVecs.pltLocY.begin(),
        			pltInfoVecs.pltLocZ.begin(),
        		 	pltInfoVecs.pltForceX.begin(),
        		 	pltInfoVecs.pltForceY.begin(),
        		 	pltInfoVecs.pltForceZ.begin())),
        	thrust::make_zip_iterator(
        		thrust::make_tuple(
					counter,
    				auxVecs.idPlt_bucket.begin(),
        		 	pltInfoVecs.pltLocX.begin(),
        		 	pltInfoVecs.pltLocY.begin(),
        		 	pltInfoVecs.pltLocZ.begin(), 
        		 	pltInfoVecs.pltForceX.begin(),
        		 	pltInfoVecs.pltForceY.begin(),
        		 	pltInfoVecs.pltForceZ.begin())) + generalParams.maxPltCount,
         //save plt forces
         thrust::make_zip_iterator(
        	 thrust::make_tuple(
				 //DOES NOT RESET FORCE
        		 pltInfoVecs.pltForceX.begin(),
        		 pltInfoVecs.pltForceY.begin(),
        		 pltInfoVecs.pltForceZ.begin())),
             PltonNodeForceFunctor(
                 generalParams.plt_other_intrct,
                 generalParams.pltRForce,
                 generalParams.pltForce,
                 generalParams.pltR,

                 generalParams.maxPltCount,
                 generalParams.fiberDiameter,
				 generalParams.maxNodeCount,

                 thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
                 thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
                 thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
                 thrust::raw_pointer_cast(pltInfoVecs.nodeUnreducedForceX.data()),
                 thrust::raw_pointer_cast(pltInfoVecs.nodeUnreducedForceY.data()),
                 thrust::raw_pointer_cast(pltInfoVecs.nodeUnreducedForceZ.data()),

                 thrust::raw_pointer_cast(pltInfoVecs.nodeUnreducedId.data()),

                 thrust::raw_pointer_cast(auxVecs.id_value_expanded.data()),//network neighbors
                 thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
                 thrust::raw_pointer_cast(auxVecs.keyEnd.data()) ) );

        //now call a sort by key followed by a reduce by key to figure out which nodes are have force applied.
        //then make a functor that takes the id and force (4 tuple) and takes that force and adds it to the id'th entry in nodeInfoVecs.nodeForceX,Y,Z
        thrust::sort_by_key(pltInfoVecs.nodeUnreducedId.begin(), pltInfoVecs.nodeUnreducedId.end(),
        			thrust::make_zip_iterator(
        				thrust::make_tuple(
        					pltInfoVecs.nodeUnreducedForceX.begin(),
        					pltInfoVecs.nodeUnreducedForceY.begin(),
        					pltInfoVecs.nodeUnreducedForceZ.begin())), thrust::less<unsigned>());

 //   thrust::copy(pltInfoVecs.nodeUnreducedId.begin(),pltInfoVecs.nodeUnreducedId.end(), pltInfoVecs.nodeImagingConnection.begin());

 //   pltInfoVecs.numConnections = thrust::count_if(
  //      pltInfoVecs.nodeImagingConnection.begin(),
  //      pltInfoVecs.nodeImagingConnection.end(), is_less_than(generalParams.maxNodeCount) );


		//reduce network force
 		unsigned endKey = thrust::get<0>(
 			thrust::reduce_by_key(
 				pltInfoVecs.nodeUnreducedId.begin(),
 				pltInfoVecs.nodeUnreducedId.end(),
 			thrust::make_zip_iterator(
 				thrust::make_tuple(
 					pltInfoVecs.nodeUnreducedForceX.begin(),
 					pltInfoVecs.nodeUnreducedForceY.begin(),
 					pltInfoVecs.nodeUnreducedForceZ.begin())),
 			pltInfoVecs.nodeReducedId.begin(),
 			thrust::make_zip_iterator(
 				thrust::make_tuple(//need t check
 					pltInfoVecs.nodeReducedForceX.begin(),
 					pltInfoVecs.nodeReducedForceY.begin(),
 					pltInfoVecs.nodeReducedForceZ.begin())),
 			thrust::equal_to<unsigned>(), CVec3Add())) - pltInfoVecs.nodeReducedId.begin();//binary_pred, binary_op

		//apply force to network
        thrust::for_each(
        	thrust::make_zip_iterator(//1st begin
        		thrust::make_tuple(
        			pltInfoVecs.nodeReducedId.begin(),
        			pltInfoVecs.nodeReducedForceX.begin(),
        			pltInfoVecs.nodeReducedForceY.begin(),
        			pltInfoVecs.nodeReducedForceZ.begin())),
        	thrust::make_zip_iterator(//1st end
        		thrust::make_tuple(
        			pltInfoVecs.nodeReducedId.begin(),
        			pltInfoVecs.nodeReducedForceX.begin(),
        			pltInfoVecs.nodeReducedForceY.begin(),
        			pltInfoVecs.nodeReducedForceZ.begin())) + endKey,
        	AddPltonNodeForceFunctor(
        		thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
        		thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
        		thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data())));
};