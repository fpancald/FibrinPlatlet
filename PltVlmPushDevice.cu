#include "PltVlmPushDevice.h"
#include "PltForceDevice.h"
#include "NodeSystemDevice.h"


//Call the plt force on nodes functor
//This functor applies force to nodes from platelets 
//as well as the self platelet from other platelets to create volume exclusion
//The interaction count is: plt_other_intrct. No imaging is done here. 
void PltVlmPushOnDevice(
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


		thrust::counting_iterator<unsigned> counter(0);

		for (unsigned i = 0; i < auxVecs.idPlt_bucket.size(); i++)
			std::cout<<"plt bucketvol: "<<auxVecs.idPlt_bucket[i] << std::endl;

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
             PltVlmPushForceFunctor(
                generalParams.plt_other_intrct,
                generalParams.pltRForce,
                generalParams.pltForce,
                generalParams.pltR,

                generalParams.maxPltCount,
                generalParams.fiberDiameter,
		        generalParams.maxNodeCount,
                generalParams.maxNeighborCount,

                thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
                thrust::raw_pointer_cast(pltInfoVecs.nodeUnreducedForceX.data()),
                thrust::raw_pointer_cast(pltInfoVecs.nodeUnreducedForceY.data()),
                thrust::raw_pointer_cast(pltInfoVecs.nodeUnreducedForceZ.data()),

                thrust::raw_pointer_cast(pltInfoVecs.nodeUnreducedId.data()),

                thrust::raw_pointer_cast(auxVecs.id_value_expanded.data()),
                thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
                thrust::raw_pointer_cast(auxVecs.keyEnd.data()),
				 
                thrust::raw_pointer_cast(auxVecs.idPlt_value_expanded.data()),
                thrust::raw_pointer_cast(auxVecs.keyPltBegin.data()),
                thrust::raw_pointer_cast(auxVecs.keyPltEnd.data()),

                thrust::raw_pointer_cast(pltInfoVecs.pltLocX.data()),
                thrust::raw_pointer_cast(pltInfoVecs.pltLocY.data()),
                thrust::raw_pointer_cast(pltInfoVecs.pltLocZ.data())) );

        //now call a sort by key followed by a reduce by key to figure out which nodes are have force applied.
        //then make a functor that takes the id and force (4 tuple) and takes that force and adds it to the id'th entry in nodeInfoVecs.nodeForceX,Y,Z
        thrust::sort_by_key(pltInfoVecs.nodeUnreducedId.begin(), pltInfoVecs.nodeUnreducedId.end(),
        			thrust::make_zip_iterator(
        				thrust::make_tuple(
        					pltInfoVecs.nodeUnreducedForceX.begin(),
        					pltInfoVecs.nodeUnreducedForceY.begin(),
        					pltInfoVecs.nodeUnreducedForceZ.begin())), thrust::less<unsigned>());


		//reduce and apply force
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
	}