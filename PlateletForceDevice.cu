#include "PlateletForceDevice.h"
#include "NodeSystemDevice.h"

//Force field-like mode
void PltForceOnDevice(
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

		//fill for image sort
    	thrust::fill(pltInfoVecs.nodeUnreducedId.begin(),pltInfoVecs.nodeUnreducedId.end(), generalParams.maxNodeCount);
	/*	unsigned begin = auxVecs.keyBegin[auxVecs.idPlt_bucket[0]];
		unsigned end = auxVecs.keyEnd[auxVecs.idPlt_bucket[0]];
		for (unsigned i = begin; i < end; i++) {
			unsigned id = auxVecs.id_bucket_expanded[i];
			std::cout<<id<< std::endl;
		}*/

        //Call the plt force on nodes functor
        thrust::transform(
        	thrust::make_zip_iterator(
        		thrust::make_tuple(
   					auxVecs.idPlt_bucket.begin(),
   					auxVecs.idPlt_value.begin(),
        			pltInfoVecs.pltLocX.begin(),
        			pltInfoVecs.pltLocY.begin(),
        			pltInfoVecs.pltLocZ.begin())),
        	thrust::make_zip_iterator(
        		thrust::make_tuple(
        			auxVecs.idPlt_bucket.begin(),
    				auxVecs.idPlt_value.begin(),
        		 	pltInfoVecs.pltLocX.begin(),
        		 	pltInfoVecs.pltLocY.begin(),
        		 	pltInfoVecs.pltLocZ.begin())) + generalParams.maxPltCount,
         //save plt forces
         thrust::make_zip_iterator(
        	 thrust::make_tuple(
				 //reset's forces
        		 pltInfoVecs.pltForceX.begin(),
        		 pltInfoVecs.pltForceY.begin(),
        		 pltInfoVecs.pltForceZ.begin())),
             PltonNodeForceFunctor(
                 generalParams.pltMaxConn,
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
                 thrust::raw_pointer_cast(pltInfoVecs.pltImagingConnection.data()),

                 thrust::raw_pointer_cast(auxVecs.id_bucket_expanded.data()),
                 thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
                 thrust::raw_pointer_cast(auxVecs.keyEnd.data()) ) );

        //now call a sort by key followed by a reduce by key to figure out which nodes are have force applied.
        //then make a functor that takes the id and force (4 tuple) and takes that force and adds it to the id'th entry in nodeInfoVecs.nodeForceX,Y,Z
        thrust::sort_by_key(pltInfoVecs.nodeUnreducedId.begin(), pltInfoVecs.nodeUnreducedId.end(),
        			thrust::make_zip_iterator(
        				thrust::make_tuple(
							pltInfoVecs.pltImagingConnection.begin(),
        					pltInfoVecs.nodeUnreducedForceX.begin(),
        					pltInfoVecs.nodeUnreducedForceY.begin(),
        					pltInfoVecs.nodeUnreducedForceZ.begin())), thrust::less<unsigned>());

    thrust::copy(pltInfoVecs.nodeUnreducedId.begin(),pltInfoVecs.nodeUnreducedId.end(), pltInfoVecs.nodeImagingConnection.begin());

    pltInfoVecs.numConnections = thrust::count_if(
        pltInfoVecs.nodeImagingConnection.begin(),
        pltInfoVecs.nodeImagingConnection.end(), is_less_than(generalParams.maxNodeCount) );

	//std::cout<<pltInfoVecs.numConnections<<std::endl;
	/*for (unsigned i = 0; i < pltInfoVecs.nodeImagingConnection.size(); i ++) {
		std::cout<<pltInfoVecs.nodeImagingConnection[i]<<std::endl;
	}*/

	//std::cout<<"num connect: "<< pltInfoVecs.numConnections<<std::endl;
	//std::cout<<"num nod imaging: "<< pltInfoVecs.nodeImagingConnection.size() <<std::endl;
	//std::cout<<"num node unreduced: "<< pltInfoVecs.nodeUnreducedId.size() <<std::endl;
	//std::cout<<"num plt imaging: "<< pltInfoVecs.pltImagingConnection.size() <<std::endl;

	/*for (unsigned i = 0; i < pltInfoVecs.nodeUnreducedId.size(); i ++) {
		std::cout<<pltInfoVecs.nodeUnreducedId[i]<<std::endl;
		std::cout<<pltInfoVecs.nodeUnreducedForceX[i]<<std::endl;
		std::cout<<pltInfoVecs.nodeUnreducedForceY[i]<<std::endl;
		std::cout<<pltInfoVecs.nodeUnreducedForceZ[i]<<std::endl;
	}*/

/*	for (unsigned i = 0; i < nodeInfoVecs.nodeLocX.size(); i ++) {
		std::cout<<"node Loc "<< nodeInfoVecs.nodeLocX[i] << " " << nodeInfoVecs.nodeLocY[i] <<" "<< nodeInfoVecs.nodeLocZ[i] << std::endl;

	}*/
//	std::cout<<"plt Loc "<< pltInfoVecs.pltLocX[0] << " " << pltInfoVecs.pltLocY[0] <<" "<< pltInfoVecs.pltLocZ[0] << std::endl;
//	std::cout<<"plt Loc "<< pltInfoVecs.pltLocX[1] << " " << pltInfoVecs.pltLocY[1] <<" "<< pltInfoVecs.pltLocZ[1] << std::endl;


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

void PltInteractionOnDevice(
  	GeneralParams& generalParams,
  	PltInfoVecs& pltInfoVecs,
  	AuxVecs& auxVecs) {

    thrust::for_each(
      	thrust::make_zip_iterator(
        	thrust::make_tuple(
        		auxVecs.idPlt_bucket.begin(),
        		auxVecs.idPlt_value.begin(),
          		pltInfoVecs.pltLocX.begin(),
          		pltInfoVecs.pltLocY.begin(),
          		pltInfoVecs.pltLocZ.begin())),
    thrust::make_zip_iterator(
        thrust::make_tuple(
          		auxVecs.idPlt_bucket.begin(),
        		auxVecs.idPlt_value.begin(),
          		pltInfoVecs.pltLocX.begin(),
          		pltInfoVecs.pltLocY.begin(),
          		pltInfoVecs.pltLocZ.begin())) + generalParams.maxPltCount,

         PltonPltForceFunctor(
             generalParams.pltMaxConn,
             generalParams.pltRForce,
             generalParams.pltForce,
             generalParams.pltR,
             generalParams.maxPltCount,
             thrust::raw_pointer_cast(pltInfoVecs.pltLocX.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltLocY.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltLocZ.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltForceX.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltForceY.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltForceZ.data()),
             thrust::raw_pointer_cast(auxVecs.idPlt_value_expanded.data()),
             thrust::raw_pointer_cast(auxVecs.keyPltBegin.data()),
             thrust::raw_pointer_cast(auxVecs.keyPltEnd.data()) ) );
};



//tendril-like force
void PltTndrlOnDevice(
  NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
	PltInfoVecs& pltInfoVecs,
	AuxVecs& auxVecs) {

  if (generalParams.currentTime==0.0){
    thrust::fill(pltInfoVecs.tndrlNodeId.begin(),pltInfoVecs.tndrlNodeId.end(), generalParams.maxNodeCount + generalParams.maxPltCount);
    }
		thrust::fill(pltInfoVecs.nodeUnreducedForceX.begin(), pltInfoVecs.nodeUnreducedForceX.end(), 0.0);
		thrust::fill(pltInfoVecs.nodeUnreducedForceY.begin(), pltInfoVecs.nodeUnreducedForceY.end(), 0.0);
		thrust::fill(pltInfoVecs.nodeUnreducedForceZ.begin(), pltInfoVecs.nodeUnreducedForceZ.end(), 0.0);

		thrust::fill(pltInfoVecs.nodeReducedForceX.begin(), pltInfoVecs.nodeReducedForceX.end(), 0.0);
		thrust::fill(pltInfoVecs.nodeReducedForceY.begin(), pltInfoVecs.nodeReducedForceY.end(), 0.0);
		thrust::fill(pltInfoVecs.nodeReducedForceZ.begin(), pltInfoVecs.nodeReducedForceZ.end(), 0.0);

		//fill for image sort
    	thrust::fill(pltInfoVecs.nodeUnreducedId.begin(),pltInfoVecs.nodeUnreducedId.end(), generalParams.maxNodeCount);
	/*	unsigned begin = auxVecs.keyBegin[auxVecs.idPlt_bucket[0]];
		unsigned end = auxVecs.keyEnd[auxVecs.idPlt_bucket[0]];
		for (unsigned i = begin; i < end; i++) {
			unsigned id = auxVecs.id_bucket_expanded[i];
			std::cout<<id<< std::endl;
		}*/

        //Call the plt force on nodes functor
        thrust::transform(
        	thrust::make_zip_iterator(
        		thrust::make_tuple(
   					auxVecs.idPlt_bucket.begin(),
   					auxVecs.idPlt_value.begin(),
        			pltInfoVecs.pltLocX.begin(),
        			pltInfoVecs.pltLocY.begin(),
        			pltInfoVecs.pltLocZ.begin())),
        	thrust::make_zip_iterator(
        		thrust::make_tuple(
        			auxVecs.idPlt_bucket.begin(),
    				auxVecs.idPlt_value.begin(),
        		 	pltInfoVecs.pltLocX.begin(),
        		 	pltInfoVecs.pltLocY.begin(),
        		 	pltInfoVecs.pltLocZ.begin())) + generalParams.maxPltCount,
         //save plt forces
         thrust::make_zip_iterator(
        	 thrust::make_tuple(
				 //reset's forces
        		 pltInfoVecs.pltForceX.begin(),
        		 pltInfoVecs.pltForceY.begin(),
        		 pltInfoVecs.pltForceZ.begin())),
             PltTndrlonNodeForceFunctor(
                 generalParams.pltMaxConn,
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
                 thrust::raw_pointer_cast(pltInfoVecs.pltImagingConnection.data()),

                 thrust::raw_pointer_cast(auxVecs.id_bucket_expanded.data()),
                 thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
                 thrust::raw_pointer_cast(auxVecs.keyEnd.data()),

                 thrust::raw_pointer_cast(pltInfoVecs.tndrlNodeId.data()),
                 thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data())
                 thrust::raw_pointer_cast(nodeInfoVecs.pltLocX.data()),
                 thrust::raw_pointer_cast(nodeInfoVecs.pltLocY.data()),
                 thrust::raw_pointer_cast(nodeInfoVecs.pltLocZ.data()),) );

        //now call a sort by key followed by a reduce by key to figure out which nodes are have force applied.
        //then make a functor that takes the id and force (4 tuple) and takes that force and adds it to the id'th entry in nodeInfoVecs.nodeForceX,Y,Z
        thrust::sort_by_key(pltInfoVecs.nodeUnreducedId.begin(), pltInfoVecs.nodeUnreducedId.end(),
        			thrust::make_zip_iterator(
        				thrust::make_tuple(
							pltInfoVecs.pltImagingConnection.begin(),
        					pltInfoVecs.nodeUnreducedForceX.begin(),
        					pltInfoVecs.nodeUnreducedForceY.begin(),
        					pltInfoVecs.nodeUnreducedForceZ.begin())), thrust::less<unsigned>());

    thrust::copy(pltInfoVecs.nodeUnreducedId.begin(),pltInfoVecs.nodeUnreducedId.end(), pltInfoVecs.nodeImagingConnection.begin());

    pltInfoVecs.numConnections = thrust::count_if(
        pltInfoVecs.nodeImagingConnection.begin(),
        pltInfoVecs.nodeImagingConnection.end(), is_less_than(generalParams.maxNodeCount) );

	//std::cout<<pltInfoVecs.numConnections<<std::endl;
	/*for (unsigned i = 0; i < pltInfoVecs.nodeImagingConnection.size(); i ++) {
		std::cout<<pltInfoVecs.nodeImagingConnection[i]<<std::endl;
	}*/

	//std::cout<<"num connect: "<< pltInfoVecs.numConnections<<std::endl;
	//std::cout<<"num nod imaging: "<< pltInfoVecs.nodeImagingConnection.size() <<std::endl;
	//std::cout<<"num node unreduced: "<< pltInfoVecs.nodeUnreducedId.size() <<std::endl;
	//std::cout<<"num plt imaging: "<< pltInfoVecs.pltImagingConnection.size() <<std::endl;

	/*for (unsigned i = 0; i < pltInfoVecs.nodeUnreducedId.size(); i ++) {
		std::cout<<pltInfoVecs.nodeUnreducedId[i]<<std::endl;
		std::cout<<pltInfoVecs.nodeUnreducedForceX[i]<<std::endl;
		std::cout<<pltInfoVecs.nodeUnreducedForceY[i]<<std::endl;
		std::cout<<pltInfoVecs.nodeUnreducedForceZ[i]<<std::endl;
	}*/

/*	for (unsigned i = 0; i < nodeInfoVecs.nodeLocX.size(); i ++) {
		std::cout<<"node Loc "<< nodeInfoVecs.nodeLocX[i] << " " << nodeInfoVecs.nodeLocY[i] <<" "<< nodeInfoVecs.nodeLocZ[i] << std::endl;

	}*/
//	std::cout<<"plt Loc "<< pltInfoVecs.pltLocX[0] << " " << pltInfoVecs.pltLocY[0] <<" "<< pltInfoVecs.pltLocZ[0] << std::endl;
//	std::cout<<"plt Loc "<< pltInfoVecs.pltLocX[1] << " " << pltInfoVecs.pltLocY[1] <<" "<< pltInfoVecs.pltLocZ[1] << std::endl;


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

/*void PltTndrlInteractionOnDevice(
  	GeneralParams& generalParams,
  	PltInfoVecs& pltInfoVecs,
  	AuxVecs& auxVecs) {

    thrust::for_each(
      	thrust::make_zip_iterator(
        	thrust::make_tuple(
        		auxVecs.idPlt_bucket.begin(),
        		auxVecs.idPlt_value.begin(),
          		pltInfoVecs.pltLocX.begin(),
          		pltInfoVecs.pltLocY.begin(),
          		pltInfoVecs.pltLocZ.begin())),
    thrust::make_zip_iterator(
        thrust::make_tuple(
          		auxVecs.idPlt_bucket.begin(),
        		auxVecs.idPlt_value.begin(),
          		pltInfoVecs.pltLocX.begin(),
          		pltInfoVecs.pltLocY.begin(),
          		pltInfoVecs.pltLocZ.begin())) + generalParams.maxPltCount,

         PltonPltForceFunctor(
             generalParams.pltMaxConn,
             generalParams.pltRForce,
             generalParams.pltForce,
             generalParams.pltR,
             generalParams.maxPltCount,
             thrust::raw_pointer_cast(pltInfoVecs.pltLocX.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltLocY.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltLocZ.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltForceX.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltForceY.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltForceZ.data()),
             thrust::raw_pointer_cast(auxVecs.idPlt_value_expanded.data()),
             thrust::raw_pointer_cast(auxVecs.keyPltBegin.data()),
             thrust::raw_pointer_cast(auxVecs.keyPltEnd.data()) ) );
};*/
