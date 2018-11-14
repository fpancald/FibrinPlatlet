
#include "LinkNodesOnDevice.h"
#include "NodeSystemDevice.h"


void LinkNodesOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

		//Default fill values at 0.
		thrust::fill(nodeInfoVecs.linksThreadMade.begin(),
			nodeInfoVecs.linksThreadMade.end(), 0);
		
		thrust::fill(nodeInfoVecs.delinksThreadMade.begin(),
			nodeInfoVecs.delinksThreadMade.end(), 0);

		thrust::fill(nodeInfoVecs.idMadeTempLeft.begin(),
				nodeInfoVecs.idMadeTempLeft.end(), 0);

		thrust::fill(nodeInfoVecs.idMadeTempRight.begin(),
				nodeInfoVecs.idMadeTempRight.end(), 0);


		//unsigned globalcount = thrust::count_if(wlcInfoVecs.globalNeighbors.begin(),wlcInfoVecs.globalNeighbors.end(),is_less_than(generalParams.maxNodeCount));
		//std::cout<<"currentEdgeCount varpre: "<< generalParams.currentEdgeCount << std::endl;
		//std::cout<<"currentEdgeCount globalpre: "<< globalcount/2 << std::endl;

		thrust::counting_iterator<unsigned> counter(0);
		thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						counter,
						auxVecs.id_bucket.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						counter,
						auxVecs.id_bucket.begin())) + generalParams.maxNodeCount,
				nodeInfoVecs.linksThreadMade.begin(),//output
				LinkNodesFunctor(
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),

					thrust::raw_pointer_cast(auxVecs.id_value_expanded.data()),
					thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
					thrust::raw_pointer_cast(auxVecs.keyEnd.data()),

					generalParams.fiberDiameter,
					generalParams.maxNeighborCount,
					generalParams.maxNodeCount,

					generalParams.maxLinksPerIteration,
					thrust::raw_pointer_cast(nodeInfoVecs.idMadeTempLeft.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.idMadeTempRight.data()) ) );



		/*	for (unsigned i = 0; i < idMadeTempLeft.size(); i++) {
				unsigned varL = idMadeTempLeft[i];
				unsigned varR = idMadeTempRight[i];

				if ((varL != 0) || (varR != 0))
					std::cout<< varL << " " <<varR << std::endl;
			}
		unsigned begin = 479 * generalParams.maxNeighborCount;
		unsigned end = begin + generalParams.maxNeighborCount;
		for (unsigned i = begin; i < end; i++){
			unsigned id = wlcInfoVecs.globalNeighbors[i];
			if (id < generalParams.maxNodeCount){
				std::cout<<" 479: "<< id <<std::endl;
			}
		}
		begin = 1004 * generalParams.maxNeighborCount;
		end = begin + generalParams.maxNeighborCount;
		for (unsigned i = begin; i < end; i++){
			unsigned id = wlcInfoVecs.globalNeighbors[i];
			if (id < generalParams.maxNodeCount){
				std::cout<<" 1004: "<< id <<std::endl;
			}
		}*/
	/*	thrust::counting_iterator<unsigned> counterDeLink(0);

		thrust::transform(
						counterDeLink,
						counterDeLink + generalParams.maxNodeCount,
				nodeInfoVecs.delinksThreadMade.begin(),
				DeLinkCopiesFunctor(
					thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
					generalParams.maxNeighborCount,
					generalParams.maxNodeCount ) );
*/
		//sort increasing
		//issue here
		thrust::sort_by_key(
			nodeInfoVecs.idMadeTempLeft.begin(),nodeInfoVecs.idMadeTempLeft.end(),
			nodeInfoVecs.idMadeTempRight.begin(),thrust::greater<unsigned>() );

		thrust::stable_sort_by_key(
			nodeInfoVecs.idMadeTempRight.begin(),nodeInfoVecs.idMadeTempRight.end(),
			nodeInfoVecs.idMadeTempLeft.begin(), thrust::greater<unsigned>() );


		unsigned idL_init = nodeInfoVecs.idMadeTempLeft[0];
		unsigned idR_init = nodeInfoVecs.idMadeTempRight[0];

		unsigned count = 0;
		for (unsigned i = 1; i < nodeInfoVecs.idMadeTempLeft.size(); i++) {
			//add extra edges and preferred lengths. Notice the lower and upper must be added since each imparts force to one single node and
			//not the neighboring node to the edge. This is b/c edges are solved per node and not per edge
			unsigned idL = nodeInfoVecs.idMadeTempLeft[i];
			unsigned idR = nodeInfoVecs.idMadeTempRight[i];

			if ((idL == idL_init) && (idR == idR_init)){
				count +=1;
			}
			else {
				count = 0;
			}
			//reset initial id's
			idL_init = idL;
			idR_init = idR;


			if ( ((idL != 0) || (idR != 0) ) && (count == 1)) {

				//count edges
				//std::cout<<"placing id: "<< idL<<" " << idR<<std::endl;

				
				nodeInfoVecs.deviceEdgeLeft[generalParams.currentEdgeCount] = (idL);
				nodeInfoVecs.deviceEdgeRight[generalParams.currentEdgeCount] = (idR);
				generalParams.currentEdgeCount += 1;
			} 

		} 
 

	/*	unsigned globalcount = thrust::count_if(wlcInfoVecs.globalNeighbors.begin(), wlcInfoVecs.globalNeighbors.end(), is_less_than(generalParams.maxNodeCount));

		unsigned linksmade = *(thrust::max_element(linksThreadMade.begin(), linksThreadMade.end() ));
		unsigned delinksmade = *(thrust::max_element(delinksThreadMade.begin(), delinksThreadMade.end() ));
		std::cout<<"max links made this iteration: "<< linksmade << std::endl;
		std::cout<<"max unlinks made this iteration: "<< delinksmade << std::endl;

		std::cout<<"currentEdgeCount var: "<< generalParams.currentEdgeCount << std::endl;
		std::cout<<"currentEdgeCount global "<< globalcount/2 << std::endl;

		unsigned temp= thrust::reduce(	wlcInfoVecs.currentNodeEdgeCountVector.begin(),
			wlcInfoVecs.currentNodeEdgeCountVector.end());
		std::cout<<"currentEdgeCount dev: "<< temp << std::endl;
	*/




};
