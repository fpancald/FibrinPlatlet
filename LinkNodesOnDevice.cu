
#include "LinkNodesOnDevice.h"
#include "NodeSystemDevice.h"


void LinkNodesOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

		thrust::device_vector<unsigned> linksThreadMade;
		linksThreadMade.resize(generalParams.maxNodeCount);
		thrust::fill(linksThreadMade.begin(),linksThreadMade.end(), 0);
		thrust::device_vector<unsigned> delinksThreadMade;
		delinksThreadMade.resize(generalParams.maxNodeCount);
		thrust::fill(delinksThreadMade.begin(),delinksThreadMade.end(), 0);

		thrust::fill(nodeInfoVecs.idEdgesMadeTemp.begin(),
				nodeInfoVecs.idEdgesMadeTemp.end(), 0);

		thrust::device_vector<unsigned> id;
		thrust::device_vector<unsigned> idMadeTempLeft;
		thrust::device_vector<unsigned> idMadeTempRight;
		id.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);
		idMadeTempLeft.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);
		idMadeTempRight.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);

		thrust::fill(id.begin(),
				id.end(), 0);
		thrust::fill(idMadeTempLeft.begin(),
				idMadeTempLeft.end(), 0);
		thrust::fill(idMadeTempRight.begin(),
				idMadeTempRight.end(), 0);


		unsigned globalcount = thrust::count_if(wlcInfoVecs.globalNeighbors.begin(),wlcInfoVecs.globalNeighbors.end(),is_less_than(generalParams.maxNodeCount));
	//	std::cout<<"currentEdgeCount varpre: "<< generalParams.currentEdgeCount << std::endl;
	//	std::cout<<"currentEdgeCount globalpre: "<< globalcount/2 << std::endl;

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
				linksThreadMade.begin(),//output
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
					thrust::raw_pointer_cast(idMadeTempLeft.data()),
					thrust::raw_pointer_cast(idMadeTempRight.data()) ) );

			thrust::sort_by_key(
				idMadeTempLeft.begin(),idMadeTempLeft.end(),
				idMadeTempRight.begin(),thrust::greater<unsigned>() );

			thrust::stable_sort_by_key(
				idMadeTempRight.begin(),idMadeTempRight.end(),
				idMadeTempLeft.begin(), thrust::greater<unsigned>() );


		/*	for (unsigned i = 0; i < idMadeTempLeft.size(); i++) {
				unsigned varL = idMadeTempLeft[i];
				unsigned varR = idMadeTempRight[i];

				if ((varL != 0) || (varR != 0))
					std::cout<< varL << " " <<varR << std::endl;
			}*/
	/*	unsigned begin = 479 * generalParams.maxNeighborCount;
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
		thrust::counting_iterator<unsigned> counterDeLink(0);

		thrust::transform(
						counterDeLink,
						counterDeLink + generalParams.maxNodeCount,
				delinksThreadMade.begin(),
				DeLinkCopiesFunctor(
					thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
					generalParams.maxNeighborCount,
					generalParams.maxNodeCount ) );

    	unsigned endKey = thrust::get<0>(
    	    thrust::unique_by_key(
    	        idMadeTempRight.begin(), idMadeTempRight.end(),
				idMadeTempLeft.begin(),
    	    thrust::equal_to<unsigned>() )) - idMadeTempRight.begin();//binary_pred


		for (unsigned i = 0; i < endKey; i++) {
			//add extra edges and preferred lengths. Notice the lower and upper must be added since each imparts force to one single node and
			//not the neighboring node to the edge. This is b/c edges are solved per node and not per edge
			unsigned idL = idMadeTempLeft[i];
			unsigned idR = idMadeTempRight[i];


			if ((idL != 0) || (idR != 0) ) {
				//count edges
				//std::cout<<"placing id: "<< idL<<" " << idR<<std::endl;

				//idEdgesMadeHost contains id's in matrix format
				nodeInfoVecs.deviceEdgeLeft[generalParams.currentEdgeCount] = (idL);
				nodeInfoVecs.deviceEdgeRight[generalParams.currentEdgeCount] = (idR);
				generalParams.currentEdgeCount += 1;
			} 

		} 


		globalcount = thrust::count_if(wlcInfoVecs.globalNeighbors.begin(), wlcInfoVecs.globalNeighbors.end(), is_less_than(generalParams.maxNodeCount));

		unsigned linksmade = *(thrust::max_element(linksThreadMade.begin(), linksThreadMade.end() ));
		unsigned delinksmade = *(thrust::max_element(delinksThreadMade.begin(), delinksThreadMade.end() ));
	/*	std::cout<<"max links made this iteration: "<< linksmade << std::endl;
		std::cout<<"max unlinks made this iteration: "<< delinksmade << std::endl;

		std::cout<<"currentEdgeCount var: "<< generalParams.currentEdgeCount << std::endl;
		std::cout<<"currentEdgeCount global "<< globalcount/2 << std::endl;
*/




};
