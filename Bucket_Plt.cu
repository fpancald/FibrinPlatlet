//This file sets the grid for network interaction with plts. 
#include "Bucket_Plt.h"
#include "System.h"

#include "functor_neighbor.h"
#include "functor_bucket_indexer.h"
#include "function_extend.h"

void init_plt_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	//always set bucket count. Update total if different. 
	domainParams.XBucketCount_plt_intc = ceil((domainParams.maxX - domainParams.minX) / domainParams.gridSpacing_plt_intc) + 1;
	domainParams.YBucketCount_plt_intc = ceil((domainParams.maxY - domainParams.minY) / domainParams.gridSpacing_plt_intc) + 1;
	domainParams.ZBucketCount_plt_intc = ceil((domainParams.maxZ - domainParams.minZ) / domainParams.gridSpacing_plt_intc) + 1;

	if ( (domainParams.XBucketCount_plt_intc * domainParams.YBucketCount_plt_intc * domainParams.ZBucketCount_plt_intc) != domainParams.totalBucketCount_plt_intc) {
		std::cout<<"resetting plt intct"<< std::endl;
        std::cout<<"x-bucket: "<< domainParams.XBucketCount_plt_intc<<std::endl;
		std::cout<<"y-bucket: "<< domainParams.YBucketCount_plt_intc<<std::endl;
		std::cout<<"z-bucket: "<< domainParams.ZBucketCount_plt_intc<<std::endl;
		//double amount of buckets in case of resizing networks
		domainParams.totalBucketCount_plt_intc = domainParams.XBucketCount_plt_intc * domainParams.YBucketCount_plt_intc * domainParams.ZBucketCount_plt_intc;
		std::cout<<"grid: "<< domainParams.gridSpacing_plt_intc << std::endl;
		std::cout<<"total bucket count: "<< domainParams.totalBucketCount_plt_intc<<std::endl;

		auxVecs.keyBegin_plt_intc.resize(domainParams.totalBucketCount_plt_intc);
		auxVecs.keyEnd_plt_intc.resize(domainParams.totalBucketCount_plt_intc);
		
        //platelets
		auxVecs.keyPltBegin.resize(domainParams.totalBucketCount_plt_intc); 
		auxVecs.keyPltEnd.resize(domainParams.totalBucketCount_plt_intc); 
 
	}

	thrust::fill(auxVecs.keyBegin_plt_intc.begin(),auxVecs.keyBegin_plt_intc.end(),0);
	thrust::fill(auxVecs.keyEnd_plt_intc.begin(),auxVecs.keyEnd_plt_intc.end(),0);
	//platelets
	thrust::fill(auxVecs.keyPltBegin.begin(),auxVecs.keyPltBegin.end(),0);
	thrust::fill(auxVecs.keyPltEnd.begin(),auxVecs.keyPltEnd.end(),0);

};

//convert buckets into neighboring scheme
void extend_plt_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	//memory is already allocated.
	unsigned endIndexExpanded = (auxVecs.endIndexBucketKeys_plt_intc) * 27;
	//platelets
	unsigned endIndexPltExpanded = (auxVecs.endIndexBucketPltKeys_plt_intc) * 27;

	//test for removing copies.
	unsigned valuesCount = auxVecs.id_value_plt_intc.size();
	thrust::fill(auxVecs.id_bucket_expanded_plt_intc.begin(),auxVecs.id_bucket_expanded_plt_intc.end(),0);
	thrust::fill(auxVecs.id_value_expanded_plt_intc.begin(),auxVecs.id_value_expanded_plt_intc.end(),0);

	thrust::fill(auxVecs.idPlt_bucket_expanded.begin(),auxVecs.idPlt_bucket_expanded.end(),0);
	thrust::fill(auxVecs.idPlt_value_expanded.begin(),auxVecs.idPlt_value_expanded.end(),0);




	/*
	* beginning of constant iterator
	*/
	thrust::constant_iterator<unsigned> first(27);
	/*
	* end of constant iterator.
	* the plus sign only indicate movement of position, not value.
	* e.g. movement is 5 and first iterator is initialized as 9
	* result array is [9,9,9,9,9];
	*/
	thrust::constant_iterator<unsigned> last = first + (auxVecs.endIndexBucketKeys_plt_intc); // this is NOT numerical addition!

	expand(first, last,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_plt_intc.begin(),
				auxVecs.id_value_plt_intc.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded_plt_intc.begin(),
				auxVecs.id_value_expanded_plt_intc.begin())));

	thrust::counting_iterator<unsigned> countingBegin(0);
 
	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded_plt_intc.begin(),
				countingBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded_plt_intc.begin(),
				countingBegin)) + endIndexExpanded,
		
		auxVecs.id_bucket_expanded_plt_intc.begin(),
		functor_neighbor(
			domainParams.XBucketCount_plt_intc,
			domainParams.YBucketCount_plt_intc,
			domainParams.ZBucketCount_plt_intc)); 

	thrust::stable_sort_by_key(auxVecs.id_bucket_expanded_plt_intc.begin(),
		auxVecs.id_bucket_expanded_plt_intc.end(),
		auxVecs.id_value_expanded_plt_intc.begin());


	thrust::counting_iterator<unsigned> search_begin(0);

	thrust::lower_bound(auxVecs.id_bucket_expanded_plt_intc.begin(),
		auxVecs.id_bucket_expanded_plt_intc.end(), search_begin,
		search_begin + domainParams.totalBucketCount_plt_intc,
		auxVecs.keyBegin_plt_intc.begin());

	thrust::upper_bound(auxVecs.id_bucket_expanded_plt_intc.begin(),
		auxVecs.id_bucket_expanded_plt_intc.end(),search_begin,
		search_begin + domainParams.totalBucketCount_plt_intc,
		auxVecs.keyEnd_plt_intc.begin());

	/*
	* beginning of constant iterator
	*/
	thrust::constant_iterator<unsigned> pltfirst(27);
	/**
	* end of constant iterator.
	* the plus sign only indicate movement of position, not value.
	* e.g. movement is 5 and first iterator is initialized as 9
	* result array is [9,9,9,9,9];
	*/


	thrust::constant_iterator<unsigned> pltlast = pltfirst + (auxVecs.endIndexBucketPltKeys_plt_intc); // this is NOT numerical addition!

	expand(pltfirst, pltlast,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.idPlt_bucket.begin(),
				auxVecs.idPlt_value.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.idPlt_bucket_expanded.begin(),
				auxVecs.idPlt_value_expanded.begin())));


	thrust::counting_iterator<unsigned> pltcountingBegin(0);

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.idPlt_bucket_expanded.begin(),
				pltcountingBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.idPlt_bucket_expanded.begin(),
				pltcountingBegin)) + endIndexPltExpanded,
			auxVecs.idPlt_bucket_expanded.begin(),
		functor_neighbor(
			domainParams.XBucketCount_plt_intc,
			domainParams.YBucketCount_plt_intc,
			domainParams.ZBucketCount_plt_intc));



	//unsigned pltnumberOfOutOfRange = thrust::count_if(auxVecs.idPlt_bucket_expanded.begin(),
	//	auxVecs.idPlt_bucket_expanded.end(), is_greater_than(domainParams.totalBucketCount) );
	//unsigned pltnumberInsideRange = endIndexPltExpanded - pltnumberOfOutOfRange;

	//unsigned endIndexPltSearch = endIndexPltExpanded - pltnumberOfOutOfRange;

	thrust::sort_by_key(auxVecs.idPlt_bucket_expanded.begin(),
		auxVecs.idPlt_bucket_expanded.end(),
		auxVecs.idPlt_value_expanded.begin());
	
	thrust::counting_iterator<unsigned> pltsearch_begin(0);

	thrust::lower_bound(auxVecs.idPlt_bucket_expanded.begin(),
		auxVecs.idPlt_bucket_expanded.end(), pltsearch_begin,
		pltsearch_begin + domainParams.totalBucketCount_plt_intc,
		auxVecs.keyPltBegin.begin());

	thrust::upper_bound(auxVecs.idPlt_bucket_expanded.begin(),
		auxVecs.idPlt_bucket_expanded.end(),pltsearch_begin,
		pltsearch_begin + domainParams.totalBucketCount_plt_intc,
		auxVecs.keyPltEnd.begin());

	/*
	unsigned choice = 0;

	unsigned bucket = auxVecs.idPlt_bucket[choice];
	std::cout<<"bucketplt 0: "<< bucket<<std::endl;
	std::cout<<"plt pos: "<<pltInfoVecs.pltLocX[0]<<" "<<pltInfoVecs.pltLocY[0]<<" "<<pltInfoVecs.pltLocZ[0]<<std::endl;
	std::cout<<"key len: "<< auxVecs.keyBegin.size() << std::endl;
	unsigned begin = auxVecs.keyBegin[bucket];
	unsigned end = auxVecs.keyEnd[bucket];
	
	std::cout<<"from bucket scheme:"<<std::endl;
	for (unsigned i = begin; i < end; i++) {
		
		unsigned nbr = auxVecs.id_value_expanded[i];
		unsigned buck = auxVecs.id_bucket[nbr];
		double x_dist = pltInfoVecs.pltLocX[choice] - nodeInfoVecs.nodeLocX[nbr];
		double y_dist = pltInfoVecs.pltLocY[choice] - nodeInfoVecs.nodeLocY[nbr];
		double z_dist = pltInfoVecs.pltLocZ[choice] - nodeInfoVecs.nodeLocZ[nbr];
		double dist = std::sqrt(std::pow(x_dist,2.0)+std::pow(y_dist,2.0)+std::pow(z_dist,2.0));
		if (dist < 1.0){
			std::cout<<"dist: "<< dist<< " between: "<< choice << " and nbr: "<< nbr<<std::endl; 
			std::cout<<"nbr: "<< nbr<< " is in bucket: "<< buck <<std::endl;
		}
	}*/

	/*
	std::cout<<"from all plt:"<<std::endl;
	for (unsigned i = 0; i < generalParams.maxNodeCount; i++) {
		unsigned nbr = i;//auxVecs.id_value_expanded[i];
		unsigned buck = auxVecs.id_bucket[nbr];
		double x_dist = pltInfoVecs.pltLocX[choice] - nodeInfoVecs.nodeLocX[nbr];
		double y_dist = pltInfoVecs.pltLocY[choice] - nodeInfoVecs.nodeLocY[nbr];
		double z_dist = pltInfoVecs.pltLocZ[choice] - nodeInfoVecs.nodeLocZ[nbr];
		double dist = std::sqrt(std::pow(x_dist,2.0)+std::pow(y_dist,2.0)+std::pow(z_dist,2.0));
		if (dist < 1.0){
			std::cout<<"dist: "<< dist<< " between: "<< choice << " and nbr: "<< nbr<<std::endl; 
			std::cout<<"nbr: "<< nbr<< " is in bucket: "<< buck <<std::endl;
		} 
	}*/


}


//takes nodes and places in buckets.
void build_plt_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {


	thrust::counting_iterator<unsigned> indexBucketBegin(0);
	// takes counting iterator and coordinates
	// return tuple of keys and values
	// transform the points to their bucket indices

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.nodeLocX.begin(),
				nodeInfoVecs.nodeLocY.begin(),
				nodeInfoVecs.nodeLocZ.begin(),
				indexBucketBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.nodeLocX.begin(),
				nodeInfoVecs.nodeLocY.begin(),
				nodeInfoVecs.nodeLocZ.begin(),
				indexBucketBegin)) + generalParams.maxNodeCount,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_plt_intc.begin(),
				auxVecs.id_value_plt_intc.begin())),
		functor_bucket_indexer(
			domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.XBucketCount_plt_intc,
            domainParams.YBucketCount_plt_intc,
            domainParams.ZBucketCount_plt_intc,
			domainParams.gridSpacing_plt_intc));

//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.id_value_plt_intc.begin(),
		auxVecs.id_value_plt_intc.begin() + generalParams.maxNodeCount,
		auxVecs.id_bucket_plt_intc.begin());

auxVecs.endIndexBucketKeys_plt_intc = generalParams.maxNodeCount;

	//platelets
	//std::cout<<"bucket platelet"<<std::endl;
	thrust::counting_iterator<unsigned> indexBucketBegin1(0);
	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				pltInfoVecs.pltLocX.begin(),
				pltInfoVecs.pltLocY.begin(),
				pltInfoVecs.pltLocZ.begin(),
				indexBucketBegin1)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				pltInfoVecs.pltLocX.begin(),
				pltInfoVecs.pltLocY.begin(),
				pltInfoVecs.pltLocZ.begin(),
				indexBucketBegin1)) + generalParams.maxPltCount,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.idPlt_bucket.begin(),
				auxVecs.idPlt_value.begin())),
		functor_bucket_indexer(
			domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.XBucketCount_plt_intc,
            domainParams.YBucketCount_plt_intc,
            domainParams.ZBucketCount_plt_intc,
			domainParams.gridSpacing_plt_intc));


	//std::cout<<"end bucket platelet"<<std::endl;
//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.idPlt_value.begin(),
		auxVecs.idPlt_value.end(),
		auxVecs.idPlt_bucket.begin());

auxVecs.endIndexBucketPltKeys_plt_intc = generalParams.maxPltCount;
 
};
