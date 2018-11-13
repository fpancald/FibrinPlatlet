
#include "BucketSchemeOnDevice.h"
#include "NodeSystemDevice.h"

//take domain and discretize into square buckets of size gridspace
void initDimensionBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	double minXTemp = (*(thrust::min_element(nodeInfoVecs.nodeLocX.begin(), nodeInfoVecs.nodeLocX.end())));
	double maxXTemp = (*(thrust::max_element(nodeInfoVecs.nodeLocX.begin(), nodeInfoVecs.nodeLocX.end())));
	double minYTemp = (*(thrust::min_element(nodeInfoVecs.nodeLocY.begin(), nodeInfoVecs.nodeLocY.end())));
	double maxYTemp = (*(thrust::max_element(nodeInfoVecs.nodeLocY.begin(), nodeInfoVecs.nodeLocY.end())));
	double minZTemp = (*(thrust::min_element(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end())));
	double maxZTemp = (*(thrust::max_element(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end())));

	//platelets
	domainParams.pltminX = (*(thrust::min_element(pltInfoVecs.pltLocX.begin(), pltInfoVecs.pltLocX.end())));
	domainParams.pltmaxX = (*(thrust::max_element(pltInfoVecs.pltLocX.begin(), pltInfoVecs.pltLocX.end())));
	domainParams.pltminY = (*(thrust::min_element(pltInfoVecs.pltLocY.begin(), pltInfoVecs.pltLocY.end())));
	domainParams.pltmaxY = (*(thrust::max_element(pltInfoVecs.pltLocY.begin(), pltInfoVecs.pltLocY.end())));
	domainParams.pltminZ = (*(thrust::min_element(pltInfoVecs.pltLocZ.begin(), pltInfoVecs.pltLocZ.end())));
	domainParams.pltmaxZ = (*(thrust::max_element(pltInfoVecs.pltLocZ.begin(), pltInfoVecs.pltLocZ.end())));

	double space = 0.0;
	domainParams.minX = min(minXTemp, domainParams.pltminX) - space;
	domainParams.maxX = max(maxXTemp, domainParams.pltmaxX) + space;
	domainParams.minY = min(minYTemp, domainParams.pltminY) - space;
	domainParams.maxY = max(maxYTemp, domainParams.pltmaxY) + space;
	domainParams.minZ = min(minZTemp, domainParams.pltminZ) - space;
	domainParams.maxZ = max(maxZTemp, domainParams.pltmaxZ) + space;

	domainParams.XBucketCount = (ceil(domainParams.maxX - domainParams.minX) / domainParams.gridSpacing) + 1;
	domainParams.YBucketCount = (ceil(domainParams.maxY - domainParams.minY) / domainParams.gridSpacing) + 1;
	domainParams.ZBucketCount = (ceil(domainParams.maxZ - domainParams.minZ) / domainParams.gridSpacing) + 1;

	if ( (domainParams.XBucketCount * domainParams.YBucketCount * domainParams.ZBucketCount) != domainParams.totalBucketCount	) {

		//double amount of buckets in case of resizing networks
		domainParams.totalBucketCount = domainParams.XBucketCount * domainParams.YBucketCount * domainParams.ZBucketCount;
		std::cout<<"grid: "<< domainParams.gridSpacing << std::endl;
		std::cout<<"total bucket count: "<< domainParams.totalBucketCount<<std::endl;

		auxVecs.keyBegin.resize(domainParams.totalBucketCount);
		auxVecs.keyEnd.resize(domainParams.totalBucketCount);
		//platelets
		auxVecs.keyPltBegin.resize(domainParams.totalBucketCount);
		auxVecs.keyPltEnd.resize(domainParams.totalBucketCount);

	}

	thrust::fill(auxVecs.keyBegin.begin(),auxVecs.keyBegin.end(),0);
	thrust::fill(auxVecs.keyEnd.begin(),auxVecs.keyEnd.end(),0);
	//platelets
	thrust::fill(auxVecs.keyPltBegin.begin(),auxVecs.keyPltBegin.end(),0);
	thrust::fill(auxVecs.keyPltEnd.begin(),auxVecs.keyPltEnd.end(),0);

};

//convert buckets into neighboring scheme
void extendBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs) {

	//memory is already allocated.
	unsigned endIndexExpanded = (auxVecs.endIndexBucketKeys) * 27;
	//platelets
	unsigned endIndexPltExpanded = (auxVecs.endIndexBucketPltKeys) * 27;

	//test for removing copies.
	unsigned valuesCount = auxVecs.id_value.size();
	thrust::fill(auxVecs.id_bucket_expanded.begin(),auxVecs.id_bucket_expanded.end(),0);
	thrust::fill(auxVecs.id_value_expanded.begin(),auxVecs.id_value_expanded.end(),0);

	thrust::fill(auxVecs.idPlt_bucket_expanded.begin(),auxVecs.idPlt_bucket_expanded.end(),0);
	thrust::fill(auxVecs.idPlt_value_expanded.begin(),auxVecs.idPlt_value_expanded.end(),0);




	/*
	* beginning of constant iterator
	*/
	thrust::constant_iterator<unsigned> first(27);
	/**
	* end of constant iterator.
	* the plus sign only indicate movement of position, not value.
	* e.g. movement is 5 and first iterator is initialized as 9
	* result array is [9,9,9,9,9];
	*/
	thrust::constant_iterator<unsigned> last = first + (auxVecs.endIndexBucketKeys); // this is NOT numerical addition!

	expand(first, last,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket.begin(),
				auxVecs.id_value.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded.begin(),
				auxVecs.id_value_expanded.begin())));


	thrust::counting_iterator<unsigned> countingBegin(0);

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded.begin(),
				countingBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded.begin(),
				countingBegin)) + endIndexExpanded,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded.begin(),
				countingBegin)),
		NeighborFunctor(
			domainParams.XBucketCount,
			domainParams.YBucketCount,
			domainParams.ZBucketCount));



	unsigned numberOfOutOfRange = thrust::count_if(auxVecs.id_bucket_expanded.begin(),
		auxVecs.id_bucket_expanded.end(), is_greater_than(domainParams.totalBucketCount) );
	unsigned numberInsideRange = endIndexExpanded - numberOfOutOfRange;

	//unsigned endIndexSearch = endIndexExpanded - numberOfOutOfRange;

	thrust::stable_sort_by_key(auxVecs.id_bucket_expanded.begin(),
		auxVecs.id_bucket_expanded.begin() + endIndexExpanded,
		auxVecs.id_value_expanded.begin());

	numberInsideRange =
		thrust::get<0>(thrust::unique_by_key(auxVecs.id_value_expanded.begin(),
			auxVecs.id_value_expanded.begin() + endIndexExpanded,
			auxVecs.id_bucket_expanded.begin())) - auxVecs.id_value_expanded.begin();

	auxVecs.id_bucket_expanded.erase(
			auxVecs.id_bucket_expanded.begin() + numberInsideRange,
			auxVecs.id_bucket_expanded.end());

	auxVecs.id_value_expanded.erase(
			auxVecs.id_value_expanded.begin() + numberInsideRange,
			auxVecs.id_value_expanded.end());



	thrust::counting_iterator<unsigned> search_begin(0);

	thrust::lower_bound(auxVecs.id_bucket_expanded.begin(),
		auxVecs.id_bucket_expanded.end(), search_begin,
		search_begin + domainParams.totalBucketCount,
		auxVecs.keyBegin.begin());

	thrust::upper_bound(auxVecs.id_bucket_expanded.begin(),
		auxVecs.id_bucket_expanded.end(),search_begin,
		search_begin + domainParams.totalBucketCount,
		auxVecs.keyEnd.begin());

	//platelets
	/*unsigned valuesPltCount = auxVecs.idPlt_value.size();
	thrust::fill(auxVecs.idPlt_bucket_expanded.begin(),auxVecs.idPlt_bucket_expanded.end(),0);
	thrust::fill(auxVecs.idPlt_value_expanded.begin(),auxVecs.idPlt_value_expanded.end(),0);
*/



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


	thrust::constant_iterator<unsigned> pltlast = pltfirst + (auxVecs.endIndexBucketPltKeys); // this is NOT numerical addition!

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
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.idPlt_bucket_expanded.begin(),
				pltcountingBegin)),
		NeighborFunctor(
			domainParams.XBucketCount,
			domainParams.YBucketCount,
			domainParams.ZBucketCount));



	unsigned pltnumberOfOutOfRange = thrust::count_if(auxVecs.idPlt_bucket_expanded.begin(),
		auxVecs.idPlt_bucket_expanded.end(), is_greater_than(domainParams.totalBucketCount) );
	unsigned pltnumberInsideRange = endIndexPltExpanded - pltnumberOfOutOfRange;

	//unsigned endIndexPltSearch = endIndexPltExpanded - pltnumberOfOutOfRange;

	thrust::sort_by_key(auxVecs.idPlt_bucket_expanded.begin(),
		auxVecs.idPlt_bucket_expanded.begin() + endIndexPltExpanded,
		auxVecs.idPlt_value_expanded.begin());

	pltnumberInsideRange =
		thrust::get<0>(thrust::unique_by_key(auxVecs.idPlt_value_expanded.begin(),
			auxVecs.idPlt_value_expanded.begin() + endIndexExpanded,
			auxVecs.idPlt_bucket_expanded.begin())) - auxVecs.idPlt_value_expanded.begin();

	auxVecs.idPlt_bucket_expanded.erase(
			auxVecs.idPlt_bucket_expanded.begin() + pltnumberInsideRange,
			auxVecs.idPlt_bucket_expanded.end());

	auxVecs.idPlt_value_expanded.erase(
			auxVecs.idPlt_value_expanded.begin() + pltnumberInsideRange,
			auxVecs.idPlt_value_expanded.end());




	thrust::counting_iterator<unsigned> pltsearch_begin(0);

	thrust::lower_bound(auxVecs.idPlt_bucket_expanded.begin(),
		auxVecs.idPlt_bucket_expanded.end(), pltsearch_begin,
		pltsearch_begin + domainParams.totalBucketCount,
		auxVecs.keyPltBegin.begin());

	thrust::upper_bound(auxVecs.idPlt_bucket_expanded.begin(),
		auxVecs.idPlt_bucket_expanded.end(),pltsearch_begin,
		pltsearch_begin + domainParams.totalBucketCount,
		auxVecs.keyPltEnd.begin());

}


//takes nodes and places in buckets.
void buildBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {


	thrust::counting_iterator<unsigned> indexBucketBegin(0);
	// takes counting iterator and coordinates
	// return tuple of keys and values
	// transform the points to their bucket indices

	//std::cout<<"bucket nodes"<<std::endl;
	thrust::for_each(
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
		BucketIndexer(
			domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.gridSpacing,
			thrust::raw_pointer_cast(auxVecs.id_bucket.data()),
			thrust::raw_pointer_cast(auxVecs.id_value.data())));

//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.id_value.begin(),
		auxVecs.id_value.begin() + generalParams.maxNodeCount,
		auxVecs.id_bucket.begin());
unsigned numberOutOfRange = thrust::count(auxVecs.id_bucket.begin(),
			auxVecs.id_bucket.begin() + generalParams.maxNodeCount, ULONG_MAX);

	auxVecs.endIndexBucketKeys = generalParams.maxNodeCount - numberOutOfRange;

	//platelets
	//std::cout<<"bucket platelet"<<std::endl;
	thrust::counting_iterator<unsigned> indexBucketBegin1(0);
	thrust::for_each(
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
		BucketIndexer(
			domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.gridSpacing,
			thrust::raw_pointer_cast(auxVecs.idPlt_bucket.data()),
			thrust::raw_pointer_cast(auxVecs.idPlt_value.data())));


	//std::cout<<"end bucket platelet"<<std::endl;
//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.idPlt_value.begin(),
		auxVecs.idPlt_value.begin() + generalParams.maxPltCount,
		auxVecs.idPlt_bucket.begin());

unsigned numberPltOutOfRange = thrust::count(auxVecs.idPlt_bucket.begin(),
			auxVecs.idPlt_bucket.begin() + generalParams.maxPltCount, ULONG_MAX);

	auxVecs.endIndexBucketPltKeys = generalParams.maxPltCount - numberPltOutOfRange;


};
