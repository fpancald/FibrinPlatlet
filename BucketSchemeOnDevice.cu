
#include "BucketSchemeOnDevice.h"
#include "NodeSystemDevice.h"

//take domain and discretize into square buckets of size gridspace
void initDimensionBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	domainParams.minX = (*(thrust::min_element(nodeInfoVecs.nodeLocX.begin(), nodeInfoVecs.nodeLocX.end())));
	domainParams.maxX = (*(thrust::max_element(nodeInfoVecs.nodeLocX.begin(), nodeInfoVecs.nodeLocX.end())));
	domainParams.minY = (*(thrust::min_element(nodeInfoVecs.nodeLocY.begin(), nodeInfoVecs.nodeLocY.end())));
	domainParams.maxY = (*(thrust::max_element(nodeInfoVecs.nodeLocY.begin(), nodeInfoVecs.nodeLocY.end())));
	domainParams.minZ = (*(thrust::min_element(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end())));
	domainParams.maxZ = (*(thrust::max_element(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end())));

	//platelets
	domainParams.pltminX = (*(thrust::min_element(pltInfoVecs.pltLocX.begin(), pltInfoVecs.pltLocX.end())));
	domainParams.pltmaxX = (*(thrust::max_element(pltInfoVecs.pltLocX.begin(), pltInfoVecs.pltLocX.end())));
	domainParams.pltminY = (*(thrust::min_element(pltInfoVecs.pltLocY.begin(), pltInfoVecs.pltLocY.end())));
	domainParams.pltmaxY = (*(thrust::max_element(pltInfoVecs.pltLocY.begin(), pltInfoVecs.pltLocY.end())));
	domainParams.pltminZ = (*(thrust::min_element(pltInfoVecs.pltLocZ.begin(), pltInfoVecs.pltLocZ.end())));
	domainParams.pltmaxZ = (*(thrust::max_element(pltInfoVecs.pltLocZ.begin(), pltInfoVecs.pltLocZ.end())));

	if (domainParams.minX > domainParams.pltminX){
		domainParams.minX = domainParams.pltminX;
	}
	if (domainParams.maxX < domainParams.pltmaxX){
		domainParams.maxX = domainParams.pltmaxX;
	}
	if (domainParams.minY > domainParams.pltminY){
		domainParams.minY = domainParams.pltminY;
	}
	if (domainParams.maxY < domainParams.pltmaxY){
		domainParams.maxY = domainParams.pltmaxY;
	}
	if (domainParams.minZ > domainParams.pltminZ){
		domainParams.minZ = domainParams.pltminZ;
	}
	if (domainParams.maxZ < domainParams.pltmaxZ){
		domainParams.maxZ = domainParams.pltmaxZ;
	}


	if (generalParams.iterationCounter == 0) {
		domainParams.originMinX = domainParams.minX;
		domainParams.originMaxX = domainParams.maxX;
		domainParams.originMinY = domainParams.minY;
		domainParams.originMaxY = domainParams.maxY;
		domainParams.originMinZ = domainParams.minZ;
		domainParams.originMaxZ = domainParams.maxZ;

		//double amount of buckets in case of resizing networks
		domainParams.XBucketCount = 2.0 * ceil(ceil(domainParams.maxX - domainParams.minX) / domainParams.gridSpacing) + 1;
		domainParams.YBucketCount = 2.0 * ceil(ceil(domainParams.maxY - domainParams.minY) / domainParams.gridSpacing) + 1;
		domainParams.ZBucketCount = 2.0 * ceil(ceil(domainParams.maxZ - domainParams.minZ) / domainParams.gridSpacing) + 1;
		domainParams.totalBucketCount = domainParams.XBucketCount * domainParams.YBucketCount * domainParams.ZBucketCount;

		if (generalParams.iterationCounter == 0 )
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

}

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
	unsigned valuesCount = auxVecs.bucketValues.size();
	thrust::fill(auxVecs.bucketKeysExpanded.begin(),auxVecs.bucketKeysExpanded.end(),0);
	thrust::fill(auxVecs.bucketValuesIncludingNeighbor.begin(),auxVecs.bucketValuesIncludingNeighbor.end(),0);




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
				auxVecs.bucketKeys.begin(),
				auxVecs.bucketValues.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				auxVecs.bucketValuesIncludingNeighbor.begin())));


	thrust::counting_iterator<unsigned> countingBegin(0);

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				countingBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				countingBegin)) + endIndexExpanded,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				countingBegin)),
		NeighborFunctor(
			domainParams.XBucketCount,
			domainParams.YBucketCount,
			domainParams.ZBucketCount));



	unsigned numberOfOutOfRange = thrust::count_if(auxVecs.bucketKeysExpanded.begin(),
		auxVecs.bucketKeysExpanded.end(), is_greater_than(domainParams.totalBucketCount) );
	unsigned numberInsideRange = endIndexExpanded - numberOfOutOfRange;

	//unsigned endIndexSearch = endIndexExpanded - numberOfOutOfRange;

	thrust::sort_by_key(auxVecs.bucketKeysExpanded.begin(),
		auxVecs.bucketKeysExpanded.begin() + endIndexExpanded,
		auxVecs.bucketValuesIncludingNeighbor.begin());

	auxVecs.bucketKeysExpanded.erase(
			auxVecs.bucketKeysExpanded.begin() + numberInsideRange,
			auxVecs.bucketKeysExpanded.end());

	auxVecs.bucketValuesIncludingNeighbor.erase(
			auxVecs.bucketValuesIncludingNeighbor.begin() + numberInsideRange,
			auxVecs.bucketValuesIncludingNeighbor.end());




	thrust::counting_iterator<unsigned> search_begin(0);

	thrust::lower_bound(auxVecs.bucketKeysExpanded.begin(),
		auxVecs.bucketKeysExpanded.end(), search_begin,
		search_begin + domainParams.totalBucketCount,
		auxVecs.keyBegin.begin());

	thrust::upper_bound(auxVecs.bucketKeysExpanded.begin(),
		auxVecs.bucketKeysExpanded.end(),search_begin,
		search_begin + domainParams.totalBucketCount,
		auxVecs.keyEnd.begin());

	//platelets
	unsigned valuesPltCount = auxVecs.bucketPltValues.size();
	thrust::fill(auxVecs.bucketPltKeysExpanded.begin(),auxVecs.bucketPltKeysExpanded.end(),0);
	thrust::fill(auxVecs.bucketPltValuesIncludingNeighbor.begin(),auxVecs.bucketPltValuesIncludingNeighbor.end(),0);




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
				auxVecs.bucketPltKeys.begin(),
				auxVecs.bucketPltValues.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketPltKeysExpanded.begin(),
				auxVecs.bucketPltValuesIncludingNeighbor.begin())));


	thrust::counting_iterator<unsigned> pltcountingBegin(0);

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketPltKeysExpanded.begin(),
				pltcountingBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketPltKeysExpanded.begin(),
				pltcountingBegin)) + endIndexPltExpanded,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketPltKeysExpanded.begin(),
				pltcountingBegin)),
		NeighborFunctor(
			domainParams.XBucketCount,
			domainParams.YBucketCount,
			domainParams.ZBucketCount));



	unsigned pltnumberOfOutOfRange = thrust::count_if(auxVecs.bucketPltKeysExpanded.begin(),
		auxVecs.bucketPltKeysExpanded.end(), is_greater_than(domainParams.totalBucketCount) );
	unsigned pltnumberInsideRange = endIndexPltExpanded - pltnumberOfOutOfRange;

	//unsigned endIndexPltSearch = endIndexPltExpanded - pltnumberOfOutOfRange;

	thrust::sort_by_key(auxVecs.bucketPltKeysExpanded.begin(),
		auxVecs.bucketPltKeysExpanded.begin() + endIndexPltExpanded,
		auxVecs.bucketPltValuesIncludingNeighbor.begin());

	auxVecs.bucketPltKeysExpanded.erase(
			auxVecs.bucketPltKeysExpanded.begin() + pltnumberInsideRange,
			auxVecs.bucketPltKeysExpanded.end());

	auxVecs.bucketPltValuesIncludingNeighbor.erase(
			auxVecs.bucketPltValuesIncludingNeighbor.begin() + pltnumberInsideRange,
			auxVecs.bucketPltValuesIncludingNeighbor.end());




	thrust::counting_iterator<unsigned> pltsearch_begin(0);

	thrust::lower_bound(auxVecs.bucketPltKeysExpanded.begin(),
		auxVecs.bucketPltKeysExpanded.end(), pltsearch_begin,
		pltsearch_begin + domainParams.totalBucketCount,
		auxVecs.keyPltBegin.begin());

	thrust::upper_bound(auxVecs.bucketPltKeysExpanded.begin(),
		auxVecs.bucketPltKeysExpanded.end(),pltsearch_begin,
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
			thrust::raw_pointer_cast(auxVecs.bucketKeys.data()),
			thrust::raw_pointer_cast(auxVecs.bucketValues.data())));

//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.bucketValues.begin(),
		auxVecs.bucketValues.begin() + generalParams.maxNodeCount,
		auxVecs.bucketKeys.begin());
	unsigned numberOutOfRange = thrust::count(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.begin() + generalParams.maxNodeCount, ULONG_MAX);

	auxVecs.endIndexBucketKeys = generalParams.maxNodeCount - numberOutOfRange;

	//platelets
	thrust::for_each(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				pltInfoVecs.pltLocX.begin(),
				pltInfoVecs.pltLocY.begin(),
				pltInfoVecs.pltLocZ.begin(),
				indexBucketBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				pltInfoVecs.pltLocX.begin(),
				pltInfoVecs.pltLocY.begin(),
				pltInfoVecs.pltLocZ.begin(),
				indexBucketBegin)) + generalParams.maxPltCount,
		BucketIndexer(
			domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.gridSpacing,
			thrust::raw_pointer_cast(auxVecs.bucketPltKeys.data()),
			thrust::raw_pointer_cast(auxVecs.bucketPltValues.data())));

//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.bucketPltValues.begin(),
		auxVecs.bucketPltValues.begin() + generalParams.maxPltCount,
		auxVecs.bucketPltKeys.begin());
	unsigned numberPltOutOfRange = thrust::count(auxVecs.bucketPltKeys.begin(),
			auxVecs.bucketPltKeys.begin() + generalParams.maxPltCount, ULONG_MAX);

	auxVecs.endIndexBucketPltKeys = generalParams.maxPltCount - numberPltOutOfRange;


};
