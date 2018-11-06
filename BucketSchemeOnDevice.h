#ifndef BUCKETSCHEMEONDEVICE_H_
#define BUCKETSCHEMEONDEVICE_H_

#include "SystemStructures.h"


void initDimensionBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);

void buildBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams);

void extendBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs);




/*
* Functor to compute neighbor buckets(center bucket included) of a node.
* @param input1 bucket index of node
* @param input2 pick from the sequence, which is also global rank of the node
*
* @return output1 bucket indices of node ( all neighbors and the center bucket) of node
* @return output2 global rank of node

example with 3x3x3 grid. with node at position (1,1,1). Notice there are 27 possibilities.
0: (1,1,1)
1: (0,0,1)
2: (1,0,1)
The words left & right denote x change, top and bottom denote y change and upper & lower denote z change
 */

struct NeighborFunctor : public thrust::unary_function<Tuu, Tuu> {
	unsigned numOfBucketsInXDim;
	unsigned numOfBucketsInYDim;
	unsigned numOfBucketsInZDim;

	__host__ __device__ NeighborFunctor(
		unsigned _numOfBucketsInXDim,
		unsigned _numOfBucketsInYDim,
		unsigned _numOfBucketsInZDim ) :
		numOfBucketsInXDim(_numOfBucketsInXDim),
		numOfBucketsInYDim(_numOfBucketsInYDim),
		numOfBucketsInZDim(_numOfBucketsInZDim) {}

	__device__ Tuu operator()(const Tuu &v) {
		unsigned relativeRank = thrust::get<1>(v) % 27;	//27 = 3^3. Takes global node id for calculation

		//takes global bucket id for calculation

		unsigned area = numOfBucketsInXDim * numOfBucketsInYDim;
		unsigned volume = area * numOfBucketsInZDim;

		unsigned xPos = thrust::get<0>(v) % numOfBucketsInXDim;	//col
		unsigned xPosLeft = xPos - 1;
		unsigned xPosRight = xPos + 1;
		if (xPosLeft > volume) {
			//wraparound unsigned
			xPosLeft = numOfBucketsInXDim - 1;
		}
		if (xPosRight >= numOfBucketsInXDim) {
			xPosRight = 0;
		}


		unsigned zPos = thrust::get<0>(v) / area; //z divide by area
		unsigned zPosUp = zPos + 1;
		unsigned zPosLow = zPos - 1;
		if (zPosLow > volume ) {
			//wraparound unsigned
			zPosLow = numOfBucketsInZDim - 1;
		}
		if (zPosUp >= numOfBucketsInZDim) {
			zPosUp = 0;
		}

		unsigned yPos = (thrust::get<0>(v) - zPos * area) / numOfBucketsInXDim;	//row
		unsigned yPosTop = yPos + 1;
		unsigned yPosBottom = yPos - 1;

		if (yPosBottom > volume) {
			//wraparound unsigend
			yPosBottom = numOfBucketsInYDim - 1;
		}
		if (yPosTop >= numOfBucketsInYDim) {
			yPosTop = 0;
		}

		switch (relativeRank) {
		//middle cases
		case 0:
			return thrust::make_tuple(thrust::get<0>(v), thrust::get<1>(v));
			//break;
		case 1:{
				unsigned topLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPos * area;
			return thrust::make_tuple(topLeft, thrust::get<1>(v));
			//break;
		}
		case 2:{
				unsigned top = xPos + yPosTop * numOfBucketsInXDim + zPos * area;
			return thrust::make_tuple(top, thrust::get<1>(v));
			//break;
		}
		case 3:{
				unsigned topRight = xPosRight + yPosTop * numOfBucketsInXDim + zPos * area;
			return thrust::make_tuple(topRight, thrust::get<1>(v));
			//break;
		}
		case 4:{
				unsigned right = xPosRight + yPos * numOfBucketsInXDim + zPos * area;
			return thrust::make_tuple(right, thrust::get<1>(v));
			//break;
		}
		case 5:{
				unsigned bottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPos * area;
			return thrust::make_tuple(bottomRight, thrust::get<1>(v));
			//break;
		}
		case 6:{
				unsigned bottom = xPos + yPosBottom * numOfBucketsInXDim + zPos * area;
			return thrust::make_tuple(bottom, thrust::get<1>(v));
			//break;
		}
		case 7:{
				unsigned bottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim + zPos * area;
			return thrust::make_tuple(bottomLeft, thrust::get<1>(v));
			//break;
		}
		case 8:{
				unsigned left = xPosLeft + yPos * numOfBucketsInXDim + zPos * area;
			return thrust::make_tuple(left, thrust::get<1>(v));
			//break;
		}
		//lower Z cases
		case 9:{
				unsigned lowerCenter = xPos + yPos * numOfBucketsInXDim +  zPosLow * area;
			return thrust::make_tuple(lowerCenter, thrust::get<1>(v));
			//break;
		}
		case 10:{
				unsigned lowerTopLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPosLow* area;
			return thrust::make_tuple(lowerTopLeft, thrust::get<1>(v));
			//break;
		}
		case 11:{
				unsigned lowerTop = xPos + yPosTop * numOfBucketsInXDim + zPosLow * area;
			return thrust::make_tuple(lowerTop, thrust::get<1>(v));
			//break;
		}
		case 12:{
				unsigned lowerTopRight = xPosRight + yPosTop * numOfBucketsInXDim  + zPosLow * area;
			return thrust::make_tuple(lowerTopRight, thrust::get<1>(v));
			//break;
		}
		case 13:{
				unsigned lowerRight = xPosRight + yPos * numOfBucketsInXDim + zPosLow * area;
			return thrust::make_tuple(lowerRight, thrust::get<1>(v));
			//break;
		}
		case 14:{
				unsigned lowerBottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPosLow * area;
			return thrust::make_tuple(lowerBottomRight, thrust::get<1>(v));
			//break;
		}
		case 15:{
				unsigned lowerBottom = xPos + yPosBottom * numOfBucketsInXDim + zPosLow * area;
			return thrust::make_tuple(lowerBottom, thrust::get<1>(v));
			//break;
		}
		case 16:{
				unsigned lowerBottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim  + zPosLow * area;
			return thrust::make_tuple(lowerBottomLeft, thrust::get<1>(v));
			//break;
		}
		case 17:{
				unsigned lowerLeft = xPosLeft + yPos * numOfBucketsInXDim + zPosLow * area;
			return thrust::make_tuple(lowerLeft, thrust::get<1>(v));
			//break;
		}
		//upper Z cases
		case 18:{
				unsigned upperCenter = xPos + yPos * numOfBucketsInXDim +  zPosUp * area;
			return thrust::make_tuple(upperCenter, thrust::get<1>(v));
			//break;
		}
		case 19:{
				unsigned upperTopLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPosUp * area;
			return thrust::make_tuple(upperTopLeft, thrust::get<1>(v));
			//break;
		}
		case 20:{
				unsigned upperTop = xPos + yPosTop * numOfBucketsInXDim + zPosUp * area;
			return thrust::make_tuple(upperTop, thrust::get<1>(v));
			//break;
		}
		case 21:{
				unsigned upperTopRight = xPosRight + yPosTop * numOfBucketsInXDim  + zPosUp * area;
			return thrust::make_tuple(upperTopRight, thrust::get<1>(v));
			//break;
		}
		case 22:{
				unsigned upperRight = xPosRight + yPos * numOfBucketsInXDim + zPosUp * area;
			return thrust::make_tuple(upperRight, thrust::get<1>(v));
			//break;
		}
		case 23:{
				unsigned upperBottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPosUp * area;
			return thrust::make_tuple(upperBottomRight, thrust::get<1>(v));
			//break;
		}
		case 24:{
				unsigned upperBottom = xPos + yPosBottom * numOfBucketsInXDim + zPosUp * area;
			return thrust::make_tuple(upperBottom, thrust::get<1>(v));
			//break;
		}
		case 25:{
				unsigned upperBottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim  + zPosUp * area;
			return thrust::make_tuple(upperBottomLeft, thrust::get<1>(v));
			//break;
		}
		case 26:{
				unsigned upperLeft = xPosLeft + yPos * numOfBucketsInXDim + zPosUp * area;
			return thrust::make_tuple(upperLeft, thrust::get<1>(v));
			//break;
		}
		default:
			return thrust::make_tuple(ULONG_MAX, ULONG_MAX);

		}
	}
};


struct BucketIndexer {
	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;

	double unitLen;
	unsigned XSize;
	unsigned YSize;
	unsigned* bucketKeysAddr;
	unsigned* bucketValuesAddr;

	__host__ __device__

	BucketIndexer(
		double _minX,
		double _maxX,
		double _minY,
		double _maxY,
		double _minZ,
		double _maxZ,
		double _unitLen,
		unsigned* _bucketKeysAddr,
		unsigned* _bucketValuesAddr) :
		minX(_minX),
		maxX(_maxX),
		minY(_minY),
		maxY(_maxY),
		minZ(_minZ),
		maxZ(_maxZ),
		bucketKeysAddr(_bucketKeysAddr),
		bucketValuesAddr(_bucketValuesAddr),
		unitLen(_unitLen),
		XSize((_maxX - _minX) / unitLen),
		YSize((_maxY - _minY) / unitLen) {}

	__device__ void operator()(const Tdddu& v) {
			//static cast double into unsigned
			unsigned id = thrust::get<3>(v);
			unsigned x = static_cast<unsigned>((thrust::get<0>(v) - minX) / unitLen);
			unsigned y = static_cast<unsigned>((thrust::get<1>(v) - minY) / unitLen);
			unsigned z = static_cast<unsigned>((thrust::get<2>(v) - minZ) / unitLen);


			// return the bucket's linear index and node's global index
			//return thrust::make_tuple(z * XSize * YSize + y * XSize + x, thrust::get<4>(v));
			unsigned bucket = z * XSize * YSize + y * XSize + x;
			//try to make it so bucket does not return unsigned32Max
			if (bucket == ULONG_MAX) {
				bucket = 0;
			}
			bucketKeysAddr[id] = bucket;
			bucketValuesAddr[id] = id;

	}
};

template <typename InputIterator1,
          typename InputIterator2,
          typename OutputIterator>
OutputIterator expand(InputIterator1 first1,
                      InputIterator1 last1,
                      InputIterator2 first2,
                      OutputIterator output)
{
  typedef typename thrust::iterator_difference<InputIterator1>::type difference_type;

  difference_type input_size  = thrust::distance(first1, last1);
  difference_type output_size = thrust::reduce(first1, last1);

  // scan the counts to obtain output offsets for each input element
  thrust::device_vector<difference_type> output_offsets(input_size, 0);
  thrust::exclusive_scan(first1, last1, output_offsets.begin());

  // scatter the nonzero counts into their corresponding output positions
  thrust::device_vector<difference_type> output_indices(output_size, 0);
  thrust::scatter_if
    (thrust::counting_iterator<difference_type>(0),
     thrust::counting_iterator<difference_type>(input_size),
     output_offsets.begin(),
     first1,
     output_indices.begin());

  // compute max-scan over the output indices, filling in the holes
  thrust::inclusive_scan
    (output_indices.begin(),
     output_indices.end(),
     output_indices.begin(),
     thrust::maximum<difference_type>());

  // gather input values according to index array (output = first2[output_indices])
  OutputIterator output_end = output; thrust::advance(output_end, output_size);
  thrust::gather(output_indices.begin(),
                 output_indices.end(),
                 first2,
                 output);

  // return output + output_size
  thrust::advance(output, output_size);
  return output;
}





#endif /* BUCKETSCHEMEONDEVICE_H_*/
