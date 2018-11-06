#include "NodeSystemDevice.h"
#include "AdvancePositionOnDevice.h"

double AdvancePositionOnDevice(
	NodeInfoVecs& nodeInfoVecs,
  PltInfoVecs& pltInfoVecs,
	GeneralParams& generalParams) {


		//At this point, the previous node location is the same as the current node,
		//we can therefore use previous node locations to update nodeLoc.
		 unsigned _seed = rand();
    	thrust::device_vector<double> gaussianData;
    	gaussianData.resize(generalParams.maxNodeCount); //
		thrust::counting_iterator<unsigned> index_sequence_begin(_seed);

    	thrust::transform(thrust::device, index_sequence_begin, index_sequence_begin + (generalParams.maxNodeCount),
        gaussianData.begin(), psrunifgen(-1.0, 1.0));

		thrust::counting_iterator<unsigned> nodeIndexBegin(0);

		thrust::transform(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeIndexBegin,
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeIndexBegin,
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin())) + generalParams.maxNodeCount,
			//second vector begin
			thrust::make_zip_iterator(
				thrust::make_tuple(
					gaussianData.begin(),
					nodeInfoVecs.nodeForceX.begin(),
					nodeInfoVecs.nodeForceY.begin(),
					nodeInfoVecs.nodeForceZ.begin())),
			//save result in third vector to test values
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin(),
					nodeInfoVecs.nodeVelocity.begin())),
			SaxpyFunctorDim3(generalParams.dtTemp,
				generalParams.viscousDamp,
				generalParams.temperature,
				generalParams.kB,
				generalParams.nodeMass,
				generalParams.maxNodeCount,
				thrust::raw_pointer_cast(nodeInfoVecs.isNodeFixed.data())));

		//finally, clear the random data.
        gaussianData.clear();
        gaussianData.shrink_to_fit();

//platelets
unsigned _seedplt = rand();
 thrust::device_vector<double> gaussianPltData;
 gaussianPltData.resize(generalParams.maxPltCount); //
thrust::counting_iterator<unsigned> index_sequence_plt_begin(_seedplt);

 thrust::transform(thrust::device, index_sequence_plt_begin, index_sequence_plt_begin + (generalParams.maxPltCount),
	 gaussianPltData.begin(), psrunifgen(-1.0, 1.0));

thrust::counting_iterator<unsigned> pltIndexBegin(0);

thrust::transform(
 thrust::make_zip_iterator(
	 thrust::make_tuple(
		 pltIndexBegin,
		 pltInfoVecs.pltLocX.begin(),
		 pltInfoVecs.pltLocY.begin(),
		 pltInfoVecs.pltLocZ.begin())),
 thrust::make_zip_iterator(
	 thrust::make_tuple(
		 pltIndexBegin,
		 pltInfoVecs.pltLocX.begin(),
		 pltInfoVecs.pltLocY.begin(),
		 pltInfoVecs.pltLocZ.begin())) + generalParams.maxPltCount,
 //second vector begin
 thrust::make_zip_iterator(
	 thrust::make_tuple(
		 gaussianPltData.begin(),
		 pltInfoVecs.pltForceX.begin(),
		 pltInfoVecs.pltForceY.begin(),
		 pltInfoVecs.pltForceZ.begin())),
 //save result in third vector to test values
 thrust::make_zip_iterator(
	 thrust::make_tuple(
		 pltInfoVecs.pltLocX.begin(),
		 pltInfoVecs.pltLocY.begin(),
		 pltInfoVecs.pltLocZ.begin(),
		 pltInfoVecs.pltVelocity.begin())),
 SaxpyFunctorDim3(generalParams.dtTemp,
	 generalParams.viscousDamp,
	 generalParams.temperature,
	 generalParams.kB,
	 generalParams.pltMass,
	 generalParams.maxPltCount,
	 thrust::raw_pointer_cast(pltInfoVecs.isPltFixed.data())));

//finally, clear the random data.
	 gaussianPltData.clear();
	 gaussianPltData.shrink_to_fit();

	return generalParams.dtTemp;
		//now that nodeLoc is different, we can calculate change and then set previous location
		//to the current location.

}
