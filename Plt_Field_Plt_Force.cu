#include "Plt_Field_Plt_Force.h"
#include "functor_plt_field_plt.h"
#include "System.h"

//for a given platelet, apply force from other platelets
//Applies force to self. 
void Plt_Field_Plt_Force(
  	GeneralParams& generalParams,
  	PltInfoVecs& pltInfoVecs,
  	AuxVecs& auxVecs) {


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
    thrust::make_zip_iterator(
      	thrust::make_tuple(
				//DOES NOT RESET FORCES
        		pltInfoVecs.pltForceX.begin(),
        		pltInfoVecs.pltForceY.begin(),
        		pltInfoVecs.pltForceZ.begin())),

         functor_plt_field_plt(
             generalParams.plt_other_intrct,
             generalParams.pltRForce,
             generalParams.pltForce,
             generalParams.pltR,
             generalParams.maxPltCount,

             thrust::raw_pointer_cast(pltInfoVecs.pltLocX.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltLocY.data()),
             thrust::raw_pointer_cast(pltInfoVecs.pltLocZ.data()),

            // thrust::raw_pointer_cast(pltInfoVecs.pltForceX.data()),
            // thrust::raw_pointer_cast(pltInfoVecs.pltForceY.data()),
            // thrust::raw_pointer_cast(pltInfoVecs.pltForceZ.data()),

             thrust::raw_pointer_cast(auxVecs.idPlt_value_expanded.data()),//plt neighbors
             thrust::raw_pointer_cast(auxVecs.keyPltBegin.data()),
             thrust::raw_pointer_cast(auxVecs.keyPltEnd.data()) ) );
};