#ifndef FUNCTOR_NORM_H_
#define FUNCTOR_NORM_H_

#include "SystemStructures.h"

struct functor_norm {

		__host__ __device__
		double operator() (const CVec3& vec) {
		//divide force by fiber cross section to get stress
		double result = sqrt(thrust::get<0>(vec) * thrust::get<0>(vec) +
			thrust::get<1>(vec) * thrust::get<1>(vec) +
			thrust::get<2>(vec) * thrust::get<2>(vec));

		return result;


	}
};

#endif