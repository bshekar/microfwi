#ifdef _OPENMP
#include <omp.h>
#endif


omp_init()
{
	
	int threads , cthreads;
	
	threads = omp_get_max_threads();
	omp_set_num_threads(threads);
#pragma omp parallel

	#pragma omp master
	{ cthreads=omp_get_num_threads();
	 printf("Using %d of %d threads \n",cthreads,threads);
	}
	

	}
