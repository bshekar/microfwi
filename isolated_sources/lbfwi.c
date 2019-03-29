#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include "lbfunc.c"

#ifdef _OPENMP
#include <omp.h>
#include "ompfunc.c"
#endif

int main(int argc, char *argv[])
{
	int it,ix,iz,ind, N;
	int nz, nx, nt, ns, ng, nb, hfop, ompchunk;
	int gxbeg, gzbeg,jgx, jgz ; /* acquisiton parameters */
	int *gxz ;
  	float dx, dz, fm, dt, amp ;
	float *wlt, **vl, **vv, **ptr=NULL;
    	float *dobs;
	float start, end;
	fdm2d fdm;


#ifdef _OPENMP
   omp_init();
#endif
	start=omp_get_wtime();
	// Variable Initialization
	nz=160, nx=450, nt=3000, ng=400 ;
	gxbeg=25, gzbeg=4, jgx=1, jgz=0 ;
	dx=0.025, dz=0.025,fm=10, dt=0.002, amp=1.0, nb=80, hfop=2, ompchunk=1 ;

	// End
	
    // Memory allocation
	
    dobs=(float*)malloc(ng*nt*sizeof(float));   /* observed data */
	fdm=fd_init(nx, nz, nt, dx, dz, dt, nb, hfop, ompchunk);
	wlt=(float*)malloc(nt*sizeof(float));/* source wavelet */
    gxz=(int*)malloc(ng*sizeof(int));
	vl=alloc2d(fdm->nz,fdm->nx);
	vv=alloc2d(fdm->nzpad,fdm->nxpad);

	
	/*printf("Reading Shotgather Data\n");*/
    readbin(dobs,"dobs.bin");
	
	/*printf("Reading Velocity Model\n");*/
	readbin(vl[0],"sm_win_overt.bin");

	if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
	{ fprintf(stderr,"geophones exceeds the computing zone!\n"); exit(1);}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
	   
	/*printf("Starting LBFGS routine.........\n");*/
	/* LBFGS setup */
	int rets;
	N=fdm->nx*fdm->nz*fdm->nt;
	lbfgsfloatval_t fs;
    lbfgsfloatval_t *s = lbfgs_malloc(N);
    lbfgs_parameter_t param_s ;
    
    if (s == NULL ) {
        printf("ERROR: Failed to allocate a memory block for LBFGS variables.\n");
        return 1;
    }
    
    /*Allocate memory for instance data*/
	instanceData *instance_s = malloc(sizeof(instanceData));
    instance_s->vv=vl;
	instance_s->gxz=gxz;
	instance_s->fdm=fdm;
	instance_s->nt=nt;
	instance_s->ng=ng;
	instance_s->dt=dt;
    instance_s->dobs=dobs;
	
    lbfgs_parameter_init(&param_s);
	param_s.orthantwise_c=500.0; /*sparsity promoting parameter */
	param_s.orthantwise_start=1;
	param_s.orthantwise_end=N;
	param_s.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
	param_s.max_iterations=4;
	
   
    //Copy initial model
    	ind = 0;
	for(iz=0; iz<nz; iz++){
        	for(ix=0; ix<nx; ix++){
            		for(it=0; it<nt; it++){
                		ind=it+nt*ix+nt*nx*iz;
                		s[ind]=0.0;
            }
        }
    }
    

    /* Start the L-BFGS optimization; this will invoke the callback functions
     evaluate() and progress() when necessary. */
    
 	rets= lbfgs(N, s, &fs, evaluate_s, progress_s, instance_s, &param_s);
    /* Report the result. */
   	printf("L-BFGS optimization terminated with status code = %d\n fx = %f \n", rets,fs);

    /** Copy source parameters **/
    
	float *src1;
	src1=(float*)malloc(N*sizeof(float));

	for(it=0; it<N; it++){
        	src1[it] = s[it];
    		}
    
    
	lbfgs_free(s);
    	free(dobs); free(gxz);
    	free(*vl); free(vl);
    
    	writebin(src1,N,"inverted_src3d.bin");
   	end=omp_get_wtime(); 
	printf("time taken=%f seconds\n",(end-start));	
    	free(src1);
	
    
	return 0;
    
}
