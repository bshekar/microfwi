#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "fdutil.c"

#ifdef _OPENMP
#include <omp.h>
#include "ompfunc.c"
#endif

#define PI 3.14159265358979323846
#define EPS FLT_EPSILON
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))



int main(int argc, char *argv[])
{
	
	int is, it, i2, i1, i, j, flag;
	int ns, ng;
	int gxbeg, gzbeg, jsx, jsz, jgx, jgz ; /* acquisiton parameters */

	int *sxz, *gxz;
  	float amp, t0, t1;
	float  *dobs, *sxzt, **v0;
    int nz, nx, nt;
    float dz, dx, dt, dtz, dtx, hfop,ompchunk;
	float start, end;
    float **p0, **p1, **p2, **ptr=NULL;

#ifdef _OPENMP
   omp_init();
#endif

	// Variable Initialization
    nz=160, nx=450, nt=3000;
    dx=0.025, dz=0.025, dt=0.002;
	 ng=400 ;
    gxbeg=25, gzbeg=4, jsx=80, jsz=20, jgx=1, jgz=0 ;
	amp=1.0, hfop=2, ompchunk=1 ;
	// End

    p0=alloc2d(nz,nx);
    p1=alloc2d(nz,nx);
    p2=alloc2d(nz,nx);
	
	dobs=(float*)malloc(ng*nt*sizeof(float));

	gxz=(int*)malloc(ng*sizeof(int));

	    sxzt=(float*)malloc(nz*nx*nt*sizeof(float));
		//read in the spatio-temporal source function
    readbin(sxzt,"source3d.bin");


	v0=alloc2d(nz,nx);
    
	memset(dobs,0,ng*nt*sizeof(float));
	printf("Reading Velocity Model\n");
    readbin(v0[0],"win_overt.bin");
    


	if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
	{ fprintf(stderr,"geophones exceeds the computing zone!\n"); exit(1);}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
    
    
    dtx=dt/dx;
    dtz=dt/dz;
    


        memset(p0[0],0,nz*nx*sizeof(float));
        memset(p1[0],0,nz*nx*sizeof(float));
        memset(p2[0],0,nz*nx*sizeof(float));
        

		start=omp_get_wtime();
		/* forward modeling */
        flag=0;
        for(it=0; it<nt; it++){
            add_spatiotemporal_source(p1, sxzt, nx, nz, nt, it);
            //Recording Data
            record_seis(&dobs[it*ng], gxz, p0, ng, nz);

            step_forward(p0, p1, p2, v0, dtz, dtx, nz, nx);
            ptr=p0; p0=p1; p1=p2; p2=ptr;
			}
			end = omp_get_wtime();
			printf("modelling finished: %f seconds\n", (end-start));
        
			printf("Writing shots in binary file\n");
			writebin(dobs,ng*nt,"dobs.bin");
			
    
    free(*v0); free(v0);
    free(*p0); free(p0);
    free(*p1); free(p1);
    free(*p2); free(p2);
	
    exit(0);
				
	}
