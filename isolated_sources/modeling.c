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
	int *csxbeg, *cszbeg;
	int *sxz, *gxz;
  	float amp, t0, t1, *wfld ;
	float *wlt, *dobs, **v0;
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
	ns=4, ng=400 ;
    gxbeg=25, gzbeg=4, jsx=80, jsz=20, jgx=1, jgz=0 ;
	amp=1.0, hfop=2, ompchunk=1 ;
	// End

    p0=alloc2d(nz,nx);
    p1=alloc2d(nz,nx);
    p2=alloc2d(nz,nx);
	wlt=(float*)malloc(nt*sizeof(float));
	dobs=(float*)malloc(ng*nt*sizeof(float));
	sxz=(int*)malloc(ns*sizeof(int));
	gxz=(int*)malloc(ng*sizeof(int));
	csxbeg=(int*)malloc(ns*sizeof(int));
	cszbeg=(int*)malloc(ns*sizeof(int));

	v0=alloc2d(nz,nx);
    wfld=(float*)malloc(nt*nz*nx*sizeof(float));
    
	memset(dobs,0,ng*nt*sizeof(float));
	printf("Reading Velocity Model\n");
    readbin(v0[0],"win_overt.bin");
    
	//custom source
    int dsx[4]={170,270,300,220};
    int dsz[4]={55,120,70,65};
    //int dsx[1] = {220};
    //int dsz[1] = {80};

	cszbeg=dsz;
	csxbeg=dsx;

	cst_sg_init(sxz, cszbeg, csxbeg, ns, nz);
	if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
	{ fprintf(stderr,"geophones exceeds the computing zone!\n"); exit(1);}
	sg_init(gxz, gzbeg, gxbeg, jgz, jgx, ng, nz);
    
    memset(wfld,0,nt*nz*nx*sizeof(float));
    
    dtx=dt/dx;
    dtz=dt/dz;
    
	for(is=0; is<ns; is++){

        memset(p0[0],0,nz*nx*sizeof(float));
        memset(p1[0],0,nz*nx*sizeof(float));
        memset(p2[0],0,nz*nx*sizeof(float));

		switch(is){
			case 0 :
				printf("Reading Source 1 \n");
				readbin(wlt,"source1.bin");
				break;
			case 1 :
				readbin(wlt,"source2.bin");
				printf("Reading Source 2 \n");
				break;
			case 2 :
				readbin(wlt,"source3.bin");
				printf("Reading Source 3 \n");
				break;
			case 3 :
				readbin(wlt,"source4.bin");
				printf("Reading Source 4 \n");
				break;
			}
		start=omp_get_wtime();
		/* forward modeling */
        flag=0;
        for(it=0; it<nt; it++){
            add_source(p1, &wlt[it], &sxz[is], 1, nz, true);
            //Recording Data
            record_seis(&dobs[it*ng], gxz, p0, ng, nz);
            //Recording wavefield
            /*
            for (i=0; i<nx; i++) {
                for (j=0; j<nz; j++){
                    wfld[flag] = wfld[flag]+p0[i][j];
                    flag=flag+1;
                }
            }
            */
            step_forward(p0, p1, p2, v0, dtz, dtx, nz, nx);
            ptr=p0; p0=p1; p1=p2; p2=ptr;
			}
			end = omp_get_wtime();
			printf("modelling finished: %f seconds\n", (end-start));
        }
			printf("Writing shots in binary file\n");
			writebin(dobs,ng*nt,"dobs.bin");
           /* writebin(wfld,nt*nz*nx,"fswfld.bin");*/
			
    free(wlt);
    free(*v0); free(v0);
    free(*p0); free(p0);
    free(*p1); free(p1);
    free(*p2); free(p2);
	
    exit(0);
				
	}
