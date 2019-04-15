#include <stdlib.h>
#include <float.h>
#include <lbfgs.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.14159265358979323846
#define EPS FLT_EPSILON 

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


/*----------------Data types and intialization----------------*/



typedef struct fdm2 {
	
	int nb;
    int nz,nzpad;
    int nx,nxpad;
    int nt;
    float dz,dx,dt;
    /*bool free;*/
    int hfop,ompchunk;
    
	} *fdm2d ;


fdm2d fd_init( int nx, int nz, int nt, float dx, float dz, float dt, int nb, int hfop, int ompchunk)
/*< Initialize the fdm structure >*/
{	fdm2d fdm ;
	fdm = (fdm2d)malloc(sizeof(*fdm));
	
	fdm->nb=nb;
	
	fdm->nz=nz;
	fdm->nx=nx;
    
    fdm->nt=nt;
	
	fdm->dz=dz;
	fdm->dx=dx;
    
    fdm->dt=dt;
			
	fdm->nzpad=nz+fdm->nb;
	fdm->nxpad=nx+2*fdm->nb;
	

	fdm->hfop = hfop;
	fdm->ompchunk = ompchunk ;
	
	return fdm ;
	}


/*<-----------------General Utilities-------------------------->*/

float **alloc2d(int nz, int nx)
/*< Allocate 2D array in a contiguous memory >*/
{
	int i;
	float **mat = (float**)malloc(nx*sizeof(float*));
	mat[0]= (float *)malloc(nz*nx*sizeof(float *));
	for (i=0;i<nx;i++){
        /*mat[i]=(float*)malloc(nz*sizeof(float));*/
        mat[i]=mat[0]+i*nz;
    }
	return mat;
}

void dealloc2d(float **mat, int nx)
/*< Free 2D array in contiguous memory >*/
{
	int i;
	for (i=0;i<nx;i++){
        free(mat[i]);
    }
    free(mat);
}

float ***alloc3d(int nz, int nx, int nt)
/*< Allocate 3D array in a contiguous memory >*/
{
	int i;
	float ***cell= (float ***)malloc(nt*sizeof(float*));
    cell[0] = alloc2d(nz,nx*nt);
    for (i=0; i< nt; i++) {
		cell[i] = cell[0]+i*nx;
        /*cell[i] = alloc2d(nz,nx);*/
    }
    return cell;
}

void dealloc3d(float ***mat, int nt, int nx)
/*< Free 3D array in contiguous memory >*/
{
	int i;
	for (i=0;i<nt;i++){
		dealloc2d(mat[i],nx);
    }
	free(mat);
}

void readbin(float *data,char *name)
/*< Read Binary File >*/
{
	FILE *file;
	long fsize;
	size_t result;
	file=fopen(name,"rb");
	if (file==NULL) {fputs ("File error \n",stderr); exit(1);}
	// obtain file size
	fseek(file , 0 , SEEK_END);
	fsize = ftell(file);
	rewind(file);
	//Memory check
	if (data==NULL) { fputs("Memory error \n",stderr); exit(2);}
	// copy file into pointers
	result = fread(data,1,fsize,file);
	if (result != fsize){fputs("Reading error \n",stderr); exit(2);}
	/*printf("Loading file complete \n");*/
	fclose(file);	
	
	}
	
	
void writebin(float *data,int size, char *name)
/*< Write Binary File  >*/
{
	FILE *file;
	file=fopen(name,"wb");
	fwrite(data,sizeof(float),size,file);
	/*printf("Writing complete \n");*/
	fclose(file);
	}	
	
/*** Finite difference functions ***/
void add_source(float **u, float *source, int *sxz, int nz, int ns, bool add)
/*< add/subtract seismic sources >*/
{
    int is, sx, sz;
    if(add){
        for(is=0;is<ns; is++){
            sx=sxz[is]/nz;
            sz=sxz[is]%nz;
            u[sx][sz]+=source[is];
        }
    }else{
        for(is=0;is<ns; is++){
            sx=sxz[is]/nz;
            sz=sxz[is]%nz;
            u[sx][sz]-=source[is];
        }
    }
}

void add_spatiotemporal_source(float **p, float *source, int nx, int nz, int nt, int tstep)
/*< add seismic sources >*/
{
    int ix, iz;
    /* #ifdef _OPENMP
     #pragma omp parallel for default(none)        \
     private(ix,iz)                \
     shared(p,source,nz,nx,nt,nb,tstep)
     #endif */
    
    for (ix=0;ix<nx;ix++) {
        for (iz=0;iz<nz;iz++) {
            p[ix][iz]=p[ix][iz]+source[tstep+nt*ix+nt*nx*iz];
        }
    }
}

void record_seis(float *seis_it, int *gxz, float **u, int ng, int nz)
/*< record seismogram at time it into a vector length of ng >*/
{
    int ig, gx, gz;
    for(ig=0;ig<ng; ig++){
        gx=gxz[ig]/nz;
        gz=gxz[ig]%nz;
        seis_it[ig]+=u[gx][gz];
    }
}

void sg_init(int *sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns, int nz)
/*< shot/geophone position initialize >*/
{
    int is, sz, sx;
    for(is=0; is<ns; is++){
        sz=szbeg+is*jsz;
        sx=sxbeg+is*jsx;
        sxz[is]=sz+nz*sx;
    }
}

void cst_sg_init(int *sxz, int* szbeg, int*sxbeg, int ns, int nz)
/*< custom add source>*/
{
    int is, sz, sx;
    for(is=0; is<ns; is++){
        sz=szbeg[is];
        sx=sxbeg[is];
        /*printf("sz and sx is %d and %d \n",sz,sx);*/
        sxz[is]=sz+nz*sx;
        /*printf("sxz is %d\n",sxz[is] );*/
    }
}

void step_forward(float **p0, float **p1, float **p2, float **vv, float dtz, float dtx, int nz, int nx)
/*< forward modeling step, Clayton-Enquist ABC incorporated >*/
{
    int ix,iz;
    float v1,v2,diff1,diff2;
    
    for (ix=0; ix < nx; ix++)
        for (iz=0; iz < nz; iz++)
        {
            v1=vv[ix][iz]*dtz; v1=v1*v1;
            v2=vv[ix][iz]*dtx; v2=v2*v2;
            diff1=diff2=-2.0*p1[ix][iz];
            diff1+=(iz-1>=0)?p1[ix][iz-1]:0.0;
            diff1+=(iz+1<nz)?p1[ix][iz+1]:0.0;
            diff2+=(ix-1>=0)?p1[ix-1][iz]:0.0;
            diff2+=(ix+1<nx)?p1[ix+1][iz]:0.0;
            diff1*=v1;
            diff2*=v2;
            p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
        }
    
    
     // top boundary
    /* 
     iz=0;
     for (ix=1; ix < nx-1; ix++) {
     v1=vv[ix][iz]*dtz;
     v2=vv[ix][iz]*dtx;
     diff1=    (p1[ix][iz+1]-p1[ix][iz])-
     (p0[ix][iz+1]-p0[ix][iz]);
     diff2=    p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
     diff1*=v1;
     diff2*=0.5*v2*v2;
     p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
     }
     */
    /* bottom boundary */
    iz=nz-1;
    for (ix=1; ix < nx-1; ix++) {
        v1=vv[ix][iz]*dtz;
        v2=vv[ix][iz]*dtx;
        diff1=-(p1[ix][iz]-p1[ix][iz-1])+(p0[ix][iz]-p0[ix][iz-1]);
        diff2=p1[ix-1][iz]-2.0*p1[ix][iz]+p1[ix+1][iz];
        diff1*=v1;
        diff2*=0.5*v2*v2;
        p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
    }
    
    /* left boundary */
    ix=0;
    for (iz=1; iz <nz-1; iz++){
        v1=vv[ix][iz]*dtz;
        v2=vv[ix][iz]*dtx;
        diff1=p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
        diff2=(p1[ix+1][iz]-p1[ix][iz])-(p0[ix+1][iz]-p0[ix][iz]);
        diff1*=0.5*v1*v1;
        diff2*=v2;
        p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
    }
    
    /* right boundary */
    ix=nx-1;
    for (iz=1; iz <nz-1; iz++){
        v1=vv[ix][iz]*dtz;
        v2=vv[ix][iz]*dtx;
        diff1=p1[ix][iz-1]-2.0*p1[ix][iz]+p1[ix][iz+1];
        diff2=-(p1[ix][iz]-p1[ix-1][iz])+(p0[ix][iz]-p0[ix-1][iz]);
        diff1*=0.5*v1*v1;
        diff2*=v2;
        p2[ix][iz]=2.0*p1[ix][iz]-p0[ix][iz]+diff1+diff2;
    }
}

/*----------------FWI routines--------------------*/


void max_energy( float **u, int *mx, int *mz, float *e, fdm2d fdm)
/*< maximum energy condition >*/
{
    int ix,iz;
    float e1, eloc1, eloc2=0;
    for(ix=0;ix<fdm->nx;ix++){
        for(iz=0;iz<fdm->nz;iz++){
            e1+=u[ix+fdm->nb][iz]*u[ix+fdm->nb][iz];
            eloc1=u[ix+fdm->nb][iz]*u[ix+fdm->nb][iz];
            if(eloc1>eloc2){eloc2=eloc1;mx[0]=ix;mz[0]=iz;}
        }
    }
    if(e1>= *e){*e=e1;}
}


void cal_residuals(float *dcal, float *dobs, float *derr, int ng)
/*< calculate residual wavefield at the receiver positions
 dcal: d_{cal}
 dobs: d_{obs}
 derr: d_{err}=d_{cal}-d_{obs} >*/
{
	int ig;
	for(ig=0;ig<ng;ig++){
        derr[ig]=dcal[ig]-dobs[ig];
    }
}

void cal_objective(float *obj,float *err, int ng)
/*< calculate the L2 norm objective function >*/
{
	int ig;
	float result=0.0;
#ifdef _OPENMP
#pragma omp parallel for shared(err)  	\
private(ig) reduction(+:result)
#endif
	for(ig=0;ig<ng;ig++){
		result+=(err[ig]*err[ig]);
    }
	*obj=result ;
	
}


/*<---------------Smoothing routines--------------------------->*/

/*--- 2D Gaussian function ---*/

void gauss2dfunc(float *ker, int oplx, int oply, float sigx, float sigy, fdm2d fdm)
/*** Populate a 2D Gaussian function ***/
{
    int hoplx, hoply, ix, iy;
    hoplx = (oplx+1)/2; hoply = (oply+1)/2;
    float xsq, ysq, dx, dy;
    dx = fdm->dx; dy = fdm->dz;
    for (ix = 0; ix < oplx; ix++){
        for (iy = 0; iy < oply; iy++){
            xsq = (ix-hoplx)*(ix-hoplx)*dx*dx/(sigx*sigx);
            ysq = (iy-hoply)*(iy-hoply)*dy*dy/(sigy*sigy);
            ker[ix+iy*oplx]=exp(-(xsq+ysq));
        }
    }
}

/* ----- 2D convolution with Gaussian ----- */

void conv2d(float **conv, float **data, int oplx, int oply, float sigx, float sigy, fdm2d fdm)
/*** Implements 2D convolution for the purpose of smoothing with a Gaussian ***/
{
    int hoplx, hoply, startx, starty, endx, endy, ix, iy, i, j, k, l, nx, ny;
    float dumr, *ker;
    
    ker=(float*)malloc(oplx*oply*sizeof(float));
    
    /*** retreive Gaussian operator ***/
    
    gauss2dfunc(ker, oplx, oply, sigx, sigy, fdm);
    
    hoplx = (oplx+1)/2; hoply = (oply+1)/2; /** half operator lengths **/
    
    nx = fdm->nx;
    ny = fdm->nz;
    
    for (iy = 0; iy < ny; iy++){
        
        starty = MAX(iy-hoply+1, 0);
        endy = MIN(iy+hoply, ny);
        
        for (ix = 0; ix < nx; ix++) {
            
            startx = MAX(ix-hoplx+1, 0);
            endx = MIN(ix+hoplx, nx);
            
            
            /* convolution with the operator */
            dumr = 0.0;
            k = MAX(hoply-1-iy, 0);
            for (i = starty; i < endy; i++) {
                l = MAX(hoplx-1-ix, 0);
                for (j = startx; j < endx; j++) {
                    dumr += data[i][j]*ker[k*oplx+l];
                    l++;
                }
                k++;
            }
            conv[ix][iy] = dumr;
        }
    }
}

/*--- 1D Gaussian function ---*/

void gauss1dfunc(float *ker, int opl, float sig, fdm2d fdm)
/*** Populate a 1D Gaussian function ***/
{
    int hopl,it;
    hopl = (opl+1)/2;
    float tsq, dt;
    dt = fdm->dt;
    for (it = 0; it < opl; it++){
        tsq = (it-hopl)*(it-hopl)*dt*dt/(sig*sig);
        ker[it]=exp(-tsq);
    }
}

/* ----- 1D convolution along time axis ----- */
void conv1d(float *conv, float *data, int opl, float sig, fdm2d fdm)
/*** Implements 1D convolution for the purpose of smoothing with a Gaussian ***/
{
    int hopl, starti, endi, ind, ix, iz, it, i, j, k, l, nx, nz, nt;
    float dumr, *ker;
    
    ker=(float*)malloc(opl*sizeof(float));
    
    /*** retreive Gaussian operator ***/
    
    gauss1dfunc(ker, opl, sig, fdm);
    
    hopl = (opl+1)/2;  /** half operator lengths **/
    
    nx = fdm->nx;
    nz = fdm->nz;
    nt = fdm->nt;
    
    for (it = 0; it < nt; it++) {
        starti = MAX(it-hopl+1, 0);
        endi = MIN(it+hopl, nt);
        
        /* convolution with the operator */
        
        ind=0;
        for (ix=0; ix<nx; ix++){
            for (iz=0; iz<nz; iz++){
                dumr = 0.0;
                l = MAX(hopl-1-it, 0);
                for (i=starti; i<endi; i++){
                    ind=i+nt*ix+nt*nx*iz;
                    dumr += data[ind]*ker[l];
                    l++;
                }
                ind = it+nt*ix+nt*nx*iz;
                conv[ind] = dumr;
            }
        }
    }
    
}


/*<----------------Scaling and Filtering-----------------> */

void src_spat_gradient(float **fxz, float **u, int it, fdm2d fdm)
/* source spatial gradient */
{
	int ix,iz;
	for(ix=0;ix<fdm->nx;ix++){
		for(iz=0;iz<fdm->nz;iz++){
			fxz[ix][iz]=u[ix+fdm->nb][iz];
        }
    }
}

void scale_spat_gradient(float **fxz, fdm2d fdm)
/* scale source gradient */
{
	int ix,iz;
	float a,fmax=0;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
private(ix,iz,a) shared(fxz,fdm) reduction(max:fmax)
#endif
	for(ix=0;ix<fdm->nx;ix++){
		for(iz=0;iz<fdm->nz;iz++){
			a=fabsf(fxz[ix][iz]);
			fmax=MAX(a,fmax);
        }
    }
    
	//printf("fmax is %f \n",fmax);
	for(ix=0;ix<fdm->nx;ix++){
		for(iz=0;iz<fdm->nz;iz++){
			fxz[ix][iz]/=(fmax+EPS);
        }
    }
	
}

void scale_3D_gradient(float *gs, int N)
/*  source gradient */
{
	int ix;
	float a,gsmax=0;
#ifdef _OPENMP
#pragma omp parallel for default(none) \
private(ix,a) shared(gs,N) reduction(max:gsmax)
#endif
	for(ix=0;ix<N;ix++){
        a=fabsf(gs[ix]);
        gsmax=MAX(a,gsmax);
    }
    
	/* printf("gsmax is %f \n",gsmax); */ 
	for(ix=0;ix<N;ix++){
        gs[ix]/=(gsmax+EPS);
    }
	
}

void scalegrad(float *gs, fdm2d fdm, int oplx, int oply, float sigx, float sigy)
/* compute a 2D function to do illumination compensation */
{
    int ind,it,ix,iz,nx,nz,nt;
    nx = fdm->nx;
    nz = fdm->nz;
    nt = fdm->nt;
    
    float **gradscale, **gradsmooth;
    gradscale=alloc2d(nz,nx);
    gradsmooth=alloc2d(nz,nx);
    
    for(ix = 0; ix < nx; ix++){
        for(iz=0; iz<nz; iz++){
            gradscale[ix][iz] = 0.0;
            for(it=0; it<nt; it++){
                ind=it+nt*ix+nt*nx*iz;
                gradscale[ix][iz] += gs[ind]*gs[ind];
            }
        }
    }
    
    /* Smooth out scaling function using 2D Gaussian */
    conv2d(gradsmooth, gradscale, oplx, oply, sigx, sigy, fdm);
    
    /* Now scale the gradient */
    for(ix = 0; ix < nx; ix++){
        for(iz=0; iz<nz; iz++){
            for(it=0; it<nt; it++){
                ind=it+nt*ix+nt*nx*iz;
                gs[ind]=gs[ind]/gradsmooth[ix][iz];
            }
        }
    }
}





/*------------------ LBFGS routines --------------------------*/

typedef struct instanceStruct {
	/*float ***mod; /*source spatial and temporal position  */
	/*float *wlt;   /*wavelet */
	float **vv;  /*velocity */
	int *gxz;   /*receiver location */
	fdm2d fdm;    /*fdm structure */
	int nt;       /*time step*/
	int ng;       /*no of receivers*/
	float dt;	  /* sampling time*/
    float *dobs;  /* Observed data */
	
} instanceData;


instanceData instance_init(float **vv, int *gxz, fdm2d fdm, int nt, int ng, int dt,float *dobs)
/* initialize instance structure */
{
	instanceData *instance;
	instance=(instanceData*)malloc(sizeof(instanceData));
	
	/*instance->mod=mod;*/
	instance->vv=vv;
	/*instance->wlt=wlt;*/
	instance->gxz=gxz;
	instance->fdm=fdm;
	instance->nt=nt;
	instance->ng=ng;
	instance->dt=dt;
	instance->dobs=dobs;
	/* printf("structure successfully allocated \n"); */
	
	return *instance;
	
}







static lbfgsfloatval_t evaluate_s( instanceData *instance,    /** user data sent by lbfgs to client **/
                                  const lbfgsfloatval_t *s,  /** current values of variables **/
                                  lbfgsfloatval_t *gs,        /** gradient vector **/
                                  const int n,               /** number of variables **/
                                  const lbfgsfloatval_t step /** current step of line search routine **/)
/*< source gradient >*/
{
	
	int it,ix,iz,ind,nz,nx,ng,nt;
	float **u0=NULL, **u1=NULL, **u2=NULL, **gp0=NULL, **gp1=NULL, **gp2=NULL, **vv=NULL, **ptr=NULL;
	float   *dobs=NULL, *dcal=NULL, *derr=NULL, *src;
	float obj, dt, dx, dz, dtz, dtx;
	int *gxz=NULL;
	lbfgsfloatval_t fs = 0.0;
    
    //Read model parameters from instance
	
	nx=instance->fdm->nx;
	nz=instance->fdm->nz;
	ng=instance->ng;
	nt=instance->nt;
	dt=instance->dt;
    dx=instance->fdm->dx;
    dz=instance->fdm->dz;
    
	//Memory Allocation: 1D

	dcal=(float*)malloc(ng*nt*sizeof(float));/* calculated data */
	derr=(float*)malloc(ng*nt*sizeof(float));//residuals
    src=(float*)malloc(nz*nx*nt*sizeof(float));//residuals
    
    //Memory Allocation: 2D & 3D
    
    //wavefield
    u0=alloc2d(nz,nx);
	u1=alloc2d(nz,nx);
    u2=alloc2d(nz,nx);
    //adjoint wavefield
	gp0=alloc2d(nz,nx);
	gp1=alloc2d(nz,nx);
    gp2=alloc2d(nz,nx);
    
    memset(u0[0],0,nz*nx*sizeof(float));
	memset(u1[0],0,nz*nx*sizeof(float));
    memset(u2[0],0,nz*nx*sizeof(float));
    
    memset(gp0[0],0,nz*nx*sizeof(float));
    memset(gp1[0],0,nz*nx*sizeof(float));
    memset(gp2[0],0,nz*nx*sizeof(float));
  
    //Read data, velocity model and receiver positions from instance
	
    dobs=instance->dobs;
    vv=instance->vv;
    gxz=instance->gxz;
  
	float smax;
	smax=0.0;	
    //Copy source model and cast into float from lbfgs_float
    for (it=0; it<nx*nz*nt; it++){
        src[it] = (float) s[it];
	/*smax=MAX(fabs(src[it]),smax);*/
	}
	/*printf("smax=%f\n",smax);*/
 
    dtx=dt/dx;
    dtz=dt/dz;
    
	for(it=0; it<nt; it++){
		add_spatiotemporal_source(u1,src,nx,nz,nt,it);
        step_forward(u0, u1, u2, vv, dtz, dtx, nz, nx);
        record_seis(&dcal[it*ng], gxz, u0, ng, nz);
        ptr=u0; u0=u1; u1=u2; u2=ptr;
		cal_residuals(&dcal[it*ng],&dobs[it*ng],&derr[it*ng],ng); //residuals
    }
   /* 
	printf("Writing residual wavefield \n");
	writebin(derr,ng*nt,"resd.bin");*/
    
    
    float *gss,*gssmooth;
    gss=(float*)malloc(nx*nz*nt*sizeof(float));
    gssmooth=(float*)malloc(nz*nx*nt*sizeof(float));
    
    ind = 0;
    
    int oplx, oply;
    float sigx, sigy;
 
    oplx=10; oply=10;
    sigx=1.0*instance->fdm->dx; sigy=1.0*instance->fdm->dz;
	
	for(it=nt-1; it>-1; it--){
        
		/* extrapolate residual wavefield */
		add_source(gp1, &derr[it*ng],gxz,nz,ng,true);
        	step_forward(gp0,gp1,gp2,vv,dtz, dtx, nz, nx);
        	for(iz=0;iz<nz;iz++){
            		for(ix=0; ix<nx; ix++){
                	ind=it+nt*ix+nt*nx*iz;
                	gss[ind] = gp0[ix][iz];
            		}
        	}
        	ptr=gp0; gp0=gp1; gp1=gp2; gp2=ptr;
    	}
    
    /** Scale and smooth gradient in time-space dimensions **/
   /* 
    int opl=10;
    float sig=2.0*instance->fdm->dt;
    conv1d(gssmooth, gss, opl, sig, instance->fdm);
    scale_3D_gradient(gssmooth, nx*nz*nt);
    */
  	/*scalegrad(gss, instance->fdm, oplx, oply, sigx, sigy);*/ 
    
    // copy gradient, this is to make sure that float value types are same
    
    ind = 0;
	for(iz=0;iz<nz;iz++){
	        for(ix=0; ix<nx; ix++){
		    for(it=0;it<nt;it++){
		        ind=it+nt*ix+nt*nx*iz;
		        gs[ind] = gss[ind];
			}
		}
	}
   /* 
     printf("Writing gradient  \n");
     writebin(gss,nz*nx*nt,"gss.bin");
     */
    
	cal_objective(&obj,derr,ng*nt);
	/*printf("Objective Function Value is %f \n",obj);*/
	fs=obj;
    
    free(gss); free(gssmooth);
    free(src);
    free(*u0); free(u0);
    free(*u1); free(u1);
    free(*u2); free(u2);
    free(*gp0); free(gp0);
    free(*gp1); free(gp1);
    free(*gp2); free(gp2);
    
	return fs;
}



static int progress_s(   /*** progress of the optimization process **/
                      instanceData *instance,
                      const lbfgsfloatval_t *s,
                      const lbfgsfloatval_t *gs,
                      const lbfgsfloatval_t fs,
                      const lbfgsfloatval_t snorm,
                      const lbfgsfloatval_t gsnorm,
                      const lbfgsfloatval_t step,
                      int n,
                      int k,
                      int ls
                      )
{
    /*
    printf("Iteration %d:\n", k);
    printf("  fx = %f \n", fs);
    printf("  xnorm = %f, gnorm = %f, step = %f \n", snorm, gsnorm, step);
     */
    int it;
    float sl1norm;
    sl1norm=0.0;
    for(it=0; it<n; it++){
	sl1norm = sl1norm + fabs(s[it]);
    }
    printf(" %d", k);
    printf(" %f ", fs);
    printf(" %f %f %f %f \n", snorm, sl1norm, gsnorm, step);
    return 0;
}
