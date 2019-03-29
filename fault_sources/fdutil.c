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
    printf("Writing complete \n");
    fclose(file);
}


/*** Finite difference functions ***/

void add_source(float **u, float *source, int *sxz, int ns, int nz, bool add)
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
    
    /*
     // top boundary
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
