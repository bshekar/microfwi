import numpy as np
from scipy import ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter

ng=450;nt=3000;nb=80;nz=160;nx=450;
ns=2

#Reading Binary File

f=open('inverted_src3d.bin','rb')
src=np.fromfile(f,dtype='float32')
src=src.reshape((nt,nx,nz),order='f')
src=gaussian_filter(src,sigma=[5,0,0]);

#Spatial source s(x,z)
sxz=np.sum(np.square(src),axis=0).T;
sxz=np.sqrt(sxz);
sxz=gaussian_filter(sxz,sigma=3.0);
perclip=np.percentile(sxz,99.0)
blobs=sxz>perclip ### connected region
labels, nlabels = ndimage.label(blobs)


#source model

fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)
im=ax.imshow(sxz,extent=[0,11.25,4.0,0])
cbar=fig.colorbar(im, cax=cax, orientation='vertical')
cbar.ax.tick_params(labelsize=8)
ax.set_xlabel("x (km)")
ax.set_ylabel("z (km)")
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
fig.savefig("source_power.pdf", bbox_inches='tight')



#isolate the sources according to labelled regions

sizes = ndimage.sum(sxz, labels, range(nlabels+1))
sources=np.zeros(sxz.shape)

for i in range(nlabels+1):
            a=(labels==(i+1))
            sources = sources + sxz*a


#isolated sources

fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)
im=ax.imshow(sources,extent=[0,11.25,4.0,0])
cbar=fig.colorbar(im, cax=cax, orientation='vertical')
cbar.ax.tick_params(labelsize=8)
ax.set_xlabel("x (km)")
ax.set_ylabel("z (km)")
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
fig.savefig("isolated_sources.pdf", bbox_inches='tight')



labels1, nlabels1 = ndimage.label(sources)

r, c = np.vstack(ndimage.center_of_mass(sources, labels1, np.arange(nlabels1) + 1)).T
r = r.astype(int);
c = c.astype(int);

zcords= ndimage.measurements.maximum_position(sources, labels=labels1, index=np.arange(1, nlabels1 + 1))

maxcords = np.asarray(zcords)

sy=maxcords[:,0]
sx=maxcords[:,1]

output_file = open('source_coords.bin','wb')
maxcords.tofile(output_file)
output_file.close()
    
#wavelet extraction
wav=np.zeros((nlabels1,nt))

print(maxcords)

string0 = 'invsrc'

for i in range(nlabels1):

    a = np.zeros(shape=(nx,nz))
    a[sx[i],sy[i]] = 1.0;
    bsr=np.tile(a,(nt,1,1));
    s=bsr*src;
    wav[i,:]=np.sum(s,axis=(1,2));
    wlt=wav[i,:]/np.max(abs(wav[i,:]))
    #wlt=gaussian_filter(wlt,sigma=5);
    wlt=wlt.astype('float32')
    plt.figure(i+5);
    plt.plot(wlt);
    
    string=string0+str(i+1)+'.bin'
    output_file = open(string, 'wb')
    wlt.tofile(output_file)
    output_file.close()

### save the wavelet plots as pdf files:

x=np.linspace(0,1.0,nt/6)

f1=open('source1.bin','rb')
f2=open('invsrc1.bin','rb')
isrc=np.fromfile(f2,dtype='float32')
isrc[0:50]=0
src=np.fromfile(f1,dtype='float32')
src=src/np.max(abs(src))

fig, ax =plt.subplots()
plt.plot(x,src[0:(nt/6)],label='Actual',color='k',linewidth=2)
plt.plot(x,isrc[0:(nt/6)],label='Inverted',color='r',linestyle='--',linewidth=2)
plt.legend(loc='upper right')
plt.xlim(0,1.0)
plt.xlabel("t (s)")
plt.ylabel("w(t)")
figstring='wavelet1.pdf'
fig.savefig(figstring, bbox_inches='tight')

f1=open('source2.bin','rb')
f2=open('invsrc4.bin','rb')
isrc=np.fromfile(f2,dtype='float32')
isrc[0:50]=0
src=np.fromfile(f1,dtype='float32')
src=src/np.max(abs(src))

fig, ax =plt.subplots()
plt.plot(x,src[0:(nt/6)],label='Actual',color='k',linewidth=2)
plt.plot(x,isrc[0:(nt/6)],label='Inverted',color='r',linestyle='--',linewidth=2)
plt.legend(loc='upper right')
plt.xlim(0,1.0)
plt.xlabel("t (s)")
plt.ylabel("w(t)")
figstring='wavelet2.pdf'
fig.savefig(figstring, bbox_inches='tight')

f1=open('source3.bin','rb')
f2=open('invsrc3.bin','rb')
isrc=np.fromfile(f2,dtype='float32')
isrc[0:50]=0
src=np.fromfile(f1,dtype='float32')
src=src/np.max(abs(src))

fig, ax =plt.subplots()
plt.plot(x,src[0:(nt/6)],label='Actual',color='k',linewidth=2)
plt.plot(x,isrc[0:(nt/6)],label='Inverted',color='r',linestyle='--',linewidth=2)
plt.legend(loc='upper right')
plt.xlim(0,1.0)
plt.xlabel("t (s)")
plt.ylabel("w(t)")
figstring='wavelet3.pdf'
fig.savefig(figstring, bbox_inches='tight')

f1=open('source4.bin','rb')
f2=open('invsrc2.bin','rb')
isrc=np.fromfile(f2,dtype='float32')
isrc[0:50]=0
src=np.fromfile(f1,dtype='float32')
src=src/np.max(abs(src))

fig, ax =plt.subplots()
plt.plot(x,src[0:(nt/6)],label='Actual',color='k',linewidth=2)
plt.plot(x,isrc[0:(nt/6)],label='Inverted',color='r',linestyle='--',linewidth=2)
plt.legend(loc='upper right')
plt.xlim(0,1.0)
plt.xlabel("t (s)")
plt.ylabel("w(t)")
figstring='wavelet4.pdf'
fig.savefig(figstring, bbox_inches='tight')


