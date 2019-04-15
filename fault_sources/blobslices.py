import numpy as np
from scipy import ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter

nz=160; nx=450;nt=3000;dt=0.002;ds=0.025


#Reading Binary File

f=open('inverted_src3d.bin','rb')
src=np.fromfile(f,dtype='float32')
src=src.reshape((nt,nx,nz),order='f')

#Spatial source s(x,z)
sxz=np.sum(np.square(src),axis=0).T;
sxz=np.sqrt(sxz);

fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)
im=ax.imshow(sxz,extent=[0,11.25,4.0,0], cmap='jet')
cbar=fig.colorbar(im, cax=cax, orientation='vertical')
cbar.ax.tick_params(labelsize=8)
ax.scatter(x=np.linspace(5.4, 5.775, 16), y=np.linspace(1.5, 1.125, 16), c='w',s=40);
ax.set_xlim(0, 11.25)
ax.set_ylim(4.0, 0.0)
ax.set_xlabel("x (km)")
ax.set_ylabel("z (km)")
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
fig.savefig("power_blobs.pdf", bbox_inches='tight',pad_inches=0)
