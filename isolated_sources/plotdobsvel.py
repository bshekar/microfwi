import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import rcParams

nz=160 ;nx=450;nt=3000;dt=0.002;ds=0.025;ng=400;

#font = {'size'   : 18}
#plt.rc('font', **font)

SMALL_SIZE = 8
MEDIUM_SIZE = 14
BIGGER_SIZE = 18

plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize

rcParams['figure.figsize'] = 12, 9

#smooth velocity model
f1=open('sm_win_overt.bin','rb')
init=np.fromfile(f1,dtype='float32')
init=init.reshape(nz,nx,order='f');

fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)
im=ax.imshow(init,extent=[0,11.25,4.0,0],cmap='jet')
cbar=fig.colorbar(im, cax=cax, orientation='vertical')

ax.set_xlabel("x (km)") 
ax.set_ylabel("z (km)")
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
cbar.ax.set_title("V (km/s)")
cbar.ax.tick_params(labelsize=BIGGER_SIZE)
fig.savefig("smooth_vel.pdf", bbox_inches='tight')




#original velocity model
f1=open('win_overt.bin','rb')
vel=np.fromfile(f1,dtype='float32')
vel=vel.reshape(nz,nx,order='f');

fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)
im=ax.imshow(vel,extent=[0,11.25,4.0,0],cmap='jet')
cbar=fig.colorbar(im, cax=cax, orientation='vertical')
#ax.set_title("Actual Velocity Model")
ax.set_xlabel("x (km)") 
ax.set_ylabel("z (km)")
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
cbar.ax.set_title("V (km/s)")
cbar.ax.tick_params(labelsize=BIGGER_SIZE)
fig.savefig("orig_vel.pdf", bbox_inches='tight')


#observed data
f1=open('dobs.bin','rb')
dobs=np.fromfile(f1,dtype='float32')
dobs=dobs.reshape(ng,nt,order='f');
dobs=dobs.T;
fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)
im=ax.imshow(dobs,extent=[0,11.25,6.0,0],cmap='gray',vmin=-1.0,vmax=1.0)
ax.set_xlabel("x (km)")
ax.set_ylabel("t (s)")
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
cbar=fig.colorbar(im, cax=cax, orientation='vertical')
fig.savefig("figures/dobs.pdf", bbox_inches='tight')


#plt.show()
