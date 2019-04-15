import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert


'''
Types of source wavelet

0 : Ricker wavelet
1 : 1st derivative of Gaussian distribution
2 : Fucher-Mullers wavelet
3 : sine cube 
4 : Bandlimited spike 

phase rotation done  using Hilbert 

change type parameter in dictionary accordingly 

'''

def wavelet(amp,nt,dt,src):
	
	t=np.linspace(0,dt*nt,nt);
	fc=src['fc'];
	ts=src['tshift'];
	t0=1.0/fc;
	phase=src['phase'];
	
	if(src['type'] == 0):
		
		#Ricker Wavelet
		tmp = (np.pi*fc*(t-ts-t0))**2;
		wt = amp*(1.0-2.0*tmp)*np.exp(-tmp);
		if(phase):
			phase=phase*np.pi/180; #radians
			hlb=hilbert(wt);
			wt=np.cos(phase)*np.real(hlb)-np.sin(phase)*np.imag(hlb);
		return wt;
	
	elif(src['type'] == 1):
		
		# First derivative Gaussian
		ts=1.2*ts;
        tmp=(np.pi*fc)**2;
        wt = -2.0*amp*tmp*(t-ts-t0)*np.exp(-tmp*(t-ts)**2);
        if(phase):
			phase=phase*np.pi/180; #radians
			hlb=hilbert(wt);
			wt=np.cos(phase)*np.real(hlb)-np.sin(phase)*np.imag(hlb);
			return wt;
        
	elif(src['type'] == 2):	
		
		# Fucher-Mullers wavelet 
		wt=amp*(np.sin(2*np.pi*(t-ts)*fc)-0.5*np.sin(4.0*np.pi*(t-ts)*fc));
		wt[:int(ts/dt)]=0.0;
		wt[int((t0+ts)/dt):]=0.0;
		if(phase):
			phase=phase*np.pi/180; #radians
			hlb=hilbert(wt);
			wt=np.cos(phase)*np.real(hlb)-np.sin(phase)*np.imag(hlb);
		return wt;
		
	elif(src['type'] == 3):
		
		# Sine raised to power 3
		wt=amp*np.sin(np.pi*(t-ts)/t0)**3;
		wt[:int(ts/dt)]=0.0;
		wt[int((t0+ts)/dt):]=0.0;
		if(phase):
			phase=phase*np.pi/180; #radians
			hlb=hilbert(wt);
			wt=np.cos(phase)*np.real(hlb)-np.sin(phase)*np.imag(hlb)
		return wt
		

'''
#source dictionary format

src={ 'fc' : 10 , #peak frequency
	  'tshift' : 0.5, #time shift seconds
	  'type' : 0, #source type
	  'phase' : 0, #phase shift in degrees
	}

'''



#modeling parameters
nz=160; nx=450; nt=3000; dt=0.002;
amp=1;



#Main

#source 1
src1= {'fc': 10,'tshift': 0.2,'type': 0, 'phase':0 }
wt1=wavelet(amp,nt,dt,src1);
wt1.astype('float32').tofile('source1.bin');

#source 2
src2= {'fc': 10,'tshift': 0.4,'type': 2, 'phase':0 }
wt2=wavelet(amp,nt,dt,src2);
wt2.astype('float32').tofile('source2.bin');

#source 3
src3= {'fc': 10,'tshift': 0.2,'type': 3, 'phase':0 }
wt3=wavelet(amp,nt,dt,src3);
wt3.astype('float32').tofile('source3.bin');

#source 4
src4= {'fc': 10,'tshift': 0.3,'type': 0, 'phase':180}
wt4=wavelet(amp,nt,dt,src4);
wt4.astype('float32').tofile('source4.bin');


#Plot
plt.figure(0);
plt.plot(wt1,color='k'); 
plt.figure(1);
plt.plot(wt2,color='r');
plt.figure(3);
plt.plot(wt3,color='b');
plt.figure(4);
plt.plot(wt4,color='g');

plt.show()
