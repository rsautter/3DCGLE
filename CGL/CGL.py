import numpy as np
from numpy.fft import fftn,ifftn,fftfreq
import math
import random
import itertools
import tqdm as tqdm
from collections import namedtuple
from scipy.signal import convolve2d

class CGL():
	'''
	CGL - Complex Ginzburg-Landau
	
	Wrote by: Rubens Andreas Sautter (2022)
	
	Adapted from Aranson, et.al.(1997)
	https://arxiv.org/abs/patt-sol/9709005
	
	Adittional References:

	de Franciscis, d’Onofrio (2012) (Tsallis-Borland model)
	https://journals.aps.org/pre/abstract/10.1103/PhysRevE.86.021118
	
	
	
	Complex Ginzburg-Landau equation solved with Fourier pseudospectral methods, and integrated with RK45.

	'''

	def __init__(self, c1=1.0, c2=1.0,h=1.0, msize = 128, ic='r',a0=1.0,b0=0.0,dim=2):
		'''
		Spatial parameters:
			ic = initial condition('r', 'g')
			h - grid spacing
			dim - dimension of the data (integer)

		GL parameters:
			c1 - diffusion parameter - (1+ib)Nabla A
			c2 - reaction parameter  - (1+ic)(|A|^2)A
			
		Noise Parameters:
			noiseSpeed - ]0,1[ - the speed (relative to the number of iterations) which the noise moves
			sigma_r - reactive noise 'strenght'
			noiseArgs - Colored noise parameters {'beta':2,std = 0.01}
		'''
		
		self.c1, self.c2 = c1,c2
		self.a0 = a0
		self.b0 = b0
		
		self.h = h
		self.ic = ic
		self.msize = msize
		self.dim = dim

	def __getRandom(self,n,dim):
		newShape = tuple([n for i in range(dim)])
		return np.random.rand(n**dim).reshape(newShape)
		
	def __getGaussian(self,n,dim):
		out = np.zeros(np.repeat(n,dim))
		c = n/2
		squareDists = np.sum((np.indices(out.shape)-c)**2,axis=0)
		return 2*(np.exp(-squareDists/(20*n))-0.5)
		
	def getInitialCondition(self):
		if self.ic=='r':
			self.a = 2*self.a0*((self.__getRandom(self.msize,self.dim)-0.5)+1j*(self.__getRandom(self.msize,self.dim)-0.5))+self.b0
		else:
			self.a = self.a0*(self.__getGaussian(self.msize,self.dim)+1j*self.__getGaussian(self.msize,self.dim))
			
		return np.array(self.a)
		
	
	def getChainedSingleReaction(self,a0=None,dt=0.1, nit=3000):
		'''
		Returns the iteration of a single amplitude (spatial part is ignored)
		
		The function integrates with rk4 method
		'''
		states = []
		delta = 1e-6*(np.random.rand()-0.5)
		if a0 is None:
			at = self.a0+delta
		else:
			at = a0
				
		for i in range(nit):
			states.append(at)
			
			t = i*dt
			k1 = self.reaction(at,t)
			k2 = self.reaction((at+dt*k1/2), t+dt/2.)
			k3 = self.reaction((at+dt*k2/2), t+dt/2.)
			k4 = self.reaction((at+dt*k3), t+dt)
			at = at + dt*(k1+2*k2+2*k3+k4)/6.
		return np.array(states)
	
		
	def reaction(self, a, t):
		a1 = a - (1+1j*self.c2)*(np.abs(a)**2)*a
		return np.array(a1)
		
	def solveRKF45(self,dt,ntimes,stepsave,dtTolerace=1e-5):
		state = self.getInitialCondition()
		times = []
		states = [state]	
			
		w = np.array([	[					0,0,0,0,0,0],
				[1/4,					0,0,0,0,0],
				[3/32,9/32,				0,0,0,0],
				[1932/2197,-7200/2197,7296/2197,	0,0,0],
				[439/216,-8,3680/513,-845/4104,	0,0],
				[-8/27, 2,-3544/2565,1859/4104,-11/40,0]
			])
		t = 0.0
				
		for time in tqdm.tqdm(range(ntimes)):
			while t < time*dt:
				step = dt
				
				k1 = step*self.timeDerivatives(state,								t		)
				k2 = step*self.timeDerivatives(state+k1*w[1,0], 		        				t+step/4	)
				k3 = step*self.timeDerivatives(state+k1*w[2,0]+k2*w[2,1], 					t+3*step/8	)
				k4 = step*self.timeDerivatives(state+k1*w[3,0]+k2*w[3,1]+k3*w[3,2],    				t+12*step/13	)
				k5 = step*self.timeDerivatives(state+k1*w[4,0]+k2*w[4,1]+k3*w[4,2]+k4*w[4,3],    		t+step		)
				k6 = step*self.timeDerivatives(state+k1*w[5,0]+k2*w[5,1]+k3*w[5,2]+k4*w[5,3]+k5*w[5,4],    	t+step/2	)
				
				approach4 = state + (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 -k5/5
				approach5 = state + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6
				
				error = np.max(np.abs(approach4-approach5))
				if error> dtTolerace:
					step = dt*((dtTolerace/(2*error))**.25)
				
					k1 = step*self.timeDerivatives(state,								t		)
					k2 = step*self.timeDerivatives(state+k1*w[1,0], 		        				t+step/4	)
					k3 = step*self.timeDerivatives(state+k1*w[2,0]+k2*w[2,1], 					t+3*step/8	)
					k4 = step*self.timeDerivatives(state+k1*w[3,0]+k2*w[3,1]+k3*w[3,2],    				t+12*step/13	)
					k5 = step*self.timeDerivatives(state+k1*w[4,0]+k2*w[4,1]+k3*w[4,2]+k4*w[4,3],    		t+step		)
					k6 = step*self.timeDerivatives(state+k1*w[5,0]+k2*w[5,1]+k3*w[5,2]+k4*w[5,3]+k5*w[5,4],    	t+step/2	)
					
					approach4 = state + (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 -k5/5
					
				t += step
				state = approach4 
			
			times.append(t)
			if time in stepsave:
				states.append(state)
		return np.array(states), np.array(times)
		
	def timeDerivatives(self,state,time):
		rFtState = fftn(np.real(state))
		iFtState = fftn(np.imag(state))
		
		# spectral shift variables
		specR = np.zeros(state.shape)
		specI = np.zeros(state.shape)
		
		# for every frequency dimension, updates the shift variable
		for i in range(len(state.shape)):
			
			# changing the dimension of the measured frequency since numpy's fftfreq is always 1D
			seq = np.ones(len(state.shape)).astype(int)
			fx  = 2*np.pi*fftfreq(state.shape[i])
			seq[i] = len(fx)
			fx = fx.reshape(*tuple(seq))
			
			#makes the shift
			specR = specR - (fx**2)*rFtState
			specI = specI - (fx**2)*iFtState
			
		# measures the laplacian from pseudospectral method
		lap  =    np.real(ifftn(specR)) + 1j*np.real(ifftn(specI))/(self.h**2)
		return (1+self.c1*1j)*np.array(lap) + self.reaction(state,time)
