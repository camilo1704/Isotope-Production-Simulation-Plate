#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.integrate as integrate
#-------------------------

RHO=5.67	#g/cm3



def E(E_inicial,dx):			
	#PASAR dx en cm
	#E sale en MeV
	f=np.loadtxt('SP.txt')
	E=f[:,0]
	SP=f[:,1]
	interp=interpolate.interp1d(E,SP,'slinear')
	ENER=np.array(E_inicial)
	i=0
	r=ENER-interp(ENER)*dx*RHO
	ENER=np.append(ENER,r)
	#print ENER.size
	while r>1E-9:
		r=ENER[ENER.size-1]-interp(ENER[ENER.size-1])*dx*RHO
		if r>0:
			ENER=np.append(ENER,r)
		i=i+1
		
	rang=i*dx
	return ENER
