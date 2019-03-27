#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.integrate as integrate
#-----------------------------------
#devuelve SE de E en MeV
def SE(E):
	f=np.loadtxt('/home/cuchuflito/Documents/PI/123Te(p,n)123I/123Te(p,n)123I_takacs.txt')
	SE=f[:,1]/1000	#en barns
	E1=f[:,0]
	interp=interpolate.interp1d(E1,SE,'slinear')
	if E>=E1[0]:
		return interp(E)
	else:
		return 0

	
	
