#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.integrate as integrate
from Energia_x import *
from SE import *
#-------------------------

def Q(Ener,dx):
	#n es nucleos en xi
	#dx es el paso
	#calculamos Q debido a dE/dx
	QdE=np.array([0])
	Energias=E(Ener,dx)
	for i in np.arange(0,Energias.shape[0]-1,1):
			QdE=np.append(QdE,Energias[i]-Energias[i+1])
		
	return QdE



