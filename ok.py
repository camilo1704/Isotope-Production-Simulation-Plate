#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.integrate as integrate
from Energia_x import *
#-------------------------
with open('prueba.txt') as f:
    A=np.loadtxt((x.replace(b':',b' ') for x in f))
print A
#interp=interpolate.interp1d(tiempo,volts,'slinear')
#	ENER=np.array(E_inicial)
