#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.integrate as integrate
from Energia_x import *
from SE import *
import matplotlib.mlab as mlab
import math
from mpl_toolkits.mplot3d import Axes3D
#-------------------------
RHO=5.67 	#g/cm3
m=122.9 	#g/cm3
N=RHO*6.02E-1/m	#E23 de avogadro por E-24 de barns, densidad atomica
carga=6.2415E12 	# transformacion de cargas/seg a micro ampere



def nucleos(Ener,dx):
	#tener en cuenta q es por cada micro ampere de particulas incidentes
	Energias=E(Ener,dx)  #me devuelve la energia de cada posicion con un espaciamiento dx
	dy=np.array([])
	for i in np.arange(0,np.where(Energias>1E-8)[0].shape[0]-1,1):
		dy=np.append(dy,SE(Energias[i+1])*np.exp(-dx*N*SE(Energias[i+1])*i)*N*carga*dx)
	#print np.sum(dy)/3.7E7 #Actividad de saturacion
	return dy
def QdE(Ener,dx):
	#n es nucleos en xi
	#dx es el paso
	#calculamos Q debido a dE/dx
	# tener en cuenta que es por cada particula incidente
	QdE=np.array([])
	Energias=E(Ener,dx)
	for i in np.arange(0,Energias.shape[0]-1,1):
		QdE=np.append(QdE,Energias[i]-Energias[i+1])
	
	return QdE
	
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
	while r>1E-9:
		r=ENER[ENER.size-1]-interp(ENER[ENER.size-1])*dx*RHO
		if r>0:
			ENER=np.append(ENER,r)
		i=i+1
	print ENER
		
	
	
	return ENER
def QNucl(E,dx):
	#devuelve la energía liberada ene MeV en cada xi debido a reaccion nuclear 
	#por cada micro ampere de corriente incidente
	print nucleos
	return nucleos(E,dx)*2.0108
def Q(E,dx,xmax):
	#Devuelve Qtot por cada micro ampere irradiado en Joule
	Q1=QNucl(E,dx)
	Q2=QdE(E,dx)*carga
	index=np.where(np.arange(0,Q2.shape[0]*dx,dx)<xmax)[0]
	Qnuevo=np.zeros(index.shape[0])
	for i in index:
		if (Q1[i]>1E-8) & (Q2[i]>1E-8):
			Qnuevo[i]=Q1[i]+Q2[i]
		else:
			if Q1[i]>1E-8:
				Qnuevo[i]=Q1[i]
			else: 
				if Q2[i]>1E-8:
					Qnuevo[i]=Q2[i]
							
	#print Qnuevo*1.60218e-13*20
	return Qnuevo*1.60218e-13
def Normaldist(x,dx,I_0):
	mu = 0
	
	sigma_x = math.sqrt(x/3.)
	sigma_y= math.sqrt(x/3.)
	x = np.arange(mu - 3*sigma_x, mu + 3*sigma_x, dx)
	y = np.arange(mu - 3*sigma_x, mu + 3*sigma_x,dx)
	M=np.zeros(y.shape[0]*x.shape[0]).reshape(x.shape[0],y.shape[0])
	
	for i in np.arange(0,x.shape[0],1):
		for j in np.arange(0,y.shape[0],1):
			M[i,j]=mlab.normpdf(x[i], mu, sigma_x)*mlab.normpdf(y[j], mu, sigma_y)
	#mlab.normpdf(x, mu, sigma_x)
	fig = plt.figure()
	ax=fig.gca(projection='3d')
	X,Y=np.meshgrid(x,y)
	
	Normalizacion=I_0/np.sum(M)
	
	M=M*Normalizacion	#normaliza para que integral de la gaussiana sea I_0
	ax.plot_surface(X,Y,M)
	plt.plot(x,mlab.normpdf(x, mu, sigma_x))
	plt.show()
	print M[int(M[0].shape[0]/2),int(M[0].shape[0]/2)]
	#return M
def Q_matrix(x,dx,I_0,E_max,x_max):
	#x_max es donde profundidad donde termina el target 
	M=Normaldist(x,dx,I_0)
	N=np.array([0])
	for i in np.arange(0,M[:,0].shape[0],1):
		for j in np.arange(0,M[:,0].shape[0],1):
			print Q(E_max,dx,x_max).shape[0]
			N=np.vstack((N,Q(E_max,dx,x_max)))
	#M[:,0].shape[0]/2
	#plt.plot(np.arange(0,M[M[:,0].shape[0]/2,M[:,0].shape[0]/2].shape[0],dx),M[:,0].shape[0]/2,M[:,0].shape[0]/2]),'-r')
	#plt.show()	
#return M
def Temp(E,dx):
	def k(T):
		return 100*(0.009875+4.699/T+215.4/(T*T))
	Qgen=0.000108*Q(E,dx,1)
	T=np.ones(Qgen.shape[0])*100
	dx=dx/100
	k=2
	h_a=2436.53
	h_he=14.34
	T_inf_he=273
	T_inf_a=273
	T[0]=(k*200+dx*h_he*T_inf_he+Qgen[0])/(dx*h_he+k)	
	print Qgen.shape[0]
	for i in np.arange(0,100000):
		for m in np.arange(1,Qgen.shape[0]-1,1):
			T[m]=(T[m-1]+T[m+1]+Qgen[m]/k)/2
		T[Qgen.shape[0]-1]=(k*T[Qgen.shape[0]-2]+h_a*dx*T_inf_a+Qgen[Qgen.shape[0]-1])/(k+h_a*dx)
		T[0]=(k*T[1]+dx*h_he*T_inf_he+Qgen[0])/(dx*h_he+k)
#	for m in np.arange(0,L,1)
	#	for n in np.arange(0,L,1)
	#		for z in np.arange(0,L,1)
	#			if (n==0 & m==0):
	#				if z==0:
	#					T(m,n,z)=(k(T(m,n,z))+Qgen(m,n,z)*dx*dx*dx*dx+T_inf_he*h_he*dx)/(dx*h_he+k)
	#				else:
	#T(0,0,0)=(T(0,1,0)+T(0,0,1)+T(1,0,0)+Q(0,0,0)/(dx*k)+h_he*dx*T_inf_he/k)/(3+h_he*dx/k)
	#for i in np.arange(1,L,1):
	#	T(0,0,i)=(T(0,0,i-1)+T(1,0,i)+T(0,1,z)+T(0,0,z+1)+Q(0,0,Z)/(k*dx)+h_he*dx*T_inf_he/k)/(4+h_he*dx/k)
	#T(0,0,L)=(T(1,0,L)+T(0,1,L)+T(0,0,L-1)+Q(0,0,L)/(dx*k)+h_he*dx*T_inf_he/k)/(3+h_he*dx/k)
	
#	T(0,i,0)=(T(1,i,0)+T(0,i,1)+T(0,L-1,0)+T(0,L+1,0)+Q(0,i,0)/(dx*k)+h_he*dx*T_inf_he/k)/(4+h_he*dx/k)
#	for i in np.arange(1,L,1):
#		T(0,L,i)=(T(0,L-1,Z)+T(1,L,i)+T(0,L,i+1)+T(0,L,i-1)+Q(0,L,i)/(dx*k)+h_he*dx*T_inf_he/k)/(4+h_he*dx/k)
#	T(0,L,L)=(T(1,L,L)+T(0,L-1,L)+T(0,L,L-1)+Q(0,L,L)/(dx*k)+h_he*dx*T_inf_he/k)/(3+h_he*dx/k)
#	for i in np.arange(1,L,1):
#		T(0,i,L)=(T(0,I-1,L)+T(0,i,L-1)+T(1,i,L)+T(0,i+1,L)+Q(0,i,L)/(dx*k)+h_he*dx*T_inf_he/k)/(4+h_he*dx/k)
#	for i in np.arange(1,L,1):
#		for j in np.arange(1,L,1):
#			T(0,i,j)=(T(0,i+1,j)+T(0,i-1,j)+T(0,i,j-1)+T(0,i,j+1)+T(1,i,j)+Q(0,i,j)/(dx*k)+h_he*T_inf*dx/k)/(5+h_he*dx/k)
#	for i in np.arange(1,L,1):
#		T(i,0,0)=(T(i+1,0,0)+T(i-1,0,0)+T(i,1,0)+T(i,0,1)+Q(i,0,0)/(dx*k))/4.
#	T(L,0,0)=(T(L-1,0,0)+T(L,1,0)+T(L,0,1)+Q(L,0,0)/(dx*k)+h_a*dx*T_inf_a/dx)/(4+h_a*dx/k)
#	for i in np.arange(1,L,1):
#		T(L,i,0)=(T(L-1,i,0)+T(L,i-1,0)+T(L,i+1,0)+T(L,i,1)+Q(L,i,0)/(dx*k)+h_a*dx*T_inf_a/k)/(4+h_a*dx/k)
#	for i in np.arange(1,L,1):
#		for j in np.arange(1,L,1):
#			T(L,L,0)=(T(L-1,L,0)+T(L,L-1,0)+T(L,L,1)+Q(L,L,0)/(dx*k))/(3+h_a*dx/k)

		
#	plt.plot(np.arange(0,Qgen.shape[0],1)*dx,T-273,'or')
#	plt.show()
Normaldist(3,0.05,20)
