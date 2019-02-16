

import numpy as np

def radius(x,y):
	return np.sqrt(x**2 + y**2)

def theta(x,y):
	return np.arctan2(y,x)

L = 1
Uinf = 1
pinf=0
rho=1

nx=200
ny=200

xmin=-5
xmax=5

ymin=-5
ymax=5

X = np.linspace(-5,5,num=nx)
Y = np.linspace(-5,5,num=ny)


import matplotlib.pyplot as plt
plt.figure()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.axes().set_aspect('equal')

#Plot stream function contours (streamlines)

def psi(r,th,L,Uinf): #Stream function for flow around cylinder
	return (Uinf*r - L/r) * np.sin(th)

PSI = np.zeros((ny,nx))
for j in range(ny):
	for i in range (nx):
		th=theta(X[i],Y[j])
		r=radius(X[i],Y[j])
		PSI[i][j]=psi(r,th,L,Uinf)

plt.contour(X,Y,np.transpose(PSI),(np.linspace(ymin,ymax,20+1)),colors='0.5',linestyles='solid')

#Plot stagnation points

xstag=[np.sqrt(L/Uinf),-np.sqrt(L/Uinf)]
ystag=[0,0]
plt.plot(xstag,ystag,'ko',fillstyle="none")

#Highlight surface stream function (stagnation streamline)

psistag=0
plt.contour(X,Y,np.transpose(PSI),[psistag],colors='black',linestyles='solid',linewidths=[2])

#plt.contour(X,Y,np.transpose(PSI),[1],colors='r',linestyles='solid',linewidths=[1])
#plt.contour(X,Y,np.transpose(PSI),[2],colors='r',linestyles='solid',linewidths=[1])

#Plot velocity potential contours

def phi(r,th,L,Uinf): #Velocity potential for flow around cylinder
	return (Uinf*r + L/r) * np.cos(th)

PHI = np.zeros((ny,nx))
for j in range(ny):
	for i in range (nx):
		th=theta(X[i],Y[j])
		r=radius(X[i],Y[j])
		PHI[i][j]=phi(r,th,L,Uinf)
plt.contour(X,Y,np.transpose(PHI),np.linspace(xmin,xmax,20+1),colors='0.5',linestyles=':')
plt.show()

'''
Velocity field
'''

def ur(r,th,L,Uinf):
	return np.cos(th)*(Uinf - L/r**2)
	
def uth(r,th,L,Uinf):
	return np.sin(th)*(-Uinf - L/r**2)

UR = np.zeros((ny,nx))
UTH= np.zeros((ny,nx))
U= np.zeros((ny,nx))
for j in range(ny):
	for i in range (nx):
		th=theta(X[i],Y[j])
		r=radius(X[i],Y[j])
		UR[i][j] =ur(r,th,L,Uinf)
		UTH[i][j]=uth(r,th,L,Uinf)
		U[i][j] = np.sqrt(UR[i][j]**2 + UTH[i][j]**2)
		
#Plot velocity field
plt.pcolormesh(X,Y,np.transpose(U/Uinf),vmin=0,vmax=2)
cb = plt.colorbar(label=r"$U/U_{\infty}$", orientation='vertical')
#Draw cylinder
R=np.sqrt(L/Uinf) #Cylinder radius
cylinder = plt.Circle((0,0), radius=R, color='black', alpha=0.5)
plt.gcf().gca().add_artist(cylinder)
plt.show()
		
'''
Pressure field
'''

P = np.zeros((ny,nx))
for j in range(ny):
	for i in range (nx):
		P[i][j]=pinf+0.5*rho*(Uinf**2 - U[i][j]**2)
cp=(P-pinf)/(0.5*rho*Uinf**2) #Pressure coefficient

#Plot pressure coefficient
plt.figure()
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.axes().set_aspect('equal')
plt.contour(X,Y,np.transpose(PSI),(np.linspace(ymin,ymax,20+1)),colors='0.5',linestyles='solid')
plt.plot(xstag,ystag,'ko',fillstyle="none")
plt.contour(X,Y,np.transpose(PSI),[psistag],colors='black',linestyles='solid',linewidths=[2])
plt.contour(X,Y,np.transpose(PHI),np.linspace(xmin,xmax,20+1),colors='0.5',linestyles=':')

plt.pcolormesh(X,Y,np.transpose(cp),vmin=-1,vmax=1)
cb = plt.colorbar(label=r"$c_p = \frac{p-p_{\infty}}{\frac{1}{2} \rho U_{\infty}^2}$", orientation='vertical')
#Draw cylinder
R=np.sqrt(L/Uinf) #Cylinder radius
cylinder = plt.Circle((0,0), radius=R, color='black', alpha=0.5)
plt.gcf().gca().add_artist(cylinder)
plt.show()

'''
Surface pressure
'''

plt.figure()
TH=np.linspace(0,-np.pi,100)
#ptop=pinf+0.5*rho*Uinf**2 * (1-4*np.sin(TH)**2) 
#ptop=(ptop-pinf)/(0.5*rho*Uinf**2)
ptop=1-4*np.sin(TH)**2
#pbot=pinf+0.5*rho*Uinf**2 * (1-4*np.sin(-TH)**2) 
#pbot=(pbot-pinf)/(0.5*rho*Uinf**2)
pbot=1-4*np.sin(-TH)**2
plt.plot(TH*180/np.pi,ptop,"k:")
plt.plot(TH*180/np.pi,pbot,"r--")
plt.xlabel(r"$\theta$ (degrees)")
plt.ylabel(r"$c_p$")
