# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 18:16:22 2019

@author: mbast
"""

import numpy as np
import matplotlib.pyplot as plt 


from functionCall import forcing 
from functionCall import exact
from functionCall import exactTimeHistory


tFinal = 4
numberOfTimeSteps = 20
lanbda = -3

deltaT = tFinal/numberOfTimeSteps

#1D arrays of zeros 
yTrapezoidal = np.zeros((numberOfTimeSteps+1))
yRungeKutta = np.zeros((numberOfTimeSteps+1))
times = np.zeros(numberOfTimeSteps+1)

#initial conditions

yTrapezoidal[0] = exact(0)
yRungeKutta[0] = exact(0)


 i = 0

for i in range(numberOfTimeSteps):
    tn = times[i]
    tnpo = tn + deltaT
    #time step for Runge-Kutta
    tnpr = tn+deltaT*0.5
    times[i+1] = tnpo
    
    yTrapezoidal[i+1] = (yTrapezoidal[i]+(deltaT/2)*(lanbda*yTrapezoidal[i]+forcing(tn)+forcing(tnpo)))/(1-deltaT*0.5*lanbda) 
    
    yHat1 = yRungeKutta[i]+ (deltaT/2)*(lanbda*yRungeKutta[i]+forcing(tn))
    yHat2 = yRungeKutta[i]+ (deltaT/2)*(lanbda*yHat1+forcing(tnpr))
    yHat3 = yRungeKutta[i]+ (deltaT)*(lanbda*yHat2+forcing(tnpr))
    
    yRungeKutta[i+1] = yRungeKutta[i]+(deltaT/6)*((lanbda*yRungeKutta[i]+forcing(tn))+2*(lanbda*yHat1+forcing(tnpr))+
               2*(lanbda*yHat2+forcing(tnpr))+(lanbda*yHat3+forcing(tnpo)))
    
    
print("Trapezoidal error =", abs(exact(tFinal)-yTrapezoidal[i+1]))
print("Runge Kutta error =", abs(exact(tFinal)-yRungeKutta[i+1]))


x,y = exactTimeHistory(tFinal,5000)
    

#plotting Trapezoidal 
plt.figure()
plt.title(r"$Trapezoidal$ $Method$")
plt.plot(times,yTrapezoidal,'o-', label = "Trapezoidal")
plt.grid()
plt.plot(x,y, label = 'exact solution')
plt.xlabel(r"$time$")
plt.ylabel(r"$Y$")
plt.legend()

#ploting Runge-Kutta
plt.figure()
plt.grid(b = None, which = 'both')
plt.title(r"$Runge-Kutta$ 4 $Method$")
plt.plot(times,yRungeKutta,'o-', label = 'Runge-Kutta')
plt.plot(x,y,label = "exact solution")
plt.xlabel(r"$time$")
plt.ylabel(r"$Y$")
plt.legend()

#ploting exact solution
plt.figure()
plt.grid()
plt.title(r"$Exact$ $solution$")
plt.plot(x,y, color = 'orange')


