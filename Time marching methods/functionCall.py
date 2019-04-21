# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 19:01:04 2019

@author: mbast
"""

import numpy as np


######################################
#          function definition       #
######################################
def exact(t):
    
    return np.cos(3*t)*np.tanh(t/3)

def forcing(t):
    
    return (1-(np.tanh(t/3)*np.tanh(t/3)))*(np.cos(3*t)/3)+3*np.tanh(t/3)*(np.cos(3*t)-np.sin(3*t))

def exactTimeHistory(tFinal,n):
    
    dt = tFinal/(n-1)
    exactTimes = np.zeros(n)
    exactSolutions = np.zeros(n)
    
    i=0
    for i in range(n):
        exactTimes[i] = dt*(i-1)
        exactSolutions[i] = exact(exactTimes[i])
    
    return exactTimes, exactSolutions
    



#########################################
    

