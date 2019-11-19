# -*- coding: utf-8 -*-

# -----------------------------------------------------
# MR Signal Python Library
# -----------------------------------------------------
#
# Basic MRI signal simulations, intended primaraly for
# learning and low-to-medium complexity MRI simulations.
#
# Derived from Stanford RAD229 Class (Matlab) functions
#
# Created on Tue Nov 19 08:57:48 2019
# author: Joshua Kaggie, Brian Hargreaves
# -----------------------------------------------------


import numpy as np

def relax(t,T1 = 4., T2 = 0.1, combine=True):
    T1 = T1 * 1.
    T2 = T2 * 1.
    t = t * 1.
    E1 = np.exp(-t/T1)
    E2 = np.exp(-t/T2)
    A = np.diag([E2, E2, E1])
    B = np.array([0,0,1-E1])        
    B = np.reshape(B,[3,1])
    if combine:
        A = np.hstack([A,B])   ### ADD TO FOURTH AXIS!  HSTACK????
        return A
    return A,B


#--------------------------------------------------------
#  By convention all rotations are left-handed
#--------------------------------------------------------

# Returns 3x3 matrix for left-handed rotation about x
def xrot(angle = 0., in_degs = True):
    if in_degs:
        angle = angle*np.pi/180.
    c = np.cos(angle)    
    s = np.sin(angle)    
    M = np.array([[1.,0.,0.],[0., c, s],[0,-s, c]])
    return M

# Returns 3x3 matrix for left-handed rotation about y
def yrot(angle = 0., in_degs = True):
    if in_degs:
        angle = angle*np.pi/180.
    c = np.cos(angle)    
    s = np.sin(angle)    
    M = np.array([[c,0.,-s],[0.,1.,0.],[s,0.,c]])
    return M


# Returns 3x3 matrix for left-handed rotation about z
def zrot(angle = 0., in_degs = True):
    'This is equivalent to help'
    if in_degs:
        angle = angle*np.pi/180.
    c = np.cos(angle)    
    s = np.sin(angle)    
    M = np.array([[c,s,0.],[-s,c,0],[0.,0.,1.]])
    return M


# Returns 3x3 matrix for rotation about axis in x-y plan phi away from x
def throt(theta = 0., phi=0., in_degs = True):
    'This is equivalent to help'
    if in_degs:
        theta = theta*np.pi/180.
        phi = phi*np.pi/180.
    
    ca = np.cos(theta)
    sa = np.sin(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)
    M = np.array([[cp*cp+sp*sp*ca, cp*sp*(1-ca), -sp*sa],
     [cp*sp-sp*cp*ca, sp*sp+cp*cp*ca, cp*sa],
     [sa*sp, -sa*cp, ca]])
    
    return M        

# Converts a time vector to the corresponding frequency vector of an FFT
def time2freq(t):
    dt = t[1]-t[0]
    num_p = len(t)
    df = 1./ num_p/dt
    f = np.arange(np.ceil((num_p-1.)/2.), np.floor((num_p-1.)/2.))*df
    return f


# abprop propagates a series operations Ai,Bi into a single A,B
# A = ... A3*A2*A1    B = ... B3 + A3*(B2 + A2*B1
#
# Example:
# TR = 1
# TI = 0.5
# TE = 0.05
# T1 = 0.5
# T2 = 0.1    
# abprop(relax(TR-TE-TI,T1,T2, combine = True),xrot(180.))
# Should return mss = [0,0,-0.42]    
def abprop(*varargin):        
    A = np.eye(3)
    B = np.array([[0],[0],[0]])
    mss = []

    argcount = 0
    nargin = len(varargin)
    
    while (argcount < nargin):
        Ai = varargin[argcount]
        argcount += 1
        
        if np.shape(Ai)[0] != 3:
            print('The number of rows is not 3.')
        if np.shape(Ai)[1] == 3:
            if argcount < nargin:
                Bi = varargin[argcount+1]
                if np.shape(Bi) != [1,3]:
                    Bi = np.array([[0],[0],[0]])
            if argcount == nargin:
                Bi = np.array([[0],[0],[0]])
        elif np.shape(Ai)[1] == 4:              
            Bi = Ai[:,3:]
            Ai = Ai[:,:3]               
        else:
            print('Arg count '+str(np.shape(Ai))+' should be 3x3 or 3x4')
            
        A = np.matmul(Ai,A)        
        B = np.matmul(Ai,B)+Bi

    mss = np.matmul(np.linalg.inv(np.eye(3)-A),B)
    return A,B, mss
     
     

import matplotlib.pyplot as plt

def magphase(x,arr):
    mag = np.abs(arr)
    phase = np.angle(arr)
    
    fig1 = plt.subplot(2,1,1)
    plt.plot(x,phase/np.pi)
    
    
    




     
     

