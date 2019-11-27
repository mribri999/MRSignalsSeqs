# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:39:31 2019

@authors: Joshua Kaggie, Brian Hargreaves
"""
import numpy as np
import matplotlib.pyplot as plt
import mrsigpy as mrs

def exampleA1_63():
    # Plot of Mz for IR Question.    
    t = np.arange(0, 1, .01)
    Mz = -1.63*np.exp(-t/0.5)+1;
    Mz[51:] = 1-np.exp(-(t[51:]-.5)/.5);
    plt.plot(t,Mz)
    mrs.lplot('Time (s)','M_z / M_0','Simple IR Signal Example',[0, 1, -1, 1]);
    mrs.setprops()




def examplenoise():
    #%
    #%	This example is to show the effects of noise in MRI.
    #%
    sigma = 0.1;		
    N = 256;
    
    #% A:  Gaussian noise in k-space
    #%
    knoise1 = randn(N,N) + i*randn(N,N);	% Complex noise, sigma=1
    kn = knoise1*sigma;
    
    #% Show noise distribution
    figure(1); 
    ghist(real(kn(:)),[],[],[],'k-space real(noise)','k-space noise');
    
    #% B:  Real-Valued Noise in image space
    #%
    imn = nft(kn);			% FFT, normalized by sqrt(N) to preserve energy
    isig = std(real(imn(:)));	% Get standard deviation to estimate sigma
    
    figure(2);
    ghist(real(imn(:)),[],[],[],'Image real(noise)','image real noise');
    
    
    #% C:  Magnitude Noise in image space
    imag = abs(imn);		% Magnitude of image.
    imsig = sqrt(2-pi/2)*isig;	% Estimate Rayleigh sqrt(variance)
    immean = sqrt(pi/2)*isig;	% Estimate Rayleigh mean
    
    figure(3); 
    ghist(imag(:),[],[],[0:0.01:1]*6*imsig,'Image Magnitude Noise','image abs(noise)');
    
    #% D:  Magnitude Noise+Signal in image space
    figure(4); hold off;
    imsignal = 5*sigma;		% SNR = 5;
    imag = abs(imn+imsignal);	% Magnitude of image.
    imsig = std(imag(:));		% Estimate of Rician sqrt(variance)
    immean = mean(imag(:));		% Estimate of Rician mean, to fit.
    
    figure(4); 
    ghist(imag(:),[],[],[],'Image abs(Signal+Noise)','abs(signal+noise)');
        
    im = zeros(N,N);		% Allocate Image
    im(N/4+1:3*N/4,N/4+1:3*N/4)=1;	% Central square with amplitude 1.
