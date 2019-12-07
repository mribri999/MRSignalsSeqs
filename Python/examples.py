# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:39:31 2019

@authors: Joshua Kaggie, Brian Hargreaves
"""
import numpy as np
import matplotlib.pyplot as plt
import mrsigpy as mrs

from numpy import mean, std, exp, diag, matmul, pi, cos, sin, zeros, ones, shape, floor, ceil
from maplotlib.pyplot import figure, plot, xlabel, ylabel, title, imshow, subplot


def conceptB4_2(sigmean = 0):
    # -- Show distributions Gaussian, Rician (& Rayleigh)
    #
    
    N=100;					# histogram bins
    Ns=500;					# samples
    
    dat = np.random.randn(Ns,Ns)+1j*np.random.randn(Ns,Ns)+sigmean;	
    #dat = dat(:);
    [yreal,x] = plt.hist(np.real(dat),N);
    [ymag]= plt.hist(abs(dat),x);
    mn = np.mean(np.abs(dat));
    sg = np.std(np.abs(dat));
    plt.hold('off');
    plt.bar(x,[ymag, ' yreal'],1);
    plt.hold('on');
    print('%g*normpdf(x,%g,%g)' % sqrt(2*pi*sg*sg)*max(ymag),mn,sg);
    fplot(tt,[min(x),max(x)],':');
    
    legend('Rician/Rayleigh (|sig|)','Gaussian (Real{sig})','Gaussian w/ Rician stats');
    plt.grid('on')
    plt.title('Signal Distributions, mean=0');
    plt.xlabel('Signal');
    plt.ylabel('Frequency');
    
    mrs.setprops();    
    
    

def exampleA1_63():
    # Plot of Mz for IR Question.    
    t = np.arange(0, 1, .01)
    Mz = -1.63*np.exp(-t/0.5)+1;
    Mz[51:] = 1-np.exp(-(t[51:]-.5)/.5);
    plt.plot(t,Mz)
    mrs.lplot('Time (s)','M_z / M_0','Simple IR Signal Example',[0, 1, -1, 1]);
    mrs.setprops()


def exampleB1_13(TR = 1, TI = 0.5, TE = 0.05, T1 = 0.5, T2 = 0.1):
    #%	Short-TR IR Signal Example
    #%
    #%	Use abprop to do this compactly.
    #%	Note we start applying A,B from TE onward as this is steady-state.
    #%
    
    #% Simple Inversion-Recovery sequence where TR is not long enough for
    #% full recovery so steady-state must be calculated.
    
    #% Times are all in seconds here.    
    #% Shorthand - often describe E1 and E2 like this:
    
    E1a = exp(-(TI/T1));
    E2a = exp(-(TI/T2));
    E1b = exp(-(TE/T1));
    E2b = exp(-(TE/T2));
    E1c = exp(-((TR-TI-TE)/T1));
    E2c = exp(-((TR-TI-TE)/T2));
    
    #% Define A,B as in lecture example
    
    A1 = diag([E2a, E2a, E1a]) * mrs.xrot(180);
    B1 = np.array([[0],[0],[1-E1a]]);               
    A2 = diag([[E2b],[E2b],[ E1b]]) * mrs.xrot(90);
    B2 = [[0],[0],[1-E1b]];
    A3 = diag([E2c, E2c, E1c]);
    B3 = [[0],[0],[1-E1c]];
    
    A = matmul(A2, matmul(A1,A3));                        
    B = B2+matmul(A2,(B1+matmul(A1,B3)));
    
    Mss = matmul(np.linalg.inv(mp.eye(3)-A),B)
    return Mss      


def exampleB1_15(TR = 1, TI = 0.5, TE = 0.05, T1 = 0.5, T2 = 0.1):
    #%	Short-TR IR Signal Example
    #%
    #%	Use abprop to do this compactly.
    #%	Note we start applying A,B from TE onward as this is steady-state.
    #%
    
    #% Simple Inversion-Recovery sequence where TR is not long enough for
    #% full recovery so steady-state must be calculated.
    
    #% Times are all in seconds here.
    
    [A,B,Mss] = mrs.abprop(mrs.relax(TR-TE-TI,T1,T2,1),mrs.xrot(180), \\
    	mrs.relax(TI,T1,T2,1),mrs.xrot(90), mrs.relax(TE,T1,T2,1));
    
    return Mss	    

def exampleB1_17():
    pass

def exampleB1_19():
    pass

def exampleB2_2(D = 1e-6, G = 0.04, T = 10, gamma = 42.58, dt = 20):
    #% Calculation of simple 1D diffusion sensitivity.
    
    #D = 1e-6;	% mm^2/ms
    #G = 0.040;		% mT/mm;
    #T = 10;			% ms.
    #gamma = 42.58	% kHz/mT
    #dt = 20;	% ms
    
    
    sig = np.sqrt(2*D*dt);#	% mm.
    
    x = np.arange(-3, 3, 0.01)*sig #[-3:.01:3]*sig;	% positions.
    
    #% == Numerical Integration to test 
    var_sum = 0;
    dx = x(2)-x(1);
    for k in range(x): #=1:length(x);
      val = cos(2*pi*gamma*G*T*x[k]) * 1./sqrt(4*pi*D*dt) * exp(-x[k]**2/(4*D*dt));  #matmul?
      var_sum = var_sum+val*dx;

    print( sum)    
    
    
    #% == Calculated Solution
    sig = exp(-(2*pi*gamma*G*T)**2*dt * D)    
    return sig



def exampleB2_23(T=1, T2 = 9.4877, T1 = 19.4932, Q = [[0],[0],[1]]):
    #%	Simple spin echo example, with M = 0.81 at 1st spin echo.
    #%
    #%	Display state at each step (after trim)
    #%
    #T = 1;
    #T2 = 9.4877;			% exp(-2T/T2)=0.81
    #T1 = 19.4932;			% exp(-T/T1)=0.95
    #Q = [0;0;1]
    Q = mrs.epg_rf(Q,pi/2,pi/2)
    Q = mrs.epg_grelax(Q,T1,T2,T);
    Q = mrs.epg_trim(Q,.001)
    Q = mrs.epg_rf(Q,pi,0);
    Q = mrs.epg_trim(Q,.001)
    Q = mrs.epg_grelax(Q,T1,T2,T);
    Q = mrs.epg_trim(Q,.001)
    Q = mrs.epg_grelax(Q,T1,T2,T);
    Q = mrs.epg_trim(Q,.001)
    Q = mrs.epg_rf(Q,pi,0);
    Q = mrs.epg_trim(Q,.001)
    Q = mrs.epg_grelax(Q,T1,T2,T);
    Q = mrs.epg_trim(Q,.001)
    return Q


def exampleB2_43():
    #% Stimulated Echo Sequence (Bloch vs EPG) example
    #%
    #%	Slides contain matrix examples
    #%
    plt.figure(1);
    
    FZ = [[0],[0],[1]]#			% Equilibrium
    mrs.epg_show(FZ,0,1,23); #disp('0 - Press Enter'); drawnow; pause;
    FZ = mrs.epg_rf(FZ,pi/2,pi/2)	% After 90y
    epg_show(FZ,0,1,23);  disp('1 - Press Enter'); drawnow; pause;
    FZ = mrs.epg_grelax(FZ,1,1,0,0,0,1) % After Gradient
    epg_show(FZ,0,1,23); disp('2 - Press Enter'); drawnow; pause;
    M = mrs.epg_FZ2spins(FZ,12);	% Check!
    FZ = mrs.epg_rf(FZ,pi/2,pi/2)	% After 90y
    mrs.epg_show(FZ,0,1,23); disp('3 - Press Enter'); drawnow; pause;
    FZ = mrs.epg_grelax(FZ,1,1,0,0,0,1) % After Gradient
    mrs.epg_show(FZ,0,1,23); disp('4 - Press Enter'); drawnow; pause;
    FZ = mrs.epg_rf(FZ,pi/2,pi/2)	% After 90y
    mrs.epg_show(FZ,0,1,23); disp('5 - Press Enter'); drawnow; pause;
    FZ = mrs.epg_grelax(FZ,1,1,0,0,0,1) % After Gradient
    mrs.epg_show(FZ,0,1,23); disp('6 - Press Enter'); draabsolute_import

def exampleB2_5(G=40, T= 1.174, gamma = 42.58, alpha = np.arange(180,60,-1)):
    #% Effect of Crushers
    
    #G = 40;		% mT/m;
    #T = 1.174;	% ms.
    #gamma = 42.58	% kHz/mT
    #alpha = [180:-1:60];
    
    #% 2 cycles/mm.
    
    
    x = np.arange(0,1,0.01)/1000;	#% m.
    
    Grot = 360*x*G*T*gamma; 	#% Rotation due to gradient at each voxel.
    
    Ms=[[1],[0],[0]];		#% Mstart
    
    Mend = []
    for alpha_i in alpha:
     for Grot_i in Grot:
      M = Ms;
      M = zrot(Grot_i)*M;		#% Crusher 1
      M = xrot(alpha_i)*M;		#% Refocusing pulse.
      M = zrot(Grot_i)*M;		#% Crusher 2
      Mend.append(M) #  Mend(:,k)=M;  I think this is what you want
     #end;
    
     #%figure(3);
     #%plot(abs(Mend(1,:)+i*Mend(2,:)));
     figure(1);
     Mxy = Mend[0]+1j*Mend[1];
     plot(x*1000,real(Mxy),'k--'); hold on;
     axis([np.min(x)*1000, np.max(x)*1000, -1.2, 1.2]);
     plot([np.min(x), np.max(x)]*1000,[1 1]*np.abs(mean(Mxy)),'b-'); 
     #hold off;
     #grid on;
     xlabel('Position (mm)'); ylabel('Signal');
     legend('M_{xy}','Avg M_{xy}');
     tt = sprintf('%d Degree Refocusing Angle',alpha(n)); title(tt);
     setprops();
     drawnow;
     #%fig2tiff('crush',n);
     print(n)
     Mse(n) = np.abs(np.mean(Mxy));
    
     
    
    figure(2);
    plot(alpha,Mse); #grid on; 
    xlabel('Refocusing Angle (deg)');
    ylabel('Spin Echo Signal'); 
    title('Spin Echo vs Refoc. Angle');
    a = plt.gca(); 
    #axis([a(1:2) 0 1]);
    mrs.setprops();

    
    
def exampleB4_14(R=2, Na = 1000, cwid = 25):
    
    #%
    #%	Example 1D SENSE reconstruction with noise.
    #%
    
    #% -- 1. SETUP PARAMETERS	
    #R=2;				% Reduction Factor
    #Na = 1000;			% #Runs for propagated noise
    #cwid = 25;			% Gaussian coil sensitivity width
    
    #% Coil noise covariance elements 
    n12=.3; 
    n13=.13*1j;
    n14=.1;
    n23=.32;
    n24=.15;
    n34=.28;
    ndiag = diag([1, .9, 1.2, 1]);
    nscale = 0.1;#			% Noise scaling
    
    
    
    #% -- 2. IMAGE SHAPE 
    x = np.arange(-64, 64) #[-64:63];
    m =[zeros(10), 1j*np.arange(20), 1j*np.arange(20,0,-1), zeros(10),
    	sqrt(900-np.arange(-30,30)**2,  zeros(7)];#	% 2 shapes
    
    askip = floor(len(x)/R);
    x=x(1:askip*R);
    m=m(1:askip*R);
    
    #% -- 3. COIL SENSITVITIES (Gaussian)
    c1 = pdf('norm',x,-45,cwid);	% Coil 1
    c2 = pdf('norm',x,-15,cwid);	% Coil 2
    c3 = pdf('norm',x,15,cwid);	% Coil 3
    c4 = pdf('norm',x,45,cwid);	% Coil 4
    
    c4 = 0*c4;
    
    if (1==0)	% Sanity check!  Uniform coil
      c1 = ones(size(x))/40;
      c2 = 0*c1;
      c3 = 0*c2;
      c4 = 0*c3;
    end;
    
    c = 35*[c1(:) c2(:) c3(:) c4(:)]';	% Coil sensitivity matrix
    c = diag([1 1.1 0.9 1])*c;		% Different Scaling per coil
    
    
    #% -- 4. CHANNEL IMAGES and K-SPACE
    cim = [1;1;1;1]*m .*c;				% Rows of Image * Coil
    ksp = ifftshift(ifft(ifftshift(cim,2),[],2),2);	% k-space
    
    
    #% -- 5. DISPLAY STUFF
    
    figure(1);
    plotc(x,m);	lplot('x','signal','Image');
    
    figure(2);
    plot(x,c);	lplot('x','Coils','Coil Sensitivities');
    legend('Coil 1','Coil 2','Coil 3','Coil 4');
    
    figure(3);
    plot(x,abs(cim));	lplot('x','Coils','Coil Images');
    title('Noise Covariance Matrix \Psi');
    
    
    
    #% -- 6. GENERATE NOISE (Covariance, then Multivariate Gaussian)
    
    psi = np.arange([[0, n12, n13, n14],[ 0, 0, n23, n24],[ 0, 0, 0, n34],[ 0, 0, 0, 0]]);
    psi = (psi+np.transpose(psi))+ndiag;
    dispim(psi);
    n = nscale*(mvnrnd([0 0 0 0],psi,length(x)).' + i*mvnrnd([0 0 0 0],psi,length(x)).');
    ksp = ksp + n;
    
    
    #% -- 7. SENSE Reconstruction (SETUP)
    ims = 0*x;			% Allocate SENSE image
    gfact=ims;			% Allocate g-factor map
    gmc = gfact;			% Allocate
    gcond = gmc;			% Allocate
    inoise = gmc;			% Allocate, for image noise.
    
    askip = floor(length(x)/R);		% Amount to skip for aliased pixel(s).
    
    ipsi = inv(psi);		% Just do inverse-Psi once!
    kr = 0*ksp;			% Allocate reduced-k-space matrix
    kr(:,1:R:end) = ksp(:,1:R:end);	% Subsampled k-space
    
    ima = R*fftshift(fft(fftshift(kr,2),[],2),2);	% FT to aliased image data.
    Cpsave = zeros(R,4,length(x));			% Save Cp values
    
    #% -- 8. SENSE (Pixel-wise Recon)
    
    for xx = 1:askip	#% Loop across pixels
    
      C = c(:,xx:askip:end);		% Make Nc x R sensitivity matrix.
      CpC = C'*ipsi*C;			% First step
      iCpC = inv(CpC);			% Second step
      Cp = iCpC*C'*ipsi;			% P-inv Reconstruction matrix.
      Cpsave(:,:,xx) = Cp;			% Save for noise sim.
    
      dCpC = diag(CpC); 			% Diagonal elements
      diCpC = diag(iCpC);
      
      gfact(xx:askip:end) = sqrt(dCpC.' .* diCpC.');	% gfactor
      gmc(xx:askip:end) = sqrt(dCpC.');	% gfactor - due to coils
      gcond(xx:askip:end) = sqrt(diCpC.');	% gfactor - due to conditioning
    
      ims(xx:askip:end) = (Cp*ima(:,xx)).';	   % Reconstruct / unalias
      inoise(xx:askip:end) = sqrt(diCpC) * nscale * sqrt(length(x));	% Store calculated noise.
    
    
      %tt=sprintf('Pixel Group %d of %d',xx,askip); disp(tt);
    end;
    
    
    #% -- 9.  NOISE SIMULATION
    
    noiseim = zeros(Na,length(x));
    
    for k = 1:Na		% Noise simulation (per channel)
      %tt=sprintf('Noise run %d of %d',k,Na); disp(tt);
      npts = ceil(length(x)/R);
      cnoise = nscale*sqrt(length(x)) * ...
    	(mvnrnd([0 0 0 0],psi,npts).' + i*mvnrnd([0 0 0 0],psi,npts).');
      for xx = 1:askip	% Loop across pixels
        noiseim(k,xx:askip:end) = squeeze(Cpsave(:,:,xx))*cnoise(:,xx);
      end;
    end; 
    nmean = mean(noiseim);
    nstd = std(real(noiseim));
    
    
    
    figure(4);
    subplot(3,1,1);
    plot(x,abs(ima.'));	lplot('x','Coils','Aliased Images (Magnitude)'); 
    legend('Coil 1','Coil 2','Coil 3','Coil 4'); 
    subplot(3,1,2);
    plot(x,real(ima.'));	lplot('x','Coils','Aliased Images (Real)');
    legend('Coil 1','Coil 2','Coil 3','Coil 4'); 
    subplot(3,1,3);
    plot(x,imag(ima.'));	lplot('x','Coils','Aliased Images (Imaginary)'); 
    legend('Coil 1','Coil 2','Coil 3','Coil 4'); 
    setprops;
    
    figure(5);
    plot(x,real(ims),'b--',x,imag(ims),'r--',x,abs(ims),'k-');	
    legend('Real','Imag','Magnitude');
    lplot('x','Signal','SENSE Image');
    setprops;
    
    figure(6);
    plot(x,real(gfact),'k-',x,real(gmc),'r--',x,real(gcond),'b--');
    legend('g-factor','g_{multicoil}','g_{condition}');
    lplot('x','g-factor','g-factor across Image',[min(x) max(x) 0 3]);
    setprops;
    
     	 
    figure(7);
    plot(x,real(nmean),'b-',x,real(nstd),'r--');
    legend('Mean','Std dev');
    lplot('x','Noise mean / std-dev','Propagated Noise Statistics');
    setprops;
    
    figure(8);
    plot(x,real(inoise));
    a = axis; a(3)=-.5; 
    lplot('x','Noise','Calculated Noise sqrt(C Psi C)',a);
    title('Calculated Noise C^H \Psi C');
    setprops;

def exampleE1_vds():
    pass

def exampleE3_spiral():
    pass

def sense1d():
    pass
    


def examplenoise():
    #%
    #%	This example is to show the effects of noise in MRI.
    #%
    sigma = 0.1;		
    N = 256;
    
    #% A:  Gaussian noise in k-space
    #%
    knoise1 = randn(N,N) + i*randn(N,N);	#% Complex noise, sigma=1
    kn = knoise1*sigma;
    
    #% Show noise distribution
    figure(1); 
    #####ghist(real(kn(:)),[],[],[],'k-space real(noise)','k-space noise');
    
    #% B:  Real-Valued Noise in image space
    #%
    imn = nft(kn);			#% FFT, normalized by sqrt(N) to preserve energy
    isig = std(real(imn));	#% Get standard deviation to estimate sigma
    
    figure(2);
    ghist(real(imn),[],[],[],'Image real(noise)','image real noise');
    
    
    #% C:  Magnitude Noise in image space
    imag = abs(imn);		#% Magnitude of image.
    imsig = sqrt(2-pi/2)*isig;#	% Estimate Rayleigh sqrt(variance)
    immean = sqrt(pi/2)*isig;#	% Estimate Rayleigh mean
    
    figure(3); 
    #####ghist(imag,[],[],[0:0.01:1]*6*imsig,'Image Magnitude Noise','image abs(noise)');
    
    #% D:  Magnitude Noise+Signal in image space
    figure(4); 
    #####hold off;
    imsignal = 5*sigma;		#% SNR = 5;
    imag = abs(imn+imsignal);#	% Magnitude of image.
    imsig = std(imag);		#% Estimate of Rician sqrt(variance)
    immean = mean(imag);	#	% Estimate of Rician mean, to fit.
    
    figure(4); 
    ghist(imag,[],[],[],'Image abs(Signal+Noise)','abs(signal+noise)');
        
    im = zeros(N,N);		#% Allocate Image
    im[N/4+1:3*N/4,N/4+1:3*N/4]=1;	#% Central square with amplitude 1.
