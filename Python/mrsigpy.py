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
# authors: Joshua Kaggie, Brian Hargreaves
# -----------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

PLOTON = True


from time import sleep





def animation_plot():
    #%matplotlib notebook  ## RUN THIS LINE IF IN JUPYTER
    import matplotlib.animation
    t = np.linspace(0,2*np.pi)
    x = np.sin(t)
    
    fig, ax = plt.subplots()
    ax.axis([0,2*np.pi,-1,1])
    l, = ax.plot([],[])
    
    def animate(i):
        l.set_data(t[:i], x[:i])
    
    ani = matplotlib.animation.FuncAnimation(fig, animate, frames=len(t))
    
    from IPython.display import HTML
    HTML(ani.to_jshtml())






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


def _3dspins():
    # This import registers the 3D projection, but is otherwise unused.
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import    
    import matplotlib.pyplot as plt
    import numpy as np    
    fig = plt.figure()
    ax = fig.gca(projection='3d')    
    # Make the grid
    x, y, z = np.meshgrid(np.arange(-1, 1, 0.1),
                          np.arange(-1, 1, 0.1),
                          np.arange(0, 1, 1))
        
    # Make the direction data for the arrows
    u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
         np.sin(np.pi * z))
    
    ax.quiver(x*0, y*0, z*0, u, v, w, length=0.1, normalize=True)
    
    plt.show()

# Converts a time vector to the corresponding frequency vector of an FFT
def time2freq(t):
    dt = t[1]-t[0]
    num_p = len(t)
    df = 1./ num_p/dt
    f = np.arange(np.ceil((num_p-1.)/2.), np.floor((num_p-1.)/2.))*df
    return f


import numpy as np
from matplotlib.pyplot import figure, plot, subplot, title, xlabel, ylabel

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


def adiabatic(peakb1 = 0.2, bw = 2000, beta = 1000., T = 0.01, Ts = 0.00001, blochsim = True):
    T = 2*np.round(T/Ts/2)*Ts
    N = T/Ts
    
    t = np.arange(Ts, T, Ts) - np.round(N/2)*Ts #time from -t/2 to t/2
    
    b1 = peakb1 * np.sech(beta*t)
    freq = bw/2 * np.tanh(beta*t)
    phase = np.cumsum(freq)*2*np.pi*Ts
    phase = phase-phase(np.round(N/2))  #zero phase half-way
    phase = np.mod(phase+np.pi,2*np.pi)-np.pi  #limit to -pi to pi
    
    if blochsim:
        t = t-t[0]  #start t at 0 for plots
        figure(1)
        print('Adiabatic silver-hoult pulse (beta = '+str(beta)+' Hz) ')
        subplot(3,1,1)
        plot(t,b1); xlable('Time(s)'); ylabel('B1(G)')
        title(tt)
        subplot(3,1,2)
        plot(t,phase); xlabel('Time(s)'); ylabel('Phase(rad)')
        subplot(3,1,3)
        plot(t,freq); xlabel('Time(s)'); ylabel('Freq(hz)')
        
        #add in 
        if 0:
            gr = 0*b1
            tp = Ts
            t1 = 0.6; t2 = 0.1
            df = np.arange(-3*bw, 3*bw, bw/20.)
            dp = 0
            mode = 0
            mx,my,mz = bloch(0)  # finish this
            
            
            
def bvalue(gradwave, T):
    gamma = 2*np.pi*42.58
    intg = np.cumsum(gradwave)*T
    b = gamma**2 * np.sum(intg**2)*T
    return b

def calcgradinfo(g,T=0.000004,k0=0,R=0.35,L=0.0014,eta=1/56, gamma = 4258):
    s = np.shape(g)    
    lg = np.max(np.size(g))
    k = k0+np.cumsum(g)*gamma*T
    t = T*(np.arange(len(g))-0.5)
    t = np.transpose(t)
    tt = t*np.ones(s[1])
    s = [[g],[g[lg]]]-[[0*g[0]],[g]]
    sz = np.shape(s)
    s = s[1:sz[0]]/T
    m1 = np.cumsum(g*tt)*gamma*T
    m2 = np.cumsum(g*(tt*tt+T**2/12))**gamma*T
    v = (1/eta)*(L*s+R*g)
    return k,g,s,m1,m2,t,v


# Use np.multivariate ???
# !!!!
_default_R = [np.arange(0.2,1,.2),np.arange(0.1,0.8,.2),np.arange(-0.2,0.6,.2)]
def corrnoise(mn=None,R=_default_R,n=100):
    if mn is None:
        mn = np.zeros(np.shape(R[0])[0])
    nvect = np.random.randn(np.shape(R)[0],n)
    R = 0.5 * (R+ np.tranpose(R))
    v,d = np.linalg.eig(R)
    
    if np.any(np.diag <=0):
        print('R must be positive definite')
    
    w = v * np.sqrt(d)
    nc = w*nvect
    Rout = (nc*np.tranpose(nc))/n        
    return nc,Rout


def cropim(im,sx=None,sy=None):
    sz = np.shape(im)
    if sx is None: sx = np.floor(sz[0]/2)
    if sy is None: sy = sx
    if sx < 1: sx = sz[0]*sx
    if sy < 1: sy = sz[0]*sy
        
    stx = np.floor(sz[1]/2-sx/2)+1
    sty = np.floor(sz[0]/2-sy/2)+1        
    
    return im[sty:sty+sy-1, stx:stx+sx-1]
        


def csens2d(klow):
    sz = np.shape(klow)
    nc = sz[2]
    
    cs = np.fft.fftshift(np.fft.fft(np.fft.fftshift(klow,axis=0),axis=0),axis=0)
    cs = np.fft.fftshift(np.fft.fft(np.fft.fftshift(cs,axis=1),axis=1),axis=1)
    cmrms = np.sqrt(np.sum(np.conj(cs)*cs,axis=-1))
    meancm = np.mean(cmrms)    
    
    f = np.argwhere(cmrms == 0)
    cs = np.reshape(cs,sz[0]*sz[1],nc)
    cmrms[f] = meancm/100
    
    cs[f] = meancm/10000
    ncs = cs / (cmrms*np.ones[0,nc])
    ncs = np.reshape(ncs,sz)
    cs = np.reshape(cs,sz)
    
    return ncs, cs, cmrms
        

def dispangle(arr):
    angarr = np.angle(arr)+np.pi
    dispim(angarr,0,2*np.pi)
    
    
def dispim(im,low=0.0,high = None):
    im = np.squeeze(im)
    
    if high is None:
        immax = np.max(np.abs(im))
        imstd = np.std(np.abs(im))
        high = immax - 0.5 * imstd
        
    scale = 256.0/(high-low)
    offset = scale*low
    im = scale * im - offset 
    plt.figure
    plt.imshow(np.abs(im),cmap='gray')
    plt.axis('square')
    plt.show
    
    
def displogim(im):
    im = np.squeeze(im)
    lowexp = 0
    im = np.log(np.abs(im))
    
    f = np.where(im<lowexp)
    im[f] = lowexp
    
    im = im - lowexp
    dispim(im)

def dispkspim(ksp = None):
# Display the (log) magnitude of k-space and the phase, then inverse FFT and
# display the magnitude and phase of the image.
  
    if ksp is None:
        ksp = zeros((256,256))      
        ksp[128,128]=1
        ksp[127,128]=0.5;
        ksp[129,128]=0.5;

    im = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ksp)))
    logksp = 10*np.log(abs(ksp)+1)
    khigh = np.max(np.abs(logksp))-0.5*np.std(np.abs(logksp))
    imhigh = np.max(np.abs(im))-0.5*np.std(np.abs(im))

    fig, axes = plt.subplots(2, 2,figsize=(16,16))
    axes[0,0].imshow(logksp,cmap="gray",vmin=0,vmax=khigh)
    axes[0,0].set_title('k-space log-magnitude')
    axes[0,0].axis("off")
    axes[0,1].imshow(np.angle(ksp),cmap="gray")
    axes[0,1].set_title('k-space phase')
    axes[0,1].axis("off")
    axes[1,0].imshow(np.abs(im),cmap="gray",vmin=0,vmax=imhigh)
    axes[1,0].set_title('Image Magnitude')
    axes[1,0].axis("off")
    axes[1,1].imshow(np.angle(im),cmap="gray")
    axes[1,1].set_title('Image Phase')
    axes[1,1].axis("off")

    plt.show
    return im


def epg_cpmg(flipangle = [np.pi/2,np.pi/2,np.pi/2], etl = None, T1 = 4, T2=.1, esp = None, plot = False):

    if etl is None: etl = len(flipangle)
    if esp is None: esp = len(flipangle)
        
    etl = int(etl)
    esp = int(esp)

    if len(flipangle)==1 and etl>1 and np.abs(flipangle).all()<np.pi:
        flipangle[1] = flipangle[0]
        flipangle[0] = np.pi*np.exp(1j*angle(flipangle[1])+np.flipangle[1])/2
        
    P = np.zeros((3,2*etl))
    P[2,0] = 1
    Pstore = np.zeros((4*etl, etl))
    Zstore = np.zeros((2*etl, etl))
    
    P = epg_rf(P,np.pi/2, np.pi/2)
    s = np.zeros(1,etl)
    
    for ech in np.arange(etl):
        P = epg_grelax(P,T1,T2,esp//2,1,0,1,1)
        P = epg_rf(P,abs(flipangle(ech)),angle(flipangle(ech)))
        P = epg_grelax(P,T1,T2, eps//2,1,0,1,1)
        
        s[ech] = P[0,0]
        Pstore[2*etl:r*etl,ech] = P[1]
        Pstore[:2*etl,ech]=np.flipud(P[0])
        Zstore[:,ech] = P(2)  #  what does .' mean in matlab??
        
    if plot:
        pass
        #room to plot things
    
    return s,phasediag,P




def epg_relax(FpFmZ, T1, T2, T):
    if (T1 < 0) or T2<0 or T<0:
        print('Your values should not be negative...  Are you a time-traveller?')
    E2 = np.exp(-T/1./T2)
    E1 = np.exp(-T/1./T1)
    
    EE = np.diag([E2,E2,E1])
    RR = 1-E1
    
    FpFmZ = np.matmul(EE,FpFmZ)
    FpFmZ[2,0] = FpFmZ[2,0]+RR
    return FpFmZ
    


# Convert from M=[[mx],[my],[mz]] (3xN) to EPG state coefficient matrix 
def epg_spins2FZ(M = [[0],[0],[1]],trim=0.01):

    #!!! Need to check M has 3 rows
    N = np.shape(M)[1] 		# Size of M
    Q= np.int(np.floor(N/2)+1)	# Max Number of columns of FZ

    # -- The following are Eqs. 2-4 from Wiegel 2010:
    # -- Note you could do ONE FFT for Fp and Fm instead, if
    # -- speed is an issue.
    M = np.fft.ifftshift(M,axes=1)	
    Mxy = M[:1,:]+ 1j*M[1:2,:]		     # Mx+jMy, to be clear
    #print("Mxy is %s" % Mxy)
    Fp = np.fft.fft(Mxy,axis=1)/N   	     # FFT to F+ 
    #Fm = np.fft.fft(np.conj(Mxy),axis=1)/N   # Could do this way...
    Z  = np.fft.fft(M[2:,:],axis=1)/N        # FFT to F+ states.

    # -- Fm coefficients from right-half of Fp, truncate before fliplr
    Fmc = np.fliplr(np.conj(np.roll(Fp[:1,:],-1,axis=1)[:1,Q-1:]))
    #print("Fp is %s" % Fp)
    #print("Fmc is %s" % Fmc)
    #print("Fm is %s" % Fm)

    # !!! Could define Fmm from right half of Fp to check

    FpFmZ = np.concatenate((Fp[:,:Q],Fmc,Z[:,:Q]),axis=0)  # Combine to 3xQ 
    #print("FpFmZ is %s" % FpFmZ)
     
    FpFmZ = epg_trim(FpFmZ,trim)            # Trim near-zero states.
    return FpFmZ




# Convert from EPG state coefficient matrix to M=[[mx],[my],[mz]] (3xN)
def epg_FZ2spins(FpFmZ = [[0],[0],[1]],N=None,frac  = 0):

    Ns = np.shape(FpFmZ)[1]-1
    if N is None: N = 2.*Ns-1

    # -- Make Fourier matrix - there may be a cleaner way!
    # -- positions effectively give fftshift after transform.
    z = (np.arange(0,N)-np.floor(N/2))/N
    #z = (np.arange(N).astype(np.float)+0.5)/N-0.5   # z = fraction across voxel
    z = np.expand_dims(z, axis=0)
    phz = (1j*2.*np.pi*z)                   # 2pi*i*z) 
    inds = np.arange(-(Ns),(Ns+1))+frac	    # State number n=[-Ns:Ns]
    inds = np.expand_dims(inds,axis=1)
    fmat = np.exp(np.matmul(inds,phz))      # Fourier matrix
    #print("fmat = %s" % fmat)

    # -- Combine Fminus and Fplus states (or Fn for n<0 and n>=0) 
    Fstates = np.hstack([np.fliplr(np.conj(FpFmZ[1:2,1:])), FpFmZ[0:1,0:]])

    # -- Get last row (Z states) 
    Zstates = 2.*FpFmZ[2:,:]		# Double to account for one-sided FT.
    Zstates[0,0]=Zstates[0,0]/2		# Don't double Z0
    #print("Z states (doubled for n>0) = %s" % Zstates)
    
    # -- Fourier transform to Mxy
    Mxy = np.matmul(Fstates,fmat)
    #print("Mxy is %s" % Mxy)

    # -- Fourier transform to Mz
    fmatz = fmat[Ns:,:]
    #print("fmatz = %s" % fmatz)
    Mz = np.real(np.matmul(Zstates,fmat[Ns:,:]))
    #print("Mz is %s " % Mz)

    # -- Extract Mx, My and Mz and return as 3xN.
    spins = np.concatenate((np.real(Mxy),np.imag(Mxy),Mz),axis=0)

    return spins

  
    


def epg_grad(FpFmZ=[[1],[1],[0]], noadd=0, positive = True):
    if not noadd:
        FpFmZ = np.hstack([FpFmZ, [[0],[0],[0]]])  
    if positive:
        FpFmZ[0][1:] = FpFmZ[0][:-1]
        FpFmZ[1][:-1] = FpFmZ[1][1:]
        FpFmZ[1][-1] = 0;
        FpFmZ[0,0] = np.conj(FpFmZ[1,0])
    else:
        FpFmZ[1][1:] = FpFmZ[1][:-1]
        FpFmZ[0][:-1] = FpFmZ[0][1:]
        FpFmZ[0][-1] = 0;
        FpFmZ[1,0] = np.conj(FpFmZ[0,0])
        
    return FpFmZ

def epg_mgrad(*kw, **kws):
    "Negative gradients"
    return epg_grad(positive = False, *kw, **kws)    


# Discard states higher than highest with a coefficient greater than thres
def epg_trim(FpFmZ, thres):
   
    a = np.max(np.argwhere(np.abs(FpFmZ)>thres),axis=0)[1]
    FpFmZ = FpFmZ[:,:a+1]
    return FpFmZ


def epg_rf(FpFmZ = [[0],[0],[1]], alpha = 90.,phi = 90, in_degs = True, return_rotation = False, 
           frames = 1, ploton = PLOTON):
    if frames > 1:
        for xi in range(frames):
            FpFmZ = epg_rf(FpFmZ, alpha/frames , phi, in_degs, return_rotation=False, frames=1, ploton = False)        
            epg_show(FpFmZ)
        return 0
        
    if in_degs:
        alpha = alpha*np.pi/180.
        phi = phi*np.pi/180.
        
    RR = [[(np.cos(alpha/2.))**2., np.exp(2.*1j*phi)*(np.sin(alpha/2.))**2., -1j*np.exp(1j*phi)*np.sin(alpha)],
      [np.exp(-2.*1j*phi)*(np.sin(alpha/2.))**2., (np.cos(alpha/2.))**2., 1j*np.exp(-1j*phi)*np.sin(alpha)],
      [-1j/2.*np.exp(-1j*phi)*np.sin(alpha), 1j/2.*np.exp(1j*phi)*np.sin(alpha),      np.cos(alpha)]];
    
    if return_rotation:
        return RR
    else:
        return np.matmul(RR,FpFmZ)



    
    



def epg_stim_calc(flips, in_degs = True):
    P = [[0],[0],[1]]
    
    if in_degs:
        flips = np.array(flips)*np.pi/180.
    
    for flip in flips:
        P = epg_rf(P, flips,np.pi/2, in_degs = in_degs)        
        P = epg_grelax(P,1,.2,0,1,0,1)
    S = P[0,0]
    return S,P



def ft(dat):
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(dat)))

def gaussian(x,mnx,sigx,y=None,mny=None,sigy=None):
    if (y is None):
      return np.exp(-(x-mnx)**2/(2*sigx**2))/np.sqrt(2*np.pi)/sigx
    else:
      xarr = np.exp(-(x-mnx)**2/(2*sigx**2))/np.sqrt(2*np.pi)/sigx
      xarr = np.expand_dims(xarr,axis=1)
      yarr = np.exp(-(y-mny)**2/(2*sigy**2))/np.sqrt(2*np.pi)/sigy
      yarr = np.expand_dims(yarr,axis=0)
      return np.matmul(xarr,yarr)




#       Function makes a histogram of a signal, and superimposes
#       a gaussian fit.
#
#       INPUT:
#               data = signal to make histogram for.
#               gmean = estimate for gaussian fit, mean [default=mean(data)]
#               gsig = estimate for gaussian fit, sigma [default=std(data)]
#               nbins = number of bins
#               gtitle = plot title
#               gleg1 = legend top line (describe data)
#               gleg2 = legend second line (describe gaussian)
#
def ghist(data,gmean = None,gsig = None,bins = 100,gtitle='Histogram and Gaussian Fit',gleg1='Data',gleg2='Gaussian'):

  if (gmean==None):
    gmean = np.mean(data)
  if (gsig==None):
    gsig = np.std(data)
  if (np.size(bins) == 2):
    bins = np.linspace(bins[0],bins[1],100)
  if (np.size(bins) < 2):
    nbins = bins		# Passed number of bins
    bins = np.linspace(gmean-3*gsig,gmean+3*gsig,nbins)

  (hist,binp) = np.histogram(data,bins,(np.min(bins),np.max(bins)))
  binp = (binp[:-1]+binp[1:])/2.0
  gfit = gaussian(binp,gmean,gsig) 	 	# Gaussian fit
  gscale = np.size(data)*(binp[1]-binp[0]);	# Scale
 
  gleg2 = gleg2 + " mn: %4g, sd: %4g" % (gmean,gsig)

  # -- Debugging! 
  #print("Data: %s" % data)
  #print("Bins: %s" % bins)
  #print("Hist: " , hist)
  #sd = np.shape(data)
  #sb = np.shape(bins)
  #sh = np.shape(hist)
  #print('Data shape,',sd)
  #print('Bins shape,',sb)
  #print('Hist shape,',sh)

  plt.plot(binp,hist,label=gleg1)
  plt.plot(binp,gfit*gscale,label=gleg2)
  plt.legend()
  lplot('value','frequency',gtitle)

  return hist

# RMS coil noise simulation
#
# Simulate the signal plus Gaussian noise, with RMS combination for
# multiple coils (or just one!) 
#
# INPUT:
#	sig = signal (Ns x 1)
#       cov = coil covariance matrix (Nc x Nc)
#	csens = coil sensitivites (Ns x Nc, or 1 x Nc if uniform)
#       Nc = number of coils (used if cov AND csens not passed
#	plotfig = figure object - if given, plots bias and stdev vs signal.
#	plottitle = title for plot (optional)
#
# RETURNS:
#       sigdist = signal distribution with noise
#	bias (mean signal difference) and standard deviation
#
# EXAMPLES:
#  	rmssig = mrs.rmscoilnoise(Nc=1)  [Simple Rayleigh/Rician noise]
#       rmssig = mrs.rmscoilnoise(csens=np.array([[1,1,1]]))  [3 coils]
#       rmssig = mrs.rmscoilnoise(cov=np.eye(3))  [3 coils, again!]

def rmscoilnoise(sig=None,cov=None,csens=None,Nn=10000,Nc=1,scale=None,plotfig=None,plottitle=None):

  if (len(sig)<2 and sig==None):   # -- If no sig, make a ramp
    sig = np.arange(0.0,10.0,0.1)	


  # -- Define csens,cov and Nc
  #    --Tricky, as user can pass any one of them!
  if (csens is None and cov is None):	# Define from Nc
    csens = np.ones((1,Nc))
  else:
    if (csens is None):
      Nc = np.shape(cov)[0]		# Assume cov is square!
    else:
      Nc = np.shape(csens)[1]

  if (cov is None):			# No coil covariance specified, use Nc
    cov = np.eye(Nc)
  else:
    Nc = np.max(np.shape(cov))	# Otherwise set Nc based on cov

  if (csens is None):		# -- Default coil sensitivities to all 1
    csens = np.ones((1,Nc))


  if (np.shape(csens)[1]!=Nc):
    print("Coil Sensitivity Matrix should be consistent with Nc and/or cov")
    print("csens is ",np.shape(csens))
    print("Number of coils is ",Nc)
    print("Covariance is ",np.shape(cov))

# -- Generate gaussian noise,
# !!! should use multivariate here.

  cnoise = np.random.multivariate_normal(np.zeros(Nc),cov,(1,Nn)) + 1j*np.random.multivariate_normal(np.zeros(Nc),cov,(1,Nn))
  cnoise = np.moveaxis(cnoise,(0,1,2),(0,2,1))
  #cnoise = np.random.normal(0.0,1.0,(1,Nc,Nn)) + 1j*np.random.normal(0.0,1.0,(1,Nc,Nn))  # Complex gaussian noise sig=1

  #print("cnoise shape is ",np.shape(cnoise))
  #print("cnoise (pt 1) is ", np.squeeze(cnoise[0,:,0]))

  sz = np.shape(sig)
  csig = np.matmul(sig,csens)  # -- Coil signals, Spread by sensitivities
  csig = np.expand_dims(csig,axis=2)    # -- Add dimension for noise samples
  cnsig = csig+cnoise          # -- Generate complex signal distributions (broadcast)
  #print("sig+cnoise (pt 1) is ", np.squeeze(cnsig[0,:,0]))
  # !!! noise is same at each sig location - not okay later!!!
  
  # -- Normalize by RMS of sensitivity.  This is done to measure the 
  # -- signal bias (otherwise the signal scales oddly.  Here the noise should
  # -- be higher with more coils.
  rmssens = np.sqrt(np.sum(csens * np.conj(csens))) 
  rmssig = np.sqrt(np.sum(cnsig * np.conj(cnsig),axis=1))/rmssens   
  #print("RMS sensitivity is ",rmssens)
  #print("RMS Signal Shape ",np.shape(rmssig))
 
  # -- Normalize by sqrt(Ncoils) - more common for actual recon!
  #rmssig = np.sqrt(np.sum(cnsig * np.conj(cnsig),axis=1)/Nc)   

  # -- Calculate signal standard deviation, bias 
  sigstd = np.std(rmssig,axis=1)
  sigbias = np.squeeze(np.mean(rmssig,axis=1)) - np.squeeze(sig)
  #print("Signal std Shape ",np.shape(sigstd))

  if (plotfig is not None):
    plotfig.plot(sig,sigbias,'b--',label='Bias')
    plotfig.plot(sig,sigstd,'r-',label='Std Dev')
    if (plottitle is None):
      plottitle = "Bias/Std.Dev"
    lplot('Signal','Bias, Std.Dev',plottitle)	# lplot labels plot.
    plotfig.legend()
 
  return (rmssig,sigbias,sigstd)

# Calculate the SENSE combination weights for coils
#
# INPUT:
#	coilsens = "Npix" x Nc x R array - !NOTE for R=1, last dimension is 1
#  	noisecov = Noise covariance matrix (Nc x Nc), defaults to Identity
#	gfactor = True to calculate g factor
#
#	NOTE:  "Npix" can be 1D, 2D, or 3D...
#
# OUTPUT:
#       SENSE weights matrix (Npix x R x Nc)   [U in Pruessmann paper]
# 	g factor ("Npix" x R) - Calculated g-factor at each pixel.
#       calcnoise ("Npix" x R) - Calculated noise at each pixel, scaled by
#		sqrt(R) which accounts for scan-time reduction
#         NOTE:  for g factor and calcnoise, the calling code knows 
#		which image dimension(s) is/are aliased so can rearrange
#		to plot.
def senseweights(coilsens, noisecov=None,gfactorcalc=False, noisecalc=False):
  
  cs = coilsens.shape			# Get shape
  imshape = cs[0:-2]			# Shape of image (1D, 2D or 3D)

  R = coilsens.shape[-1]
  Nc = coilsens.shape[-2]
  Npts = np.int(np.prod(coilsens.shape)/R/Nc) # -- Number of pixels

  # -- Reshape to pixels x coils x aliased pixels
  print("SENSE Weight Calc - %d pts, %d coils, R=%d" % (Npts,Nc,R))
  coilsens = np.reshape(coilsens,(Npts,Nc,R),order='F')
  #print("Reshaped sensitivity size",coilsens.shape)

  # -- Allocate output arrays (weights and g factor)
  senseout = np.zeros((Npts,R,Nc),'complex64')
  if gfactorcalc is True:
    gfactor = np.zeros(R)
    gfactorout = np.zeros((Npts,R))

  if noisecalc is True:
    ncalc = np.zeros(R)
    ncalcout = np.zeros((Npts,R))

  # -- Default noise covariance to identity matrix 
  if (noisecov is None):
    noisecov = np.eye(Nc)
  psiI = np.linalg.inv(noisecov)

  # -- Loop over Pixels 
  for i in range(Npts):

    C = coilsens[i,:,:]		# -- Get Nc x R sensitivity for pixel

    # -- Multiply terms, minimizing # multiplications 
    # -- Equation [2] from Pruessmann SENSE paper

    Cc = np.conj(C.T)
    ChPsiI = np.matmul(Cc,psiI)	
    ChPsiIC = np.matmul(ChPsiI, C)	
    ChPsiICI = np.linalg.inv(ChPsiIC)
    sensemat = np.matmul( ChPsiICI , ChPsiI )
    senseout[i,:,:] = sensemat

    # -- g-factor calculation at pixel
    # -- Pruessmann et al. Eq. 23
    if gfactorcalc is True:
      for g in range(0,R):	
        gfactor[g] = np.sqrt( np.real(ChPsiICI[g,g] * ChPsiIC[g,g] ))
   
        if (False):  #-- Printout for debuggin
          print("  g-factor term 1 is ",np.sqrt(ChPsiICI[g,g]))
          print("  g-factor term 2 is ",np.sqrt(ChPsiIC[g,g]))
          print("Overall gfactor is ",gfactor[g]);
      gfactorout[i,:] = gfactor

    if noisecalc is True:
      ncalc = np.sqrt( np.real(np.diag(ChPsiICI) ))
      ncalcout[i,:] = ncalc


 
    if (False):		# -- Print for debugging
      print("For pixel %d:" % i)
      print("SENSE weight matrix is ",sensemat)
      print("g-factor is ",gfactor)
      
 
  # -- Reshape to pixels
  outputshape = imshape + (R,) + (Nc,)		# !!!better way?
  senseout = np.reshape(senseout,outputshape,order='F')
  
  if gfactorcalc is True:
    outputshape = imshape + (R,) 		# !!!better way?
    gfactorout = np.reshape(gfactorout,outputshape,order='F')
    if noisecalc is True:
      ncalcout = np.sqrt(R) * np.reshape(ncalcout,outputshape,order='F')
      return(senseout,gfactorout,ncalcout)
    else:
      return(senseout,gfactorout)

  else:
    return senseout





# Apply SENSE combination to data.
#
# INPUT:
#	signal = Npix x Nc array of received signals for FOV/R.
#	sweights = SENSE weights, Npix x R x Nc array of SENSE weights.
#
#       "Npix" can be nxm or nxmxp for 2D or 3D.

# OUTPUT:
#       spix = Npix x R array of reconstructed pixels
#
#	It is up to the calling code to sort out the unaliased pixels...

def senserecon(signal,sweights):

  ws = sweights.shape
  imshape = ws[0:-2]
  Nc = ws[-1]
  R = ws[-2]
  Npts = np.int(np.prod(ws)/Nc/R)

  ss = signal.shape

  # -- Check Consistency
  if (ss[-1] != Nc):
    print("Error (senserecon):  inconsistent signal and weights matrices.")
    print("  signal size is ",ss)
    print("  weights size is ",ws)
    print("  signal should be pix x channels, weights pix x R x channels")

  signal = np.reshape(signal,(Npts,Nc),order='F')
  sweights = np.reshape(sweights,(Npts,R,Nc),order='F')
  
  # -- Allocate output
  #print("Output size will be %d x %d" % (Npts,R))
  sigout = np.zeros((Npts,R),dtype='complex64') 
  for i in range(Npts):
    sigi = signal[i,:]
    wi = sweights[i,:,:]
    sigout[i,:] = np.matmul(wi,sigi)
  
  # -- Reshape to pixels
  outputshape = imshape + (R,) 		# !!!better way?
  sigout = np.reshape(sigout,outputshape,order='F')

  return(sigout)


# RMS Coil Combination
#
# INPUT:
#	csignals = Npix x Nc array of pixels
#	csens = Npix x Nc array of coil sensitivities 
# OUTPUT:
#  	RMS combination.
def rmscombine(csignals,csens=None):

  if (csens is None):		#!! better way?  None doesn't work for arrays
    rmssig = np.sqrt(np.sum(csignals * np.conj(csignals), axis=-1 ))
  else:
    csens = np.squeeze(csens)
    rmssens = np.sqrt(np.sum(csens * np.conj(csens),axis = -1))
    #print("Size of rmssens is",rmssens.shape)
    rmssig = np.sqrt(np.sum(csignals * np.conj(csignals),axis=-1))/rmssens

  return rmssig

    





def gridmat():
    pass

def homodyneshow(im,w=16):
    m,n = np.shape(im)
    lpf = np.zeros(m)
    lpf[np.floor(m/2)-w/2:np.floor(m/2)+w/2-1] =1
    
    mf = np.cumsum(lpf)
    mf = 2*mf/np.max(mf)
    subplot(3,2,1)
    plot(lpf)
    subplot(3,2,2)
    plot(mf)
    ksp = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(im)))
    klp = np.diag(lpf)*ksp
    imlp = ft(klp)
    subplot(3,2,3)
    dispim(imlp)
    plt.axis('off')
    title('low res image')
    subplot(3,2,4)
    dispangle(imlp)
    plt.axis('off')
    title('phase')
    
    khf = np.diag(mf)*ksp
    imhf = ft(khf)
    subplot(3,2,5)
    dispim(imhf)
    title('zero filled')
    subplot(3,2,6)
    dispangle(imhf)
    title('phase')
    
    phest =  np.angle(imlp)
    imhd = imhf * np.exp(-1j*phest)
    return imhd

def ksquare(center=0, swidth=1.9, kloc=None, tsamp=0.000004, df=0):
    if kloc is None:
        kx,ky = np.meshgrid(np.arange(-128,128)/128*5, np.arange(-128,128)/128*5)
        kloc = kx*i*ky
    sdata = swidth* np.sinc(swidth*np.real(kloc))*np.sinc(swidth*np.imag(kloc))
    kdata = 0*sdata
    
    for q in np.arange(len(center)):
        thisk = np.exp(1j*2*np.pi*np.real(center(q))*np.real(kloc))*sdata
        thisk = np.exp(1j*2*np.pi*np.imag(center(q))*np.imag(kloc))*thisk
        kdata = kdata + thisk
        
    ph = np.exp(2*1j*np.pi*tsamp*np.arange(len(kloc))*df)
    kdata = np.diag(ph)*kdata
    return kdata
        

        
def lfphase(m,n = None,w = 4):
    if n is None:
        n=m
    arr = np.zeros((m,n))
    arr[np.floor[m/2-w]:np.floor[m/2-w],np.floor[n/2-w]:np.floor[n/2-w]] = np.random.randn((2*w+1,2*w+1))
    ang = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(arr)))
    ang = np.real(ang)
    ang = 2*np.pi(ang-np.min(ang))/(np.max(ang)-np.min(ang)-np.pi)  # watch the dimensions..
    ph = np.exp(1j*ang)
    return ph

# Label plot with x,y axis labels, and title.  Turn on grid by default
#
# INPUT:
#	xlab = x-axis label
#	ylab = y-axis label
#	ptitle = plot title
#	ax = plot axis limits (4-tuple)
#	gridon = True to turn on grid
#	legend = True to turn on legend
#
def lplot(xlab = None,ylab = None, ptitle = None,ax = None,gridon = True,legend=False):
    if xlab is not None:
      xlabel(xlab)
    if ylab is not None:
      ylabel(ylab)
    if ptitle is not None:
      title(ptitle)
    if ax is not None:
      plt.axis(ax)
    if (gridon is True):
      plt.grid(gridon)
    if (legend is True):
      plt.legend()



                                

def lsfatwater():
    pass


def magphase(x,arr):
    mag = np.abs(arr)
    phase = np.angle(arr)
    
    fig1 = plt.subplot(2,1,1)
    plt.plot(x,phase/np.pi)
   

# epg_showstate(ax,FZ,frac=0,Nspins=19,voxvar=0):
#
# Plots a set of vectors on a single axis, from EPG coefficients.
#
# INPUTS:
#	ax = plot figure/axis upon which to plot
#	FZ = EPG coefficient matrix to display
#	frac = Additional dephasing.  Essentially adds to n for F_n 
#              states to show the dephasing/rephasing during gradient.
#	Nspins = Number of spins to use in graphical displays
#	voxvar = 1 to vary origin of vectors along z, 2 along x and
#		 0 all from (0,0,0).  Often use 1 for F_n and 2 for Z_n.
#
# OUTPUT:
#	M = corresponding magnetization endpoints ([[Mx],[My],[Mz]])
	
def epg_showstate(ax,FZ,frac=0,Nspins=19,voxvar=0):
   
  M = epg_FZ2spins(FZ,Nspins)
  scale = 1.

  mx = M[0:1,:]
  my = M[1:2,:]
  mz = M[2:3,:]
  x = np.zeros((1,Nspins))
  y = np.zeros((1,Nspins))
  z = np.zeros((1,Nspins))

  if (voxvar==1):
    z = 2*(np.arange(0,Nspins)-np.floor(Nspins/2))/Nspins  #2x fills axis!
  if (voxvar==2):
    x = 2*(np.arange(0,Nspins)-np.floor(Nspins/2))/Nspins
    #x = 2*(np.arange(0,Nspins)-np.float(Nspins)/2.+0.5)/Nspins

  ax.quiver(x, y, z, mx,my,mz,normalize=False)
  ax.set(xlim=(-scale,scale),ylim=(-scale,scale),zlim=(-scale,scale))


  # -- Turn off grid axis with numbering, and add lines 
  axlims = np.array([-1,1])
  ax.plot(axlims,0*axlims,0*axlims,'k-')	# x axis
  ax.plot(0*axlims,axlims,0*axlims,'k-')	# y axis
  ax.plot(0*axlims,0*axlims,axlims,'k-')	# z axis
  ax.set_axis_off()

  return M



# epg_show()
# Uses subplots to graphically show EPG decomposition of magnetisation.
# Subplots are rows for F+, F- and Z coefficients.  In most epg_ functions
# a matrix of the same size, "FZ" or "FpFmZ" is used to store the 
# coefficients.
#
# INPUTS:
#	FZ = EPG coefficient matrix to display
#	Nspins = Number of spins to use in graphical displays
#	frac = Additional dephasing.  Essentially adds to n for F_n 
#              states to show the dephasing/rephasing during gradient.
#       skipfull = True will not show "full" spins state in row 2, col 1
#	twists   = True to show as "twists" and "cosines" for F/Z vs
#		   False to have all vectors start at (0,0,0).
# OUTPUT:
#	none (plot is updated)
#

def epg_show(FZ,Nspins=19,frac=0,skipfull=False,twists=True):
#    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import    
#    fig = plt.figure()
#    ax = fig.gca(projection='3d')   

# STARTING to write this!
# Basic version works... lots to do!

  m = np.shape(FZ)[0]
  n = np.shape(FZ)[1]
  slabel = ('F$_{','F$_{-','Z$_{')

  #fig = plt.figure(figsize=plt.figaspect(np.float(m)/np.float(n)))
  fig = plt.figure(figsize=(3*n,3*m))	# Note (width,height)

  # -- Go through FZ Matrix
  for mm in range(m):
    for nn in range(n):
      figax = fig.add_subplot(m,n,nn+mm*n+1, projection='3d') 
      voxelvar=0	# -- Initialize - all M vectors start at (0,0,0)
      if (twists==True):
        voxelvar=1	# -- Twists for F by default
        if (mm > 1):
          voxelvar=2    # -- Z states with variation along x

      
      if (nn==0 and mm==1):	# -- Plot of all states/spins
        if (skipfull==False):
          epg_showstate(figax,FZ,frac=frac,Nspins=Nspins,voxvar=0) # All spins
          figax.title.set_text('All Spins')	# All spins combined
        else:
          figax.set_axis_off()

      else:			# -- Plot single state.
        Q = 0*FZ 
        Q[mm,nn]=FZ[mm,nn]      # -- Just 1 non-zero basis at a time
        epg_showstate(figax,Q,frac=frac,Nspins=Nspins,voxvar=voxelvar) 

        # -- Label state
        stateval = "%s%d}$ = %4.2f +%4.2fj" % (slabel[mm], \
		nn,np.real(FZ[mm,nn]),np.imag(FZ[mm,nn]))
        figax.title.set_text(stateval)		# F+ states

      # !!! Orient plots!  Could make F states from above so you just
      #     see the circle.  Once they are in color, this could be a good
      #     option.

  return
 
# Generate an image with a white diamond
def diamond(width=100, nx=256,ny=256,cx=0,cy=0):
    mx, my = np.meshgrid(np.arange(-nx/2,nx/2),np.arange(-ny/2,ny/2))
    mx = mx - cx
    my = my - cy
    mp = mx+my
    mm = mx-my
    mout = np.zeros((ny,nx)) 
    mout[mp<width/2] = mout[mp<width/2] + 1
    mout[mp>-width/2] = mout[mp>-width/2] + 1
    mout[mm<width/2] = mout[mm<width/2] + 1
    mout[mm>-width/2] = mout[mm>-width/2] + 1
    mout = np.maximum(mout,3)-3
    return mout

# Generate an image with a white circle
def circ(radius=100, nx=256,ny=256,cx=0,cy=0):
    mx, my = np.meshgrid(np.arange(-nx/2,nx/2),np.arange(-ny/2,ny/2))
    mx = mx - cx
    my = my - cy
    mout = mx**2+my**2
    r2 = radius * radius
    mout[mout<r2] = -1000
    mout[mout>r2] = 1000
    mout = (-mout/2000)+1
    return mout

# Generate an image with points increasing with angle
def angleim(nx=256,ny=256):
    mx, my = np.meshgrid(np.arange(-nx/2,nx/2),np.arange(-ny/2,ny/2))
    mc = mx + 1j*my
    m = np.angle(mc) + np.pi
    return(m)


def makenoiseykspace(nscale = 0.001):
    im = circ(100, 256,256)
    ksp = np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(im)))
    kspace = []
    for n in np.arange(1000):
        kspace.append(ksp+1j*nscale*np.random.randn((256,256))+nscale*np.random.randn((256,256)))
    return np.array(kspace)

def mingrad(area,Gmax = 50,Smax = 200,dt = 0.004, gamma = 42.58):
    area = area*100
    Atri = gamma *Gmax**2/Smax
    
    if area <= Atri:
        tramp = np.sqrt(area/gamma/Smax)
        Nramp = np.ceil(tramp/dt)
        g = np.arange(Nramp)*Smax*dt
        g = np.array((g,np.fliplr(g)))
    else:
        tramp = Gmax/Smax
        Nramp = np.ceil(tramp/dt)
        gramp = np.arange(Nramp)/Nramp*Gmax
        Nplat = np.ceil(area/gamma/Gmax/dt - Gmax/Smax/dt)
        g = np.array([gramp, Gmax*np.ones(Nplat), np.fliplr(gramp)])
    t = np.arange(len(g))*dt
    g = g * area/gamma/np.sum(g)/dt
    return g, t



        
def nlegend():
    pass

def plotc():
    pass

def plotgradinfo():
    pass



def setprops():
    import matplotlib.style as style
    style.use('seaborn') # 'ggplot'


def sinc(x):
    return np.sinc(x)

def sweptfreqrf():
    pass

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


def whirl(N,res,fov, tsample = 0.000004,
          upsamp = 16, gmax = 3.9, smax = 14500, gamma = 4258):
     
    
    dT = tsample/upsamp
    kmax = 0.5/res
    delta = N/2./np.pi/fov
    r1 = 1./np.sqrt(5)*delta
    r2 = 3./np.sqrt(5)*delta
    
    
    Gc= np.sqrt(2*smax*r1/gamma)
    g = np.arange(0,Gc, smax*dT)
    k = np.cumsum(g)*gamma*dT
    ng = len(g)
    ng1 = ng*1  # multiply by one to create a new instance because python is funny
    
    G = g[-1]
    r = k[-1]
    kk = k[-1]
    phi = 0
    
    # pre allocate space?
    maxng = 10000*upsamp
    if maxng > len(g):
        g = np.pad(g, maxng-len(g), 'constant').astype(np.complex)
        k = np.pad(k, maxng-len(k), 'constant').astype(np.complex)
    Grec = np.array(g)
    Grec[10000] = 0
    phirec = 0*Grec
    
    done = False
    dphi = smax/np.abs(G)*dT
    
    while r<r2:
        ng = ng+1
        phi = phi+dphi
        g[ng] = G*np.exp(1j*phi)
        kk = kk+g[ng]*dT*gamma
        k[ng] = kk
        r = np.abs(kk)
    
    ng2 = ng*1  # multiply by one to create a new instance because python is funny
    
    dG = 0
    
    
    
    while r<kmax and ng<maxng:
        Gapp = G
        r = np.abs(k[ng]+gamma*g[ng]*dT/2)
        ng = ng+1
        
        dphi = gamma*Gapp/np.sqrt(r**2-delta**2)*dT
        dG = np.sqrt(smax - (Gapp*dphi/dT)**2)*dT
        
        G = G+dG
        if G>= gmax:
            G = gmax
            dG = 0
            
        phi = phi+dphi
        g[ng] = G*np.exp(1j*phi)
        
        kk = kk+g[ng]*dT*gamma
        k[ng] = kk
        r = np.abs(kk)
        
        if ng>= maxng:
            print('max points', maxng)
            
            
    k = k[:ng:upsamp]
    g = g[:ng:upsamp]
    
    g1 = g[:np.round(ng1/upsamp)]
    g2 = g[np.round(ng1/upsamp):np.round(ng2/upsamp)]
    g3 = g[np.round(ng2/upsamp)+1:len(g)]
    
    return g, g1, g2, g3
        
        

def vds(smax, gmax, T, N, Fcoeff, rmax, z=0,
        gamma = 4258, oversamp = 8):
    
    To = T*1./oversamp
    
    q0 = 0.
    q1 = 0.
    
    Nprepare = 10000
    theta = np.zeros(Nprepare)
    r = np.zeros(Nprepare)
    
    r0 = 0
    r1 = 1
    
    time = np.zeros(Nprepare)
    t = 0
    count = 1
    
    theta = np.zeros(Nprepare)
    r = np.zeros(Nprepare)
    time = np.zeros(Nprepare)
        
    while r < rmax:
        print('Error findq2r2 does not exist')
        
        q1 = q1+q2*To
        q0 = q0 + q1*To
        t = t+To
        r1 = r1+r2*To
        r0 = r0 + r1*To
        
        count = count+1
        theta[count] = q0
        r[count] = r0
        time[count] = t
        
        if np.remainder(count,100) == 0:
            print('points')
    r = r[oversamp/2:count:oversamp]
    theta = theta[oversamp/2:count:oversamp]
    time = time[oversamp/2:count:oversamp]
    
    ltheta = 4*np.floor(len(theta)/4)
    r = r[:ltheta]
    theta = theta[:ltheta]
    time = time[:ltheta]
    
    r = r*np.exp(1j*theta)
    g = 1./gamma*(np.gradient(g))/T
    
    s = np.gradient(g)
    

def qdf(a,b,c):
    return np.roots([a,b,c])


def findq2r2(smax, gmax, r, r1, T, Ts, N, Fcoeff, rmax, z):
    gamma = 4258.
    smax = smax + z*gmax
    F= 0.
    dFdr = 0
    for rind in np.arange(Fcoeff):
        F = F+Fcoeff[rind]*(r/rmax)**(rind-1)
        if rind>1:
            dFdr = dFdr + (rind-1)*Fcoeff[rind]*(r/rmax)**(rind-2)/rmax
            
    GmaxFOV = 1./gamma/F/Ts
    Gmax = np.min([GmaxFOV, gmax])
    maxr1 = np.sqrt((gamma*Gmax)**2)/(1+(2*p.pi*F*r/N)**2)

    if r1>maxr1:
        r2 = (maxr1-r1)/T
        
        twopiFoN = 2*np.pi*F/N
        twopiFoN2 = twopiFoN**2
        
        A = 1+ twopiFoN2*r*r
        
        B = 2*twopiFoN2*r*r1*r1 + 2*twopiFoN2/F*dFdr*r*r*r1*r1 + 2*z*r1 + 2*twopiFoN2*r1*r
        C1 = twopiFoN2**2*r*r*r1**4 + 4*twopiFoN2*r1**4 + (2*pi/N*dFdr)**2*r*r*r1**4 + 4*twopiFoN2/F*dFdr*r*r1**4 - (gamma)**2*smax**2
        C2 = z*(z*r1**2 + z*twopiFoN2*r1**2 + 2*twopiFoN2*r1**3*r + 2*twopiFoN2/F*dFdr*r1**3*r)
        C = C1+C2;

        theroots = np.roots((A,B,C))
        r2 = np.real(theroots[0])
        
        slew = 1./gamma*(r2*twopiFoN2*r*r1**2+1j*twopiFoN*(2*r1**2+r*r2+dFdr/F*r*r1**2))
        sr = np.abs(slew)/smax
        
        if np.abs(slew)/smax > 1.01:
            print('slew',slew)
            
    q2 = 2*np.pi/N*dFdr*r1**2+2*np.pi/F/N*r2
    
    return q2, r2


def vecdcf(g,k = None,N = None,res = None,FOV = None):
    if np.shape(g)[1]==2:
        g = g[:,0]+1j*g[:,1]
        
    if k is None:
        k = np.minimum.accumulate(g, 0.5)  # this isn't equivalent.  hmm
        k = k * 0.5/ np.max(np.abs(k))  # watch dimension collapse      
        
    Nsamps = len(k)
    
    g3 = np.array([np.real(g), np.imag(g), 0*np.real(g)])
    k3 = np.array([np.real(k), np.imag(k), 0*np.real(k)])
    
    dcf = 0*k
    
    for m in np.arange(Nsamps):
        dcf[m] = np.linalg.norm()
    # I'm afraid I have to stop for the moment.  Pull!




    

        
        
        
        
        






    







