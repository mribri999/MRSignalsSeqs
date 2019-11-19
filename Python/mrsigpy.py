# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 08:57:48 2019
@author: Joshua Kaggie, Brian Hargreaves
"""


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



import numpy as np
from matplotlib.pyplot import figure, plot, subplot, title, xlabel, ylabel


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
    
    
def dispim(im,low=0,high = None):
    im = np.squeeze(im)
    
    if high is None:
        immax = np.max(np.abs(im))
        imstd = np.std(np.abs(im))
        high = immax - 0.5 * imstd
        
    scale = 256/(high-low)
    offset = scale*low
    
    from matplotlib import cm
    colormap = cm.get_cmap()
    if colormap.is_gray():        
        print('set the colormap')
    plt.imshow(np.abs(im))
    plt.axis('square')
        
    
def displogim(im):
    im = np.squeeze(im)
    lowexp = 0
    im = np.log(np.abs(im))
    
    f = np.where(im<lowexp)
    im[f] = lowexp
    
    im = im - lowexp
    dispim(im)
    
def ft(dat):
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(dat)))

def gaussian(x,mn,sig):
    return np.exp(-(x-mn)**2/(2*sig**2))/np.sqrt(2*np.pi)/sig

def ghist(data,gmean = None,gsig = None,bins = None,
          gtitle = 'Data and Gaussian',gleg1,gleg2):
    N  = len(data)

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

def lplot(xlab = None,ylab = None,tit = None,ax = None,grid = None):
    if xlab is not None:
        xlabel(xlab)
    if ylab is not None:
        ylabel(ylab)
    if title is not None:
        title(tit)
    if ax is not None:
        plt.axis(ax)
                
                                

def lsfatwater():
    pass


def magphase(x,arr):
    mag = np.abs(arr)
    phase = np.angle(arr)
    
    fig1 = plt.subplot(2,1,1)
    plt.plot(x,phase/np.pi)
    


def circ(radius, nx,ny):
    mx, my = np.meshgrid(np.arange(-nx/2,nx/2+1),np.arange(-ny/2,ny/2+1))
    mout = mx**2+my**2
    mout[mout<radius] = -1000
    mout[mout>radius] = 1000
    mout = (-mout/2000)+1
    return mout

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



def setprops():
    pass

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


def time2freq(t):
    dt = t[1]-t[0]
    num_p = len(t)
    df = 1./ num_p/dt
    f = np.arange(np.ceil((num_p-1.)/2.), np.floor((num_p-1.)/2.))*df
    return f

def vds():
    pass

def vecdcf():
    pass

def whirl():
    pass



#left handed
def xrot(angle = 0., in_degs = True):
    if in_degs:
        angle = angle*np.pi/180.
    c = np.cos(angle)    
    s = np.sin(angle)    
    M = np.array([[1.,0.,0.],[0., c, s],[0,-s, c]])
    return M

#This is also equivalent to help
def yrot(angle = 0., in_degs = True):
    if in_degs:
        angle = angle*np.pi/180.
    c = np.cos(angle)    
    s = np.sin(angle)    
    M = np.array([[c,0.,-s],[0.,1.,0.],[s,0.,c]])
    return M



def zrot(angle = 0., in_degs = True):
    'This is equivalent to help'
    if in_degs:
        angle = angle*np.pi/180.
    c = np.cos(angle)    
    s = np.sin(angle)    
    M = np.array([[c,s,0.],[-s,c,0],[0.,0.,1.]])
    return M

    

        
        
        
        
        






    








