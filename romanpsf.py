import numpy as np
import os
from scipy.special import jv

i=complex(0,1)
pi=np.pi

def barpsf(x,y,k,w,l,x0,y0,rotation=0):
    xrot=x*np.cos(rotation)-y*np.sin(rotation)
    yrot=x*np.sin(rotation)+y*np.cos(rotation)
    x0rot=x0*np.cos(rotation)-y0*np.sin(rotation)
    y0rot=x0*np.sin(rotation)+y0*np.cos(rotation)

    #return np.sinc(k*w*xrot/2/np.pi)*np.sinc(k*l*yrot/2/np.pi)*np.exp(-i*k*x0*x)*np.exp(-i*k*y*y0)*k/2/np.pi*np.sqrt(l*w)
    return np.sinc(k*w*xrot/2/np.pi)*np.sinc(k*l*yrot/2/np.pi)*np.exp(-i*k*x0*x)*np.exp(-i*k*y*y0)*k/2/np.pi*l*w

def thickbarpsf(x,y,wave,w,l,x0,y0,z0,dzdy,rotation=0):
    xrot=x*np.cos(rotation)-y*np.sin(rotation)
    yrot=x*np.sin(rotation)+y*np.cos(rotation)
    x0rot=x0*np.cos(rotation)-y0*np.sin(rotation)
    y0rot=x0*np.sin(rotation)+y0*np.cos(rotation)

    xpart=np.exp((i*2*np.pi*xrot/wave*w)/(dzdy*w-2*z0))*(dzdy*(dzdy*(1+np.exp((-4*i*2*np.pi*xrot/wave*w*z0)/(dzdy**2*w**2-4*z0**2)))*w+2*(-1+np.exp((-4*i*2*np.pi*xrot/wave*w*z0)/(dzdy**2*w**2-4*z0**2)))*z0)-(2*i*2*np.pi*xrot/wave*z0*expi((2*i*2*np.pi*xrot/wave*z0)/(dzdy**2*w-2*dzdy*z0)))/np.exp((2*i*2*np.pi*xrot/wave*z0)/(dzdy**2*w-2*dzdy*z0))+(2*i*2*np.pi*xrot/wave*z0*expi((-2*i*2*np.pi*xrot/wave*z0)/(dzdy**2*w+2*dzdy*z0)))/np.exp((2*i*2*np.pi*xrot/wave*z0)/(dzdy**2*w-2*dzdy*z0)))/(2.*dzdy**2)
    ypart=np.sinc(2*np.pi*yrot*l/wave/z0)/2
    return xpart*ypart

def romanpsf(xp,yp,wave,a=0.7556,x0=0,y0=0,d=.008,flength=20,wbar=.085,baramp1=1,baramp2=1,baramp3=1,baramp4=1,baramp5=1,baramp6=1,bartheta1=256,bartheta2=342,bartheta3=318,bartheta4=223,bartheta5=199,bartheta6=103,ainner=.92,aouter=2.4,lbar=.75,oversample=1):

    x=xp*10E-6/oversample
    y=yp*10E-6/oversample

    theta=np.sqrt(x**2+y**2)/flength
    
    k=2*np.pi/wave/flength

    area=pi*(aouter**2-ainner**2)
    norm=k/2/np.pi/np.sqrt(area)

    mirrorpsfouter=2*jv(1,k*aouter/2*theta)/(k*aouter*theta)
    mirrorpsfinner=2*jv(1,k*ainner/2*theta)/(k*ainner*theta)

    return (mirrorpsfouter-mirrorpsfinner-
            baramp1*barpsf(x,y,k,wbar,lbar,0.416,0.728,rotation=bartheta1*pi/180)-
            baramp2*barpsf(x,y,k,wbar,lbar,0.864,-0.026,rotation=bartheta2*pi/180)-
            baramp3*barpsf(x,y,k,wbar,lbar,0.477,-0.678,rotation=bartheta3*pi/180)-
            baramp4*barpsf(x,y,k,wbar,lbar,-0.427,-0.665,rotation=bartheta4*pi/180)-
            baramp5*barpsf(x,y,k,wbar,lbar,-0.800,-0.013,rotation=bartheta5*pi/180)-
            baramp6*barpsf(x,y,k,wbar,lbar,-0.366,0.694,rotation=bartheta6*pi/180))
    
    

