import numpy as np
import os
from scipy.special import j1
from . import getromanspikes
'''
This module calculates the surface brightness of the Roman WFI PSF at a given angular offset from 
a given source. This code is optimized for large angular offsets from the source, and takes a different
approach than STPSF. Rather than numerically calculating the fourier transform of the pupil to create
a template, this code uses analytic expressions for the fourier transform of the aperture and mirror struts.
As such, the surface brightness can be calculated for any arbitrary angular offset efficiently, without having
to oversample the pupil image. The cost of this efficiency is that, due to not modeling the details of
the pupil, it is not as accurate on small scales. It is able to reproduce the diffraction spikes produced
by the mirror struts.

Calling the calculator is simple:

>p0=romanpsf.romanpsf(x_offset_arcsec,y_offset_arcsec,wave,bars=bardict)

where p0 is the normalized surface brightness of the psf at x and y offsets from the source.
The input x and y offsets are assumed to be in the ideal coordinate system (https://roman-docs.stsci.edu/data-handbook/wfi-data-levels-and-products/coordinate-systems) in arcsec
wave is the wavelength for which the PSF is calculated in meters
bardict is a list of dictionaries describing the position-dependent parameters of the mirror struts as seen by the focal plane
additional useage to include the LOLO-diffuser struts is in development

defaultbarparam is a dictionary defined in defaultbars.py that sets the default parameters of this dictionary (calculated
for a source at the center of SCA01).


In general, you want to recalculate these parameters for a given source position to account for the parallax
of the mirror struts at different points in the focal plane. The getbarparamtheta method allows for the
straightforward calculation of this. Useage is as follows

getbarparamtheta(thetax,thetay,bardict)

thetax and thetay describe the position of the source in observatory (v2, v3) coordinates (in arcsec)
the method will fill in the necessary parameters to bardict given the offset, assuming the true strut
positions (this method modifies bardict, so store a copy if necessary)

TODO:
account for elliptical distortion of primary, secondary
include diffraction spikes from LOLO-diffuser struts. They should be aligned with the mirror struts in a reference position, but are visible for an offset source
'''

i=complex(0,1)
pi=np.pi
arcsectorad=pi/180/60/60

#set the default bar parameters given thetax,thetay = (1312.94914525,-1040.78537268) (i.e. the center of chip 1) and filter = 'F146'


defaultbarparam = getromanspikes.getromanspikes('F062',1312.94914525,-1040.78537268)
    

def barpsf(x,y,k,w,l,x0,y0,rotation=0):
    #calculates the fourier transform pattern of a bar with a given center and rotation

    #Rotate x and y to a coordinate system aligned with the bar
    xrot=x*np.cos(rotation)-y*np.sin(rotation)
    yrot=x*np.sin(rotation)+y*np.cos(rotation)

    #the resulting Fourier transform is a sinc in x, and y, with phases shifted to the central x and y
    return np.sinc(k*w*xrot/2/np.pi)*np.sinc(k*l*yrot/2/np.pi)*np.exp(-i*k*x0*x)*np.exp(-i*k*y*y0)*l*w

def trappsf(x,y,k,a,b,c,x0,y0,rotation=0):
    #calculates the fourier transform pattern of a right trapezoid with base c and heights a and b
    
    #Rotate x and y to a coordinate system aligned with the base
    xrot=x*np.cos(rotation)-y*np.sin(rotation)
    yrot=x*np.sin(rotation)+y*np.cos(rotation)

    #from mathematica
    return np.sinc(c*k*xrot/2/np.pi)*c*np.sinc(k*(a*(c-xrot)+b*xrot)*yrot/2/c/np.pi)*(a*(c-xrot)+b*xrot)/c*np.exp(-i*k*(c**2*xrot+a*c*yrot+(b-a)*xrot*yrot/2/c))*np.exp(-i*k*x0*xrot)*np.exp(-i*k*yrot*y0)


def romanpsf(xoffset,yoffset,wave,barparam=defaultbarparam,x0=0,y0=0,flength=18.813,ainner=.783/2,aouter=2.553/2,pixsize=10E-6):
    #calculates the PSF of Roman's WFI

    #pixscale
    pixscale=pixsize/flength
    
    #convert to radians
    x=xoffset*arcsectorad/pixscale*pixsize
    y=yoffset*arcsectorad/pixscale*pixsize

    #theta_perp for disk fourier transform
    theta=np.sqrt(x**2+y**2)
    
    k=2*np.pi/wave/flength

    #fourier transform of disk is J1 bessel function divided by k*a*theta_perp (https://adriftjustoffthecoast.wordpress.com/2013/06/06/2d-fourier-transform-of-the-unit-disk/)
    mirrorpsfouter=np.array(j1(k*aouter*theta)/(k*aouter*theta),dtype=complex)*2*np.pi*aouter**2
    mirrorpsfinner=np.array(j1(k*ainner*theta)/(k*ainner*theta),dtype=complex)*2*np.pi*ainner**2

    #subtract inner mirror from the outer
    final=mirrorpsfouter-mirrorpsfinner

    #subtract the bars
    bararea=[]
    for j in range(len(barparam)):
        if barparam[j]['type']=='bar':
            final-=barparam[j]['amp']*barpsf(x,y,k,barparam[j]['w'],barparam[j]['l'],barparam[j]['x0'],barparam[j]['y0'],rotation=barparam[j]['theta'])
            #keep track of the area of the bars for normalization
            bararea.append(barparam[j]['w']*barparam[j]['l'])
        else:
            final-=barparam[j]['amp']*trappsf(x,y,k,barparam[j]['a'],barparam[j]['b'],barparam[j]['base'],barparam[j]['x0'],barparam[j]['y0'],rotation=barparam[j]['theta'])
            bararea.append(barparam[j]['base']*(barparam[j]['a']+barparam[j]['b'])/2)

    #alow for a phase shift to a given x0 and y0 in case this is passed through more optics
    #but don't do it by default
    if not ((x0==0) and (y0==0)):
        final*=np.exp(-i*k*x0*x)*np.exp(-i*k*y*y0)

    #normalize the final pattern by k/2/pi/sqrt(area)*pixelsize
    pixscale=pixsize/flength/arcsectorad
    final*=k/2/np.pi/np.sqrt(np.pi*(aouter**2-ainner**2)-np.sum(bararea))*pixsize/pixscale

    return final

