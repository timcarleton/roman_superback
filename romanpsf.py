import numpy as np
import os
from scipy.special import j1
import getbarparam
import copy
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
bardict is a dictionary describing the position-dependent parameters of the mirror struts as seen by the focal plane
additional useage to include the LOLO-diffuser struts is in development

defaultbarparam is a dictionary defined below that sets the default parameters of this dictionary (calculated
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

defaulttruebars={}
defaulttruebars['1']={}
defaulttruebars['2']={}
defaulttruebars['3']={}
defaulttruebars['4']={}
defaulttruebars['5']={}
defaulttruebars['6']={}

defaulttruebars['1']['theta0']=1.26321156e+00
defaulttruebars['1']['thetax1']=-1.33813468e+00
defaulttruebars['1']['thetay1']=-1.78752296e+00
defaulttruebars['1']['thetax2']=-1.95988066e+02
defaulttruebars['1']['thetay2']=1.32420061e+02
defaulttruebars['1']['thetaxy']=-1.21648756e+02

defaulttruebars['2']['theta0']=2.90956877e+00
defaulttruebars['2']['thetax1']=1.42334965e+00
defaulttruebars['2']['thetay1']4.73671855e+00
defaulttruebars['2']['thetax2']=1.26269313e+02
defaulttruebars['2']['thetay2']=-1.26721123e+02
defaulttruebars['2']['thetaxy']=7.28451245e+01

defaulttruebars['3']['theta0']=2.31524411e+00
defaulttruebars['3']['thetax1']=-2.59609541e+00
defaulttruebars['3']['thetay1']=1.15389898e-01
defaulttruebars['3']['thetax2']=1.78182896e+02
defaulttruebars['3']['thetay2']=1.31457694e+02
defaulttruebars['3']['thetaxy']=1.44438833e+02

defaulttruebars['4']['theta0']=8.31332570e-01
defaulttruebars['4']['thetax1']=1.40658521e+00
defaulttruebars['4']['thetay1']-1.22379163e+00
defaulttruebars['4']['thetax2']=-3.29482919e+01
defaulttruebars['4']['thetay2']=-2.20938807e+02
defaulttruebars['4']['thetaxy']=-5.30852335e+01

defaulttruebars['5']['theta0']=2.20435522e-01
defaulttruebars['5']['thetax1']=2.70536433e+00
defaulttruebars['5']['thetay1']=3.89843943e+00
defaulttruebars['5']['thetax2']=1.62554514e+02
defaulttruebars['5']['thetay2']=-1.03490979e+02
defaulttruebars['5']['thetaxy']=1.01776819e+02

defaulttruebars['6']['theta0']=1.89001822e+00
defaulttruebars['6']['thetax1']=-1.00974915e-01
defaulttruebars['6']['thetay1']=-2.59070366e+00
defaulttruebars['6']['thetax2']=-1.59612554e+02
defaulttruebars['6']['thetay2']=1.33837863e+02
defaulttruebars['6']['thetaxy']=-9.56383102e+01

#center of chip1
thetax0=1312.9491452484797*arcsectorad
thetay0=-1040.7853726755036*arcsectorad

defaultbarparam={}
defaultbarparam['1']={}
defaultbarparam['2']={}
defaultbarparam['3']={}
defaultbarparam['4']={}
defaultbarparam['5']={}
defaultbarparam['6']={}

#effective length of optical path in m
defaultbarparam['zref']=18.813

#width of secondary mirror struts
defaultbarparam['1']['w']=7.43934760e-02
defaultbarparam['2']['w']=7.48706892e-02
defaultbarparam['3']['w']=7.46923847e-02
defaultbarparam['4']['w']=7.51136392e-02
defaultbarparam['5']['w']=7.52385870e-02
defaultbarparam['6']['w']=7.45845916e-02

#length of secondary mirror struts
defaultbarparam['1']['l']=1.81879241
defaultbarparam['2']['l']=1.78030693
defaultbarparam['3']['l']=1.81857799
defaultbarparam['4']['l']=1.81168471
defaultbarparam['5']['l']=1.78371425
defaultbarparam['6']['l']=1.83067248

#centroid of secondary mirror struts
#just results in a phase offset
defaultbarparam['1']['x0']=1.4335505999999998
defaultbarparam['1']['y0']=2.72694435

defaultbarparam['2']['x0']=3.08062875
defaultbarparam['2']['y0']=-0.05643900000000002

defaultbarparam['3']['x0']=1.6969326
defaultbarparam['3']['y0']=-2.49930705

defaultbarparam['4']['x0']=-1.70916105
defaultbarparam['4']['y0']=-2.4588590999999997

defaultbarparam['5']['x0']=-3.09379785
defaultbarparam['5']['y0']=-0.023516249999999933

defaultbarparam['6']['x0']=-1.6470781499999998
defaultbarparam['6']['y0']=2.74011345

for i in range(1,7):
    getbarparamtheta(thetax0,thetay0,defaultbarparam,truebars=defaulttruebars)

#normalized amplitude of diffraction spikes
#investigation involving JWST PSFs finds that the diffraction spikes
#have a slightly lower amplitude than one would expect for their geometry
#this may be due to wierdness in the multi-plane diffraction or some other
#characteristic of the struts. For now, scaling the diffraction spikes by
#a certain amount seems to work.
#This technically messes with the normalization but not by a noticable
#amount because the struts are a small fraction of the area

defaultbarparam['1']['amp']=1
defaultbarparam['2']['amp']=1
defaultbarparam['3']['amp']=1
defaultbarparam['4']['amp']=1
defaultbarparam['5']['amp']=1
defaultbarparam['6']['amp']=1
    

def barpsf(x,y,k,w,l,x0,y0,rotation=0):
    #calculates the fourier transform pattern of a bar with a given center and rotation

    #rotate x and y to a coordinate system aligned with the bar
    xrot=x*np.cos(rotation)-y*np.sin(rotation)
    yrot=x*np.sin(rotation)+y*np.cos(rotation)

    #the resulting Fourier transform is a sinc in x, and y, with phases shifted to the central x and y
    return np.sinc(k*w*xrot/2/np.pi)*np.sinc(k*l*yrot/2/np.pi)*np.exp(-i*k*x0*x)*np.exp(-i*k*y*y0)*l*w


def romanpsf(xoffset,yoffset,wave,barparam=defaultbarparam,x0=0,y0=0,flength=18.813,ainner=.783,aouter=2.553,pixsize=10E-6):
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
    for j in range(6):
        final-=barparam[str(j+1)]['amp']*barpsf(x,y,k,barparam[str(j+1)]['w'],barparam[str(j+1)]['l'],barparam[str(j+1)]['x0'],barparam[str(j+1)]['y0'],rotation=barparam[str(j+1)]['theta'])
        #keep track of the area of the bars for normalization
        bararea.append(barparam[str(j+1)]['w']*barparam[str(j+1)]['l'])

    #alow for a phase shift to a given x0 and y0 in case this is passed through more optics
    #but don't do it by default
    if not ((x0==0) and (y0==0)):
        final*=np.exp(-i*k*x0*x)*np.exp(-i*k*y*y0)

    #normalize the final pattern by k/2/pi/sqrt(area)*pixelsize
    final*=k/2/np.pi/np.sqrt(np.pi*(aouter**2-ainner**2)-np.sum(bararea))*pixsize
    
    return final

def getbarparamtheta(thetaxi,thetayi,barparam,truebars=defaulttruebars):
    #get the apparent projected angle of the mirror struts
    #given the angular offset from the telescope axis
    
    #do quadratic x,y polynomial for theta instead of full geometric solution
    #the geometry must be just a bit too complicated

    dx=thetax-thetax0
    dy=thetay-thetay0
    for j in range(6):
        bari=str(j+1)
        barparam[bari]['theta']=truebars[bari]['theta0']+truebars[bari]['thetax1']*dx+truebars[bari]['thetay1']*dy+truebars[bari]['thetax2']*dx**2+truebars[bari]['thetay2']*dy**2+truebars[bari]['thetaxy']*dx*dy

    
def getbardict(thetax,thetay,truebars=defaulttruebars):
    #same as above, but returns a new dict instead of modifying an existing one
    barparam=copy.deepcopy(defaultbarparam)
    
    dx=thetax-thetax0
    dy=thetay-thetay0

    for j in range(6):
        bari=str(j+1)
        barparam[bari]['theta']=truebars[bari]['theta0']+truebars[bari]['thetax1']*dx+truebars[bari]['thetay1']*dy+truebars[bari]['thetax2']*dx**2+truebars[bari]['thetay2']*dy**2+truebars[bari]['thetaxy']*dx*dy

    return barparam
