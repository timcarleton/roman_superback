import numpy as np
from scipy.interpolate import SmoothBivariateSpline,LinearNDInterpolator
from . import combinespikes
from scipy import linalg
import matplotlib.pyplot as plt
#import plotspike
from matplotlib import patches
import os

#complication - x1, y1, x2, y2 do not move through focal plane in nice way - this is because of parallax, along with the intersection of the bar with the circular aperture
#solution - find x1real,y1real,x2real,y2real,zreal
#based on x1real, y1real,x2real,y2real,zreal find x1app,y1app,x2app,y2app from parallax
#based on x1app,y1app,x2app,y2app, find intersection with outer circle - x2int, y2int
#based x2int, y2int become new x2,y2. based on x2, y2,thetaapp,and length, get new x1, y1
#this doesn't quite work because, for some unknown reason, the theta values of the lolo diffusers change slightly across the field of view -shouldn't happen in parallax case?
#intermediate hack - for unmeasured lolo spikes, replace with center of mirror strut


def fitplane(x,y,z):
    A = np.c_[x, y, np.ones(len(x))] #ones are weights
    C,_,_,_ = linalg.lstsq(A, z)    # coefficients
    return C

def evalplane(coef,x,y):
    return coef[0]*x+coef[1]*y+coef[2]

def getx1x2(phi2,theta,l,r2=1.2):

    x2=r2*np.cos(phi2)
    y2=r2*np.sin(phi2)

    x1=x2-l*np.cos(theta)
    y1=y2-l*np.sin(theta)
    return x1,y1,x2,y2


refphi={}
reftheta={}
refw={}
refphi2={}
reftheta2={}
refw2={}
reflength={}
ref2x1={}
ref2x2={}
ref2y1={}
ref2y2={}

widths = {'F062':0.02639,'F087':0.02635,'F106':0.02629,'F129':0.02634,'F158':0.02634,'F184':0.2768,'F213':0.03695,'F146':0.02636}
refpos = np.loadtxt(os.environ['SUPERBACK_PATH']+'/data/psf/thetaxy_ref.txt')
radius=2363.12/2000

rotposx = refpos[:,0]*np.cos(-np.pi/3) - refpos[:,1]*np.sin(-np.pi/3)
rotposy = refpos[:,0]*np.sin(-np.pi/3) + refpos[:,1]*np.cos(-np.pi/3)

for filt in ['F062','F106']:
    
    refphi[filt]={}
    reftheta[filt]={}
    refw[filt]={}
    refphi2[filt]={}
    reftheta2[filt]={}
    refw2[filt]={}
    reflength[filt]={}

    ref2x1[filt]={}
    ref2x2[filt]={}
    ref2y1[filt]={}
    ref2y2[filt]={}
            
    for spike in range(1,7,1):
        dat=np.loadtxt(os.environ['SUPERBACK_PATH']+'/data/psf/bestfitphitheta_'+filt+'_'+str(spike)+'.txt')
        lengths=np.loadtxt(os.environ['SUPERBACK_PATH']+'/data/psf/fitlength_'+str(spike)+'.txt')
        refphi[filt][str(spike)]=dat[:,0]
        reftheta[filt][str(spike)]=dat[:,1]
        refw[filt][str(spike)]=dat[:,2]

        refphi2[filt][str(spike)]=dat[:,3]
        reftheta2[filt][str(spike)]=dat[:,4]
        #refw2[filt][str(spike)]=dat[:,5]
        refw2[filt][str(spike)]=np.zeros(len(dat))+widths[filt]
        reflength[filt][str(spike)]=lengths
        x1ref,y1ref,x2ref,y2ref=getx1x2(dat[:,3],dat[:,4],lengths,radius)
        ref2x1[filt][str(spike)]=x1ref
        ref2x2[filt][str(spike)]=x2ref
        ref2y1[filt][str(spike)]=y1ref
        ref2y2[filt][str(spike)]=y2ref

splphi={}
spltheta={}
splw={}
splphi2={}
spltheta2={}
splw2={}
spllength={}
coef2x1={}
coef2x2={}
coef2y1={}
coef2y2={}



interp2=True
for filt in ['F062','F106']:
    
    splphi[filt]={}
    spltheta[filt]={}
    splw[filt]={}
    splphi2[filt]={}
    spltheta2[filt]={}
    splw2[filt]={}
    spllength[filt]={}
    coef2x1[filt]={}
    coef2x2[filt]={}
    coef2y1[filt]={}
    coef2y2[filt]={}
            
    for spike in range(1,7,1):
        splphi[filt][str(spike)]=SmoothBivariateSpline(refpos[:,0],refpos[:,1],refphi[filt][str(spike)],kx=2,ky=2)
        spltheta[filt][str(spike)]=SmoothBivariateSpline(refpos[:,0],refpos[:,1],reftheta[filt][str(spike)],kx=2,ky=2)
        splw[filt][str(spike)]=SmoothBivariateSpline(refpos[:,0],refpos[:,1],refw[filt][str(spike)],kx=2,ky=2)
        spllength[filt][str(spike)]=SmoothBivariateSpline(refpos[:,0],refpos[:,1],reflength[filt][str(spike)],kx=1,ky=1)

        #because there is relatively little data for spike2 values, need to convert phi and theta into x1, x2, y1, y2. Theta2 does have some noticeable variation, so can't use constant value
        wwith2=np.where(np.isfinite(refphi2[filt][str(spike)]))[0]
        if interp2:
            if len(wwith2)>=9:
                splphi2[filt][str(spike)]=SmoothBivariateSpline(rotposx[wwith2],rotposy[wwith2],refphi2[filt][str(spike)][wwith2],kx=2,ky=2)
                spltheta2[filt][str(spike)]=SmoothBivariateSpline(rotposx[wwith2],rotposy[wwith2],reftheta2[filt][str(spike)][wwith2],kx=2,ky=2)
                splw2[filt][str(spike)]=SmoothBivariateSpline(rotposx[wwith2],rotposy[wwith2],refw2[filt][str(spike)][wwith2],kx=2,ky=2)
            else:
                #need to fix this - linear interpolator breaks for positions outside reference points, so currently filling with mean.
                #should do extrapolation
                #splphi2[filt][str(spike)]=LinearNDInterpolator(refpos[wwith2,:],refphi2[filt][str(spike)][wwith2],fill_value=np.nanmean(refphi2[filt][str(spike)][wwith2]))
                #spltheta2[filt][str(spike)]=LinearNDInterpolator(refpos[wwith2,:],reftheta2[filt][str(spike)][wwith2],fill_value=np.nanmean(reftheta2[filt][str(spike)][wwith2]))
                #splw2[filt][str(spike)]=LinearNDInterpolator(refpos[wwith2,:],refw2[filt][str(spike)][wwith2],fill_value=np.nanmean(refw2[filt][str(spike)][wwith2]))
                tosplphi = np.array(refphi2[filt][str(spike)])
                tospltheta = np.array(reftheta2[filt][str(spike)])
                tosplw2 = np.array(refw2[filt][str(spike)])

                wnan=np.where(~np.isfinite(refphi2[filt][str(spike)]))[0]
                tosplphi[wnan]=refphi[filt][str(spike)][0]
                tospltheta[wnan]=reftheta[filt][str(spike)][0]
                tosplw2[wnan]=np.nanmean(refw2[filt][str(spike)])

                splphi2[filt][str(spike)]=SmoothBivariateSpline(rotposx,rotposy,tosplphi,kx=2,ky=2)
                spltheta2[filt][str(spike)]=SmoothBivariateSpline(rotposx,rotposy,tospltheta,kx=2,ky=2)
                splw2[filt][str(spike)]=SmoothBivariateSpline(rotposx,rotposy,tosplw2,kx=2,ky=2)

                
        else:
            x1,y1,x2,y2=getx1x2(refphi2[filt][str(spike)][wwith2],reftheta2[filt][str(spike)][wwith2],reflength[filt][str(spike)][wwith2],r2=radius)
            x10,y10,x20,y20=getx1x2(refphi[filt][str(spike)][0],reftheta[filt][str(spike)][0],reflength[filt][str(spike)][0],r2=radius)
            x1=np.append(x1,x10)
            y1=np.append(y1,y10)
            x2=np.append(x2,x20)
            y2=np.append(y2,y20)
            refx = np.append(refpos[wwith2,0],refpos[0,0])
            refy = np.append(refpos[wwith2,1],refpos[0,1])
            
            x1fit = fitplane(refx,refy,x1)
            x2fit = fitplane(refx,refy,x2)
            y1fit = fitplane(refx,refy,y1)
            y2fit = fitplane(refx,refy,y2)


            #if len(wwith2)>=4:
            splw2[filt][str(spike)]=LinearNDInterpolator(refpos[wwith2,:],refw2[filt][str(spike)][wwith2],fill_value=np.nanmean(refw2[filt][str(spike)][wwith2]))
            #else:
                

            coef2x1[filt][str(spike)] = x1fit
            coef2x2[filt][str(spike)] = x2fit
            coef2y1[filt][str(spike)] = y1fit
            coef2y2[filt][str(spike)] = y2fit            
        


        

def getromanspikes(filt,thetax,thetay,plot=False):

    phi=[]
    theta=[]
    width=[]
    phi2=[]
    theta2=[]
    width2=[]
    length=[]
    x0=[]
    y0=[]
    x02=[]
    y02=[]

    for spike in range(1,7,1):
        phi.append(splphi[filt][str(spike)](thetax,thetay)[0,0])
        theta.append(spltheta[filt][str(spike)](thetax,thetay)[0,0])
        width.append(splw[filt][str(spike)](thetax,thetay)[0,0])
        length.append(spllength[filt][str(spike)](thetax,thetay)[0,0])

        x1s,y1s,x2s,y2s=getx1x2(phi[-1],theta[-1],length[-1],r2=radius)
        x0.append(np.mean([x1s,x2s]))
        y0.append(np.mean([y1s,y2s]))

        thetaxrot = thetax*np.cos(-np.pi/3) - thetay*np.sin(-np.pi/3)
        thetayrot = thetax*np.sin(-np.pi/3) + thetay*np.cos(-np.pi/3)
        
        width2s=splw2[filt][str(spike)](thetax,thetay)

        if len(width2s.shape)!=0:
            width2.append(width2s[0,0])
        else:
            width2.append(width2s)

        if interp2:

            
            phi2s=splphi2[filt][str(spike)](thetaxrot,thetayrot)
            theta2s=spltheta2[filt][str(spike)](thetaxrot,thetayrot)
        
            if len(phi2s.shape)!=0:
                phi2.append(phi2s[0,0])
            else:
                phi2.append(phi2s)
                
            if len(theta2s.shape)!=0:
                theta2.append(theta2s[0,0])
            else:
                theta2.append(theta2s)

            x1s2,y1s2,x2s2,y2s2 = getx1x2(phi2[-1],theta2[-1],length[-1],r2=radius)
            x02.append(np.mean([x1s2,x2s2]))
            y02.append(np.mean([y1s2,y2s2]))
        else:
            x1s2=evalplane(coef2x1[filt][str(spike)],thetax,thetay)
            x2s2=evalplane(coef2x2[filt][str(spike)],thetax,thetay)
            y1s2=evalplane(coef2y1[filt][str(spike)],thetax,thetay)
            y2s2=evalplane(coef2y2[filt][str(spike)],thetax,thetay)

            x02.append(np.mean([x1s2,x2s2]))
            y02.append(np.mean([y1s2,y2s2]))
            theta2.append(np.atan2(y2s2-y1s2,x2s2-x1s2))

    allspikes=[]        
    for spike in range(6):
        spike1 = {}
        spike1['type']='bar'
        spike1['x0']=x0[spike]
        spike1['y0']=y0[spike]
        spike1['l']=length[spike]
        spike1['w']=width[spike]
        spike1['theta']=theta[spike]+np.pi/2
        spike1['amp']=1

        spike2 = {}
        spike2['type']='bar'
        spike2['x0']=x02[spike]
        spike2['y0']=y02[spike]
        spike2['l']=length[spike]
        spike2['w']=width2[spike]
        spike2['theta']=theta2[spike]+np.pi/2
        spike1['amp']=1

        
        combsp = combinespikes.combinespikes(spike1,spike2)

        if type(combsp)==tuple:

            for sp in combsp:
                sp['amp']=1
            allspikes +=list(combsp)
        else:
            combsp['amp']=1
            allspikes+=[combsp]

    return allspikes

