import numpy as np
import copy

def getrot(x,y,theta):
    return x*np.cos(theta)-y*np.sin(theta),x*np.sin(theta)+y*np.cos(theta)

def fixtrap(trap):
    fixed=copy.deepcopy(trap)
    if trap['base']<0 and trap['a']>0 and trap['b']>0:
        fixed['x0']=trap['x0']+trap['base']*np.cos(trap['theta'])
        fixed['y0']=trap['y0']+trap['base']*np.sin(trap['theta'])
        fixed['base']*=-1
        fixed['a']=1.0*trap['b']
        fixed['b']=1.0*trap['a']
        
    elif trap['base']>0 and trap['a']<0 and trap['b']<0:
        fixed['theta']=trap['theta']+np.pi
        fixed['x0']=trap['x0']+trap['base']*np.cos(trap['theta'])
        fixed['y0']=trap['y0']+trap['base']*np.sin(trap['theta'])
        fixed['a']=-1.0*trap['b']
        fixed['b']=-1.0*trap['a']

    elif trap['base']>0 and trap['a']>0 and trap['b']>0:
        #print('unexpected trapezoid',trap)
        None
    return fixed
        
    
    
def combinespikes(spike1,spike2):
    #given two bar spikes that may be overlapping, returns some number of spikes accounting for their overlap
    #if they never overlap, it returns a copy of the two spikes
    #if they completely overlap, it returns one spike that has the properties of the overlapping spike
    #if they partially overlap, it returns up to three spikes describing the final set

    #assumes that spike2 is thinner than spike1 - makes the math a lot easier

    rotangle=-(spike1['theta']-np.pi/2)
    flip=2*int(spike2['theta']>spike1['theta'])-1
    #if spike2['theta']>spike1['theta'], all methods work
    #to get them to work if spike2['theta']<spike1['theta'], use this flip bit
    
    xout1=spike1['l']/2*np.cos(spike1['theta']-np.pi/2)+spike1['x0']
    yout1=spike1['l']/2*np.sin(spike1['theta']-np.pi/2)+spike1['y0']

    xout2=spike2['l']/2*np.cos(spike2['theta']-np.pi/2)+spike2['x0']
    yout2=spike2['l']/2*np.sin(spike2['theta']-np.pi/2)+spike2['y0']

    xin1=-spike1['l']/2*np.cos(spike1['theta']-np.pi/2)+spike1['x0']
    yin1=-spike1['l']/2*np.sin(spike1['theta']-np.pi/2)+spike1['y0']

    xin2=-spike2['l']/2*np.cos(spike2['theta']-np.pi/2)+spike2['x0']
    yin2=-spike2['l']/2*np.sin(spike2['theta']-np.pi/2)+spike2['y0']

    shiftx=1.0*xin1
    shifty=1.0*yin1

    #print(np.atan2(yout1-shifty,xout1-shiftx))
    #print(xin1-shiftx,yin1-shifty,xout1-shiftx,yout1-shifty)
    #print(xin2-shiftx,yin2-shifty,xout2-shiftx,yout2-shifty)
    #print('r',rotangle)
    xi1r,yi1r=getrot(xin1-shiftx,yin1-shifty,rotangle)
    xi2r,yi2r=getrot(xin2-shiftx,yin2-shifty,rotangle)

    xo1r,yo1r=getrot(xout1-shiftx,yout1-shifty,rotangle)
    xo2r,yo2r=getrot(xout2-shiftx,yout2-shifty,rotangle)


    yi2r*=flip
    yo2r*=flip
    
    #print('a',xi1r,yi1r,xo1r,yo1r)
    #print(xi2r,yi2r,xo2r,yo2r)
    
    m2=(yo2r-yi2r)/(xo2r-xi2r)


    appw2=spike2['w']/np.cos(spike2['theta']-spike1['theta'])

    ytopout0 = yo2r + appw2/2
    ytopin0 = yi2r + appw2/2
    ybotout0 = yo2r - appw2/2
    ybotin0 = yi2r - appw2/2


    xintersecttoptop=(spike1['w']/2-ytopin0+m2*xi2r)/m2 #top of spike 2 and top of spike 1
    xintersecttopbottom=(spike1['w']/2-ybotin0+m2*xi2r)/m2 #bottom of spike 2 and top of spike 1
    xintersectbottomtop=(-spike1['w']/2-ytopin0+m2*xi2r)/m2 #top of spike 2 and bottom of spike 1
    xintersectbottombottom=(-spike1['w']/2-ybotin0+m2*xi2r)/m2 #bottom of spike 2 and bottom of spike 1


    #print(xintersectbottomtop)
    xtopout = xo2r - spike2['w']/2*np.sin(spike2['theta']-spike1['theta'])
    xbotout = xo2r + spike2['w']/2*np.sin(spike2['theta']-spike1['theta'])

    xtopin = xi2r - spike2['w']/2*np.sin(spike2['theta']-spike1['theta'])
    xbotin = xi2r + spike2['w']/2*np.sin(spike2['theta']-spike1['theta'])

    ytopout = m2 * (xo1r - xtopout) + ytopout0
    ytopin = m2 * (xi1r - xtopin) + ytopin0
    ybotout = m2 * (xo1r - xbotout) + ybotout0
    ybotin = m2 * (xi1r - xbotin) + ybotin0
    

    
    #ytoptop=m2*(xintersecttoptop-xi2r)+yi2r
    #ytopbottom=m2*(xintersecttopbottom-xi2r)+yi2r
    #ybottomtop=m2*(xintersectbottomtop-xi2r)+yi2r
    #ybottombottom=m2*(xintersectbottombottom-xi2r)+yi2r
        
    #print('y',yo2r,yi2r,xo2r,xi2r,m2,appw2,appl2,ytop2i,ybottom2i)

    toptopin = xintersecttoptop>0 and xintersecttoptop<xo1r
    topbottomin = xintersecttopbottom>0 and xintersecttopbottom<xo1r
    bottomtopin = xintersectbottomtop>0 and xintersectbottomtop<xo1r
    bottombottomin = xintersectbottombottom>0 and xintersectbottombottom<xo1r

    #print(toptopin,topbottomin,bottomtopin,bottombottomin)
    if (not toptopin) and (not topbottomin) and (not bottomtopin) and (not bottombottomin):

        #print(ybotout,ybotin,ytopout,ytopin,spike1['w'])
        #no intersection - either parallel spikes or single trapezoid
        if (ybotout > spike1['w']/2 and ybotin > spike1['w']/2) or (ytopout < -spike1['w']/2 and ytopin < -spike1['w']/2):
            #checked without flip
            #checked with flip
            #two bars
            return copy.deepcopy(spike1),copy.deepcopy(spike2)
        elif (ytopout < spike1['w']/2 and ytopin < spike1['w']/2) and (ybotout > -spike1['w']/2 and ybotin > -spike1['w']/2):
            #checked without flip
            #checked with flip
            #one bar
            return copy.deepcopy(spike1)
        elif ytopout > spike1['w']/2 and ytopin > spike1['w']/2 and ybotout<spike1['w']/2 and ybotin<spike1['w']/2:
            #checked without flip
            #checked with flip
            finalspike1 = {}
            finalspike1['type']='trap'
                
            x0rs1=xi1r
            y0rs1=yi1r-flip*spike1['w']/2

            finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)
                
            finalspike1['x0']=finalxs1+shiftx
            finalspike1['y0']=finalys1+shifty
            finalspike1['base']=1.0*spike1['l']
            finalspike1['a'] = flip*(ytopin-(yi1r-spike1['w']/2))
            finalspike1['b'] = flip*(ytopout -(yo1r-spike1['w']/2))
            finalspike1['theta']=1.0*spike1['theta']-np.pi/2
            return fixtrap(finalspike1)
        
        elif ytopout < spike1['w']/2 and ytopin < spike1['w']/2 and ybotout<-spike1['w']/2 and ybotin<-spike1['w']/2:
            #checked without flip
            #checked with flip
            finalspike1 = {}
            finalspike1['type']='trap'
                
            x0rs1=xo1r
            y0rs1=yo1r+flip*spike1['w']/2

            finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)
                
            finalspike1['x0']=finalxs1+shiftx
            finalspike1['y0']=finalys1+shifty
            finalspike1['base']=1.0*spike1['l']
            finalspike1['a'] = flip*((yo1r+spike1['w']/2)-ybotout)
            finalspike1['b'] = flip*((yi1r+spike1['w']/2)-ybotin)
            finalspike1['theta']=1.0*spike1['theta']-3*np.pi/2
            return fixtrap(finalspike1)
    elif (toptopin) and (not topbottomin) and (not bottomtopin) and (not bottombottomin):
        #checked without flip
        #checked with flip
        #inner bar and outer trapezoid
        finalspike1 = {}
        finalspike1['type']='bar'
                
        x0rs1=np.mean([xi1r,xintersecttoptop])
        y0rs1=yi1r

        finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)
                
        finalspike1['x0']=finalxs1+shiftx
        finalspike1['y0']=finalys1+shifty
        finalspike1['l']=xintersecttoptop - xi1r
        finalspike1['w']=1.0*spike1['w']
        finalspike1['theta']=1.0*spike1['theta']

        finalspike2 = {}
        finalspike2['type']='trap'

        x0rs2=xintersecttoptop
        y0rs2=yi1r-flip*spike1['w']/2
                
        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)
                
        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty
        finalspike2['base']=xo1r - xintersecttoptop
        finalspike2['a'] = flip*1.0*spike1['w']
        finalspike2['b'] = flip*(ytopout0 - ( yo1r-spike1['w']/2))
        finalspike2['theta']=1.0*spike1['theta']-np.pi/2

        return finalspike1,fixtrap(finalspike2)

    elif (not toptopin) and (not topbottomin) and (not bottomtopin) and (bottombottomin):

        #checked without flip
        #checked with flip
        #inner trapezoid and outer bar
        finalspike1 = {}
        finalspike1['type']='trap'
                
        x0rs1=xintersectbottombottom
        y0rs1=yi1r+flip*spike1['w']/2
        

        finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)
                
        finalspike1['x0']=finalxs1+shiftx
        finalspike1['y0']=finalys1+shifty
        finalspike1['base']=xintersectbottombottom-xi1r
                
        finalspike1['a'] = flip*1.0*spike1['w']
        finalspike1['b'] = flip*(yi1r+spike1['w']/2-ybotin0)
        finalspike1['theta']=1.0*spike1['theta'] - 3*np.pi/2

        finalspike2 = {}
        finalspike2['type']='bar'

        x0rs2 = np.mean([xintersectbottombottom,xo1r])
        y0rs2 = yo1r

        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)

        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty

        finalspike2['w']=1.0*spike1['w']
        finalspike2['l']=xo1r-xintersectbottombottom
        finalspike2['theta']=1.0*spike1['theta']
        
        
        return fixtrap(finalspike1),finalspike2

    elif (toptopin) and (not topbottomin) and (not bottomtopin) and (bottombottomin):
        #inner trapezoid, middle bar, outer trapezoid
        #checked without flip
        #checked with flip
        #print('s')
        finalspike1 = {}
        finalspike1['type']='trap'
                
        x0rs1=xintersectbottombottom
        y0rs1=yi1r+flip*spike1['w']/2
        

        finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)
                
        finalspike1['x0']=finalxs1+shiftx
        finalspike1['y0']=finalys1+shifty
        finalspike1['base']=xintersectbottombottom-xi1r
                
        finalspike1['a'] = flip*1.0*spike1['w']
        finalspike1['b'] = flip*(yi1r+spike1['w']/2-ybotin0)
        finalspike1['theta']=1.0*spike1['theta'] + np.pi/2

        finalspike2 = {}
        finalspike2['type']='bar'

        x0rs2 = np.mean([xintersectbottombottom,xintersecttoptop])
        y0rs2 = yo1r

        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)

        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty

        finalspike2['w']=1.0*spike1['w']
        finalspike2['l']=xintersecttoptop-xintersectbottombottom
        finalspike2['theta']=1.0*spike1['theta']

        finalspike3 = {}
        finalspike3['type']='trap'
                
        x0rs3=xintersecttoptop
        y0rs3=yi1r-flip*spike1['w']/2
        

        finalxs3,finalys3=getrot(x0rs3,y0rs3,-rotangle)
                
        finalspike3['x0']=finalxs3+shiftx
        finalspike3['y0']=finalys3+shifty
        finalspike3['base']=xo2r-xintersecttoptop
                
        finalspike3['a'] = flip*1.0*spike1['w']
        finalspike3['b'] = flip*(ytopout0- (yo1r-spike1['w']/2))
        finalspike3['theta']=1.0*spike1['theta']-np.pi/2
        return fixtrap(finalspike1),finalspike2,fixtrap(finalspike3)

    elif (not toptopin) and (topbottomin) and (not bottomtopin) and (not bottombottomin):
        #inner trapezoid, outer bar (from spike 1), and outer trapezoid (the outer part of spike 2 becomes a thin trapezoid because the intersection with spike 1 is not at 90 deg)
        #checked without flip
        #checked with flip (minor misalignment may be due to imperfect xin)

        finalspike1 = {}
        finalspike1['type']='trap'

        x0rs1 = xi1r
        y0rs1 = yi1r - flip*spike1['w']/2
        
        finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)

        finalspike1['x0']=finalxs1+shiftx
        finalspike1['y0']=finalys1+shifty

        finalspike1['base']=xintersecttopbottom-xi1r
        finalspike1['a'] = flip*(ytopin0 - (yi1r - spike1['w']/2))
        finalspike1['b'] = flip*(m2*(xintersecttopbottom-xi2r)+ytopin0-(yi1r - spike1['w']/2))

        finalspike1['theta'] = 1.0*spike1['theta']-np.pi/2

        finalspike2 = {}
        finalspike2['type']='bar'

        x0rs2 = np.mean([xintersecttopbottom,xo1r])
        y0rs2 = yo1r

        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)
        
        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty

        finalspike2['l']=xo1r-xintersecttopbottom
        finalspike2['w']=1.0*spike1['w']
        finalspike2['theta']=1.0*spike1['theta']

        
        finalspike3 = {}
        finalspike3['type']='trap'

        x0rs3 = xbotout
        y0rs3 = flip*ybotout0

        finalxs3,finalys3=getrot(x0rs3,y0rs3,-rotangle)

        finalspike3['x0']=finalxs3+shiftx
        finalspike3['y0']=finalys3+shifty

        finalspike3['base']=flip*1.0*spike2['w']
        finalspike3['a']= np.sqrt((xbotout - xintersecttopbottom)**2 + (ybotout0 - (yo1r+spike1['w']/2))**2)

        ytopint = m2*(xintersecttopbottom-xi2r)+ytopin0
        finalspike3['b']= np.sqrt((xtopout - xintersecttopbottom)**2 + (ytopout0 - ytopint)**2)

        finalspike3['theta'] = spike2['theta']

        return fixtrap(finalspike1),finalspike2,fixtrap(finalspike3)

    elif (not toptopin) and (not topbottomin) and (bottomtopin) and (not bottombottomin):
        #inner trapezoid of rotated spike2, inner bar (spike 1), and outer trapezoid
        #checked without flip
        #checked with flip
        #print('t')
        finalspike1 = {}
        finalspike1['type']='trap'

        x0rs1 = xtopin
        y0rs1 = flip*ytopin0

        finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)

        finalspike1['x0']=finalxs1+shiftx
        finalspike1['y0']=finalys1+shifty

        finalspike1['base']=flip*1.0*spike2['w']
        finalspike1['a']= np.sqrt((xtopin - xintersectbottomtop)**2 + (ytopin0 - (yo1r-spike1['w']/2))**2)

        ybotint = m2*(xintersectbottomtop-xi2r)+ybotin
        finalspike1['b']= np.sqrt((xbotin - xintersectbottomtop)**2 + (ybotin0 - ybotint)**2)

        finalspike1['theta'] = spike2['theta']+np.pi

        finalspike2 = {}
        finalspike2['type']='bar'

        x0rs2 = np.mean([xi1r,xintersectbottomtop])
        y0rs2 = yi1r

        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)

        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty

        finalspike2['l']=xintersectbottomtop-xi1r
        finalspike2['w']=1.0*spike1['w']

        finalspike2['theta']=1.0*spike1['theta']

        finalspike3 = {}
        finalspike3['type']='trap'

        x0rs3 = xo1r
        y0rs3 = yo1r + flip*spike1['w']/2

        finalxs3,finalys3=getrot(x0rs3,y0rs3,-rotangle)

        finalspike3['x0']=finalxs3+shiftx
        finalspike3['y0']=finalys3+shifty

        finalspike3['base']=xo1r - xintersectbottomtop
        finalspike3['a']=flip*(yo1r + spike1['w']/2 - ybotout0)
        finalspike3['b']=flip*(yo1r + spike1['w']/2 - (m2*(xintersectbottomtop-xi2r)+ybotin0))

        finalspike3['theta']=1.0*spike1['theta']+np.pi/2
        return fixtrap(finalspike1),finalspike2,fixtrap(finalspike3)

    elif (toptopin) and (topbottomin) and (not bottomtopin) and (not bottombottomin):
        #bar 1 remains unchanged. the remnants of bar2 are a tilted trapezoid
        #checked without flip
        #print('a')

        finalspike1 = copy.deepcopy(spike1)
        

        #TODO
        #the a and b calculations of these seem to have very slight offsets
        finalspike2 = {}
        finalspike2['type']='trap'

        x0rs2 = xbotout
        y0rs2 = flip*ybotout0

        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)

        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty

        finalspike2['base']=flip*1.0*spike2['w']
        finalspike2['a'] = np.sqrt((xbotout - xintersecttopbottom)**2 + (ybotout0-(yo1r+spike1['w']/2))**2)
        finalspike2['b'] = np.sqrt((xtopout - xintersecttoptop)**2 + (ytopout0 - (yo1r+spike1['w']/2))**2)

        finalspike2['theta'] = spike2['theta']

        return finalspike1,fixtrap(finalspike2)

    elif (not toptopin) and (not topbottomin) and (bottomtopin) and (bottombottomin):
        #bar 1 remains unchanged. the remnants of bar2 are a tilted trapezoid
        #checked without flip
        #checked with flip
        #print('b')
        

        finalspike1 = copy.deepcopy(spike1)
        

        finalspike2 = {}
        finalspike2['type']='trap'

        x0rs2 = xtopin
        y0rs2 = flip*ytopin0

        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)

        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty

        finalspike2['base']=flip*1.0*spike2['w']
        finalspike2['a'] = np.sqrt((xtopin - xintersectbottomtop)**2 + (ytopin0-(yo1r-spike1['w']/2))**2)
        finalspike2['b'] = np.sqrt((xbotin - xintersectbottombottom)**2 + (ybotin0 - (yo1r-spike1['w']/2))**2)

        finalspike2['theta'] = spike2['theta']+np.pi

        return finalspike1,fixtrap(finalspike2)

    elif (toptopin) and (not topbottomin) and (bottomtopin) and (bottombottomin):
        #bar 1 is shortened to the end part, which becomes another trapezoid. the remnants of bar2 are a two tilted trapezoids
        #checked without flip
        #checked with flip
        #print('b')
        

        finalspike1 = {}
        finalspike1['type'] = 'bar'

        x0rs1 = np.mean([xi1r,xintersecttoptop])
        y0rs1 = yo1r
        
        finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)

        finalspike1['x0'] = finalxs1+shiftx
        finalspike1['y0'] = finalys1+shifty
        
        finalspike1['w'] = 1.0*spike1['w']
        finalspike1['l']=xintersecttoptop-xi1r
        finalspike1['theta']=1.0*spike1['theta']
        
        finalspike2 = {}
        finalspike2['type']='trap'

        x0rs2 = xtopin
        y0rs2 = flip*ytopin0

        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)

        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty

        finalspike2['base']=flip*1.0*spike2['w']
        finalspike2['a'] = np.sqrt((xtopin - xintersectbottomtop)**2 + (ytopin0-(yo1r-spike1['w']/2))**2)
        finalspike2['b'] = np.sqrt((xbotin - xintersectbottombottom)**2 + (ybotin0 - (yo1r-spike1['w']/2))**2)

        finalspike2['theta'] = spike2['theta']+np.pi

        finalspike3 = {}
        finalspike3['type']='trap'

        x0rs3 = xo1r
        y0rs3 = yi1r - flip*spike1['w']/2
        
        finalxs3,finalys3=getrot(x0rs3,y0rs3,-rotangle)

        finalspike3['x0']=finalxs3+shiftx
        finalspike3['y0']=finalys3+shifty
        
        finalspike3['a'] = -1*flip*(ytopout0- (yo1r-spike1['w']/2))
        finalspike3['base'] = xo1r - xintersecttoptop
        finalspike3['b'] = -1*flip*1.0*spike1['w']
        finalspike3['theta']=1.0*spike1['theta']+np.pi/2

        return finalspike1,fixtrap(finalspike2),fixtrap(finalspike3)

    elif (toptopin) and (topbottomin) and (not bottomtopin) and (bottombottomin):
        #bar 1 is shortened to the end part, which becomes another trapezoid. the remnants of bar2 are a two tilted trapezoids
        #checked without flip
        #checked with flip


        finalspike1 = {}
        finalspike1['type'] = 'bar'

        x0rs1 = np.mean([xo1r,xintersectbottombottom])
        y0rs1 = yo1r
        
        finalxs1,finalys1=getrot(x0rs1,y0rs1,-rotangle)

        finalspike1['x0'] = finalxs1+shiftx
        finalspike1['y0'] = finalys1+shifty
        
        finalspike1['w'] = 1.0*spike1['w']
        finalspike1['l']=xo1r - xintersectbottombottom
        finalspike1['theta']=1.0*spike1['theta']
        finalspike2 = {}
        finalspike2['type']='trap'

        x0rs2 = xbotout
        y0rs2 = flip*ybotout0

        finalxs2,finalys2=getrot(x0rs2,y0rs2,-rotangle)

        finalspike2['x0']=finalxs2+shiftx
        finalspike2['y0']=finalys2+shifty

        finalspike2['base']=flip*1.0*spike2['w']
        finalspike2['a'] = np.sqrt((xbotout - xintersecttopbottom)**2 + (ybotout0-(yo1r+spike1['w']/2))**2)
        finalspike2['b'] = np.sqrt((xtopout - xintersecttoptop)**2 + (ytopout0 - (yo1r+spike1['w']/2))**2)

        finalspike2['theta'] = spike2['theta']


        finalspike3 = {}
        finalspike3['type']='trap'
                
        x0rs3=xintersectbottombottom
        y0rs3=yi1r+flip*spike1['w']/2
        

        finalxs3,finalys3=getrot(x0rs3,y0rs3,-rotangle)
                
        finalspike3['x0']=finalxs3+shiftx
        finalspike3['y0']=finalys3+shifty
        finalspike3['base']=xintersectbottombottom-xi1r
                
        finalspike3['a'] = flip*1.0*spike1['w']
        finalspike3['b'] = flip*(yi1r+spike1['w']/2-ybotin0)
        finalspike3['theta']=1.0*spike1['theta'] + np.pi/2


        return finalspike1,fixtrap(finalspike2),fixtrap(finalspike3)

    else:
        print('unexpected case')
        return None
        

        
    #     xintersecttoptop=(spike1['w']/2-ytopi+m2*xi2r)/m2 #top of spike 2 and top of spike 1
    #     xintersecttopbottom=(spike1['w']/2-ybottomi+m2*xi2r)/m2 #bottom of spike 2 and top of spike 1
    #     xintersectbottomtop=(-spike1['w']/2-ytopi+m2*xi2r)/m2 #top of spike 2 and bottom of spike 1
    #     xintersectbottombottom=(-spike1['w']/2-ybottomi+m2*xi2r)/m2 #bottom of spike 2 and bottom of spike 1

    #     ytoptop=m2*(xintersecttoptop-xi2r)+yi2r
    #     ytopbottom=m2*(xintersecttopbottom-xi2r)+yi2r
    #     ybottomtop=m2*(xintersectbottomtop-xi2r)+yi2r
    #     ybottombottom=m2*(xintersectbottombottom-xi2r)+yi2r


    #     fullspikes=[]
    #     if spike2['theta']<spike1['theta']:

    #         if xintersectbottomtop<xintersecttopbottom:
    #             print('warning, weird shape, no combination')
            
    #         if xintersecttopbottom>xi1r:
    #             #two spikes of inner edges of spikes
    #             finalspike1={}
    #             finalspike1['type']='bar'

    #             x0rs1=np.mean([xi1r,xintersecttopbottom])
    #             y0rs1=spike1['y0']

    #             finalxs1=x0rs1*np.cos(-rotangle)+y0rs1*np.sin(-rotangle)
    #             finalys1=-x0rs1*np.sin(-rotangle)+y0rs1*np.cos(-rotangle)

    #             finalspike1['x0']=finalxs1
    #             finalspike1['y0']=finalys1
    #             finalspike1['l']=xintersecttopbottom-xi1r
    #             finalspike1['w']=spike1['w']
    #             finalspike1['theta']=spike1['theta']

    #             fullspikes.append(finalspike1)

    #             finalspike2={}
    #             finalspike2['type']='trap'
    #             x0rs2=xtop2i
    #             y0rs2=ytop2i

    #             finalxs2=x0rs2*np.cos(-rotangle)+y0rs2*np.sin(-rotangle)
    #             finalys2=-x0rs2*np.sin(-rotangle)+y0rs2*np.cos(-rotangle)

    #             finalspike2['x0']=finalxs2
    #             finalspike2['y0']=finalys2
    #             finalspike2['base']=spike2['w']
    #             finalspike2['a']=np.sqrt((ytop2i-ytopi)**2+(xtop2i-xintersecttoptop)**2)
    #             finalspike2['b']=np.sqrt((ybottom2i-ytopi)**2+(xbottom2i-xintersecttopbottom)**2)
    #             finalspike2['theta']=spike2['theta']-np.pi/2

    #             fullspikes.append(finalspikes2)

    #         else:
    #             if xintersecttoptop>xi1r and xintersecttoptop<xo1r:
                    
    #                 finalspike3={}
    #                 finalspike3['type']='trap'
                    
    #                 x0rs3=xi1r
    #                 y0rs3=yi1r-spike1['w']/2

    #                 finalxs3=x0rs3*np.cos(-rotangle)+y0rs3*np.sin(-rotangle)
    #                 finalys3=-x0rs3*np.sin(-rotangle)+y0rs3*np.cos(-rotangle)

    #                 finalspike3['x0']=finalxs3
    #                 finalspike3['y0']=finalys3
    #                 finalspike3['base']=xintersecttoptop-xi1r
    #                 finalspike3['a']=ytop2i-(yi1r-spike1['w']/2)
    #                 finalspike3['b']=spike1['w']
    #                 finalspike3['theta']=spike1['theta']
                    
    #                 fullspikes.append(finalspike3)

    #             elif xintersecttoptop>xo1r:
    #                 finalspike4={}
    #                 finalspike4['type']='trap'

    #                 x0rs4=xi1r
    #                 y0rs4=yi1r-spike1['w']/2

    #                 finalxs4=x0rs4*np.cos(-rotangle)+y0rs4*np.sin(-rotangle)
    #                 finalys4=-x0rs4*np.sin(-rotangle)+y0rs4*np.cos(-rotangle)

    #                 finalspike4['x0']=finalxs3
    #                 finalspike4['y0']=finalys3
    #                 finalspike4['base']=xo1r-xi1r
    #                 finalspike4['a']=ytop2i-(yi1r-spike1['w']/2)
    #                 finalspike4['b']=ytop2o-(yo1r-spike1['w']/2)
    #                 finalspike4['theta']=spike1['theta']
                    
    #                 fullspikes.append(finalspike4)
                    
    #         if xintersecttoptop<xintersectbottombottom:
    #             finalspike5={}
    #             finalspike5['type']='bar'
    #         finalspike3={}
    #         finalspike3['type']='trap'
            
    #         x0rs3=np.max([xi1r,xintersecttopbottom])
    #         y0rs3=yi1r-spike1['w']/2
            
    #         finalxs3=x0rs3*np.cos(-rotangle)+y0rs3*np.sin(-rotangle)
    #         finalys3=-x0rs3*np.sin(-rotangle)+y0rs3*np.cos(-rotangle)

            
    #         finalspike3['x0']=finalxs3
    #         finalspike3['y0']=finalys3
    #         finalspike3['base']=xintersecttoptop-xi1r
    #         finalspike3['a']=ytop2i-(yi1r-spike1['w']/2)
    #         finalspike3['b']=spike1['w']
    #         finalspike3['theta']=spike1['theta']

    #         finalspike4={}
    #         finalspike4['type']='bar'

    #         x0rs4=np.mean([xintersecttopbottom,xintersectbottomtop])
    #         y0rs4=yi1r

    #         finalxs4=x0rs4*np.cos(-rotangle)+y0rs4*np.sin(-rotangle)
    #         finalys4=-x0rs4*np.sin(-rotangle)+y0rs4*np.cos(-rotangle)
            
    #         finalspike4['x0']=finalxs4
    #         finalspike4['y0']=finalys4
    #         finalspike4['l']=xintersecttopbottom-xintersectbottomtop
    #         finalspike4['w']=spike1['w']
    #         finalspike4['theta']=spike1['theta']

    #         finalspike5={}
    #         finalspike5['type']='trap'

    #         x0rs5=np.min([xo1r,xintersectbottomtop])
    #         y0rs5=yo1r+spike1['w']/2

    #         finalxs5=x0rs5*np.cos(-rotangle)+y0rs5*np.sin(-rotangle)
    #         finalys5=-x0rs5*np.sin(-rotangle)+y0rs5*np.cos(-rotangle)

    #         finalspike5['x0']=finalxs5
    #         finalspike5['y0']=finalys5
    #         finalspike5['base']=xintersectbottomtop-xintersectbottombottom
    #         finalspike5['a']=y1ir+spike1['w']/2-ybottom2o
    #         finalspike5['b']=spike1['w']
    #         finalspike5['theta']=spike1['theta']-np.pi


                
        
    #     xintersects=[xi1r,xintersecttoptop,xintersecttopbottom,xintersectbottomtop,xintersectbottombottom,xo1r]
        
    #     intersectsort=np.argsort(xintersects)
        
    #     finalspike1=copy.deepcopy(spike1)
    #     finalspike1['x0']=np.mean([xintersects[intersectsort[0]],xintersects[intersectsort[
        
    # elif abs(yi1r-yi2r)<(spike1['w']/2+spike2['w']/2) and abs(yo1r-yo2r)<(spike1['w']/2+spike2['w']/2):
    #     #completely merged into the one spike
    #     finalspike1=copy.deepcopy(spike1)
    #     finalxir=np.min([xi1r,xi2r])
        
    #     if yi1r<yi2r:
    #         #spike 2 on top
    #         finalyir=np.mean([yi1r-spike1['w']/2,yi2r+appw2/2])
    #         winner=(yi2r+appw2/2)-(yi1r-spike1['w']/2)
    #     else:
    #         finalyir=np.mean([yi1r+spike1['w']/2,yi2r-appw2/2])
    #         winner=(yi1r+spike1['w']/2)-(yi2r-appw2)
            
    #     finalxor=np.max([xo1r,xo2r])
    #     if yo1r<yo2r:
    #         #spike 2 on top
    #         finalyor=np.mean([yo1r-spike1['w']/2,yo2r+appw2/2])
    #         wouter=(yo2r+appw2/2)-(yo1r-spike1['w']/2)
    #     else:
    #         finalyor=np.mean([yo1r+spike1['w']/2,yo2r-appw2/2])
    #         wouter=(yo1r+spike2['w']/2)-(yo2r-appw2/2)

    #     if yc1r<yc2r:
    #         #spike 2 on top
    #         finalycr=np.mean([yc1r-spike1['w']/2,yc2r+appw2/2])
    #         wc=(yc2r+appw2/2)-(yc1r-spike1['w']/2)
    #     else:
    #         finalyor=np.mean([yc1r+spike1['w']/2,yc2r-appw2/2])
    #         wc=(yc1r+spike2['w']/2)-(yc2r-appw2/2)

    #     finalxr=np.mean([finalxor,finalxir])
    #     finalyr=np.mean([finalyor,finalyir])
    #     finalw=np.mean([winner,wouter,wc,wc]) #wc in to deal with the possibility that spike 2 crosses over
    #     finall=np.sqrt((finalyor-finalyir)**2+(finalxor-finalxir)**2)
    #     finalx=finalxr*np.cos(-rotangle)+finalyr*np.sin(-rotangle)
    #     finaly=-finalxr*np.sin(-rotangle)+finalyr*np.cos(-rotangle)

    #     finaldeltatheta=np.arctan2(finalyor-finalyir,finalxor-finalyir)
        
    #     finalspike1['x0']=finalx
    #     finalspike1['y0']=finaly
    #     finalspike1['w']=finalw
    #     finalspike1['l']=finall
    #     finalspike1['theta']=rotangle+finaldeltatheta+np.pi/2

    #     return (finalspike1,)

    # else:
    #     #starts as one spike, but splits into two spikes
    #     appw2=spike2['w']/np.cos(spike2['theta']-spike1['theta'])

    #     if yo1r<yo2r:
    #         #spike2 on top of spike1
    #         xintersect=(yo1r+spike1['w']/2-(yo2r-appw2/2)+m2*(xo2r+appw2/2*np.sin(spike2['theta']-spike1['theta'])))/m2
            
    #     else:
    #         #spike1 below spike2
    #         xintersect=(yo1r-spike1['w']/2-(yo2r+appw2/2)+m2*(xo2r-appw2/2*np.sin(spike2['theta']-spike1['theta'])))/m2

    #     y2int=yi2r+m2*(xintersect-xi2r)

    #     if abs(yi1r-yi2r)<(spike1['w']/2+spike2['w']/2) and abs(yo1r-yo2r)>(spike1['w']/2+spike2['w']/2):
    #         #most likely case - 1 spike inner and 2 spikes outer
    #         finalxir=np.min([xi1r,xi2r])
    #         centerx3=np.mean([finalxir,xintersect])

    #         y3=yi2r+m2*(centerx3-xi2r)
            
    #         if yi1r<yi2r:
    #             #spike 2 on top
    #             finalyir=np.mean([yi1r-spike1['w']/2,yi2r+appw2/2])
    #             winner=(yi2r+appw2/2)-(yi1r-spike1['w']/2)
    #         else:
    #             finalyir=np.mean([yi1r+spike1['w']/2,yi2r-appw2/2])
    #             winner=(yi1r+spike1['w']/2)-(yi2r-appw2/2)
            
    #         if yo1r<y2int:
    #             #spike 2 on top
    #             finalyint=np.mean([yo1r-spike1['w']/2,y2int+appw2/2])
    #             wouter=(y2int+appw2/2)-(yo1r-spike1['w']/2)
    #         else:
    #             finalyint=np.mean([yo1r+spike1['w']/2,y2int-appw2/2])
    #             wouter=(yo1r+spike2['w']/2)-(y2int-appw2/2)

    #         if yc1r<y3:
    #             #spike 2 on top
    #             wc=(y3+appw2/2)-(yc1r-spike1['w']/2)
    #         else:
    #             wc=(yc1r+spike2['w']/2)-(y3-appw2/2)

    #         finalxr3=np.mean([xintersect,finalxir])
    #         finalyr3=np.mean([finalyint,finalyir])
    #         finalw3=np.mean([winner,wouter,wc,wc])
    #         finall3=np.sqrt((finalyint-finalyir)**2+(xintersect-finalxir)**2)
    #         finalx3=finalxr3*np.cos(-rotangle)+finalyr3*np.sin(-rotangle)
    #         finaly3=-finalxr3*np.sin(-rotangle)+finalyr3*np.cos(-rotangle)

    #         finaldeltatheta3=np.arctan2(finalyint-finalyir,xintersect-finalxir)

    #         finalspike3=copy.deepcopy(spike1)
    #         finalspike3['x0']=finalx3
    #         finalspike3['y0']=finaly3
    #         finalspike3['w']=finalw3
    #         finalspike3['l']=finall3
    #         finalspike3['theta']=rotangle+finaldeltatheta3+np.pi/2


    #         finalxr1=np.mean([xintersect,xo1r])
    #         finalyr1=yc1r
    #         finalw1=spike1['w']
    #         finall1=xo1r-xintersect
    #         finalx1=finalxr1*np.cos(-rotangle)+finalyr1*np.sin(-rotangle)
    #         finaly1=-finalxr1*np.sin(-rotangle)+finalyr1*np.cos(-rotangle)
            
    #         finalspike1=copy.deepcopy(spike1)
    #         finalspike1['x0']=finalx1
    #         finalspike1['y0']=finaly1
    #         finalspike1['w']=finalw1
    #         finalspike1['l']=finall1
    #         finalspike1['theta']=spike1['theta']

    #         finalxr2=np.mean([xintersect,xo2r])
    #         finalyr2=np.mean([y2int,yo2r])
    #         finalw2=spike2['w']
    #         finall2=np.sqrt((yo2r-y2int)**2+(xo2r-xintersect)**2)
    #         finalx2=finalxr2*np.cos(-rotangle)+finalyr2*np.sin(-rotangle)
    #         finaly2=-finalxr2*np.sin(-rotangle)+finalyr2*np.cos(-rotangle)
            
    #         finalspike2=copy.deepcopy(spike2)
    #         finalspike2['x0']=finalx2
    #         finalspike2['y0']=finaly2
    #         finalspike2['w']=finalw2
    #         finalspike2['l']=finall2
    #         finalspike2['theta']=spike2['theta']

    #         return finalspike1,finalspike2,finalspike3
        
    #     else:
    #         #less likely case - 2 spikes inner and 1 spike outer

    #         finalxor=np.max([xo1r,xo2r])
    #         centerx3=np.mean([finalxor,xintersect])
    #         y3=yi2r+m2*(centerx3-xi2r)
            
    #         if yo1r<yo2r:
    #             #spike 2 on top
    #             finalyor=np.mean([yo1r-spike1['w']/2,yo2r+appw2/2])
    #             wouter=(yo2r+appw2/2)-(yo1r-spike1['w']/2)
    #         else:
    #             finalyor=np.mean([yo1r+spike1['w']/2,yo2r-appw2/2])
    #             wouter=(yo1r+spike1['w']/2)-(yo2r-appw2/2)
            
    #         if yi1r<y2int:
    #             #spike 2 on top
    #             finalyint=np.mean([yi1r-spike1['w']/2,y2int+appw2/2])
    #             wint=(y2int+appw2/2)-(yo1r-spike1['w']/2)
    #         else:
    #             finalyint=np.mean([yo1r+spike1['w']/2,y2int-appw2/2])
    #             wint=(yo1r+spike2['w']/2)-(y2int-appw2/2)
            
    #         if yc1r<y3:
    #             #spike 2 on top
    #             wc=(y3+appw2/2)-(yc1r-spike1['w']/2)
    #         else:
    #             wc=(yc1r+spike2['w']/2)-(y3-appw2/2)

            
    #         finalxr3=np.mean([xintersect,finalxor])
    #         finalyr3=np.mean([finalyint,finalyor])
    #         finalw3=np.mean([wint,wouter,wc,wc])
    #         finall3=np.sqrt((finalyint-finalyor)**2+(xintersect-finalxor)**2)
    #         finalx3=finalxr3*np.cos(-rotangle)+finalyr3*np.sin(-rotangle)
    #         finaly3=-finalxr3*np.sin(-rotangle)+finalyr3*np.cos(-rotangle)

    #         finaldeltatheta3=np.arctan2(finalyor-finalyint,finalxor-xintersect)

    #         finalspike3=copy.deepcopy(spike1)
    #         finalspike3['x0']=finalx3
    #         finalspike3['y0']=finaly3
    #         finalspike3['w']=finalw3
    #         finalspike3['l']=finall3
    #         finalspike3['theta']=rotangle+finaldeltatheta3+np.pi/2


    #         finalxr1=np.mean([xintersect,xi1r])
    #         finalyr1=yc1r
    #         finalw1=spike1['w']
    #         finall1=xintersect-xi1r
    #         finalx1=finalxr1*np.cos(-rotangle)+finalyr1*np.sin(-rotangle)
    #         finaly1=-finalxr1*np.sin(-rotangle)+finalyr1*np.cos(-rotangle)
            
    #         finalspike1=copy.deepcopy(spike1)
    #         finalspike1['x0']=finalx1
    #         finalspike1['y0']=finaly1
    #         finalspike1['w']=finalw1
    #         finalspike1['l']=finall1
    #         finalspike1['theta']=spike1['theta']

    #         finalxr2=np.mean([xintersect,xi2r])
    #         finalyr2=np.mean([y2int,yi2r])
    #         finalw2=spike2['w']
    #         finall2=np.sqrt((yi2r-y2int)**2+(xi2r-xintersect)**2)
    #         finalx2=finalxr2*np.cos(-rotangle)+finalyr2*np.sin(-rotangle)
    #         finaly2=-finalxr2*np.sin(-rotangle)+finalyr2*np.cos(-rotangle)
            
    #         finalspike2=copy.deepcopy(spike2)
    #         finalspike2['x0']=finalx2
    #         finalspike2['y0']=finaly2
    #         finalspike2['w']=finalw2
    #         finalspike2['l']=finall2
    #         finalspike2['theta']=spike2['theta']

    #         return finalspike1,finalspike2,finalspike3
        
