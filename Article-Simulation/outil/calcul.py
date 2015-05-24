'''
Created on 18 janv. 2012

@author: Geoffrey

'''
from brian import *
import numpy as np
from scipy.interpolate import interp1d


def createFEFActivation(sim_time, a, b, Amp, Ifn, r, dec=0):
    t = np.linspace(0, sim_time, sim_time+1)
    f =  a*(b-(t-dec))
    f[0:(b+dec)] = 0
    g = Amp/(1 + np.exp(r*(Ifn-(t-dec))))
    return (TimedArray(g+f, dt=1*ms))

def createBurst(freq_max, mu, sigma, sim_time): ## mu = centre gauss en ms, compter 3sigma pour largeur, sim_time = duree de la simulation en ms
    A = freq_max
    x = np.linspace(0, sim_time, sim_time+1)
    fr = A*np.exp((-(x-mu)**2)/(sigma**2)) #+ A*exp((-(x-mu-100)**2)/(sigma**2))
    return(TimedArray(fr, dt=1*ms))

def createGaussianStim(mu_x, mu_y, sigma, resolution):
    x = np.arange(0, resolution, 1)
    y = np.arange(0, resolution, 1)
    X,Y = np.meshgrid(x, y)
    pattern = np.exp(-(  (X-mu_x)**2 + (Y-mu_y)**2.0 )/ (2.0*sigma**2.0) )
    return(pattern.reshape((resolution**2)))

def createCircleStim(map_width, c_x, c_y, xl, yl): ## better if xl and yl are even
    x = np.arange(0, map_width, 1)
    y = np.arange(0, map_width, 1)
    X,Y = np.meshgrid(x, y)
    pattern = ( (X-c_x)**2 + ((yl/xl)*(Y-c_y))**2 ) < xl**2
    return(pattern.reshape((map_width**2))) 

def createSquareStim(mu_x, mu_y, size, resolution):
    M = np.zeros((resolution, resolution))
    M[mu_y-size/2:mu_y+size/2, mu_x-size/2:mu_x+size/2] = 1.0
    return(M)

def barreVerticale(map_width, size, pos_x): ## size doit etre divisible par 2
    pattern = np.zeros((map_width, map_width))
    bord = ((map_width-size)/2)
    pattern[bord:map_width-bord, pos_x] = 1
    return(pattern.reshape((map_width**2)))
 
def creerRectangle(map_width, x, y, xl, yl):  
    pattern = np.zeros((map_width, map_width))
    pattern[y:y+yl, x:x+xl] = 1
    return(pattern.reshape((map_width**2))) 

def freq_inst(i, sim_time, bin, M):
    # on a cree [0, t1, t2, t3, ..., tf]:
    spike_timei = np.zeros((len(M[i])+2))
    spike_timei[1:-1] = M[i]
    spike_timei[-1] = sim_time
    # on cree le vecteur parallele notant le nb de spike [0, 1, 2, 3, ...,nb_f, nb_f]
    spike_count = np.arange(0, len(spike_timei))
    spike_count[-1] = len(spike_timei)-2 ## a la fin de la sim, le nombre de spike est pareil que la derniere fois.
    
    f = interp1d(spike_timei, spike_count, kind="linear")
    xnew = np.arange(0, sim_time, bin)
    freq_cum = f(xnew)
    return(np.diff(freq_cum)/bin)

def freq_group(Ms, sim_time, bin, size):
    fm = np.zeros((sim_time/bin)-1)
    for i in xrange(0,size):
        fm += freq_inst(i,sim_time,bin, Ms)
    return(fm/(size))

def c1DTo2D(i, width):
    return (i % width, int(i)/width )

def c2DTo1D(coord, width):
    return(width*coord[1]+coord[0])

def connection_map(width, A, se): ## width doit etre un int. ## A, se doivent etre des floats
    x = np.arange(0, width**2, 1)
    y = np.arange(0, width**2, 1)
    X,Y = np.meshgrid(x, y)
    
    A = 1/(se*np.sqrt(2*np.pi))
    # exp(|x-xc|/b + |y-yc|) est l'equation de base
    # on fait les transformations 2D vers 1D avec la relation: i = width.y + x
    M = A*np.exp(-(  ((X%width-Y%width)**2) + ((X/width-Y/width)**2) )/(2*se**2) )
    return(M)

def connection_map1(width, A, se): ## width doit etre un int. ## A, se doivent etre des floats
    x = np.arange(0, width**2, 1)
    y = np.arange(0, width**2, 1)
    X,Y = np.meshgrid(x, y)
    
    #A = 1/(se*np.sqrt(2*np.pi))
    # exp(|x-xc|/b + |y-yc|) est l'equation de base
    # on fait les transformations 2D vers 1D avec la relation: i = width.y + x
    M = A*np.exp(-(  ((X%width-Y%width)**2) + ((X/width-Y/width)**2) )/(2*se**2) )
    return(M)
 
def connection_map2(width, A, se1, se2): ## width doit etre un int. ## A, se doivent etre des floats
    x = np.arange(0, width**2, 1)
    y = np.arange(0, width**2, 1)
    X,Y = np.meshgrid(x, y)
     ##/(se*np.sqrt(2*np.pi))
    # exp(|x-xc|/b + |y-yc|) est l'equation de base
    # on fait les transformations 2D vers 1D avec la relation: i = width.y + x
    M = A*np.exp(-(  ((X%width-Y%width)**2)/(2*se1**2) + ((X/width-Y/width)**2)/(2*se2**2) ) ) -0.5
    return(M)   

def connection_map3(width,se1,se2, K): # morelet
    x = np.arange(0, width**2, 1)
    y = np.arange(0, width**2, 1)
    X,Y = np.meshgrid(x, y)
    G = -((X%width-Y%width)**2)/(se1**2) - ((X/width-Y/width)**2)/(se2**2)
    M=(1 + G) * np.exp(G/K)
    return(M)

def connection_map4(width, se1, se2, K, inh): ### DoG:  width, Sx, Sy, aggrandiss inhib, amplit de inh
    x = np.arange(0, width**2, 1)
    y = np.arange(0, width**2, 1)
    X,Y = np.meshgrid(x, y)
    A1 = 1.0 + inh
    A2 = inh
     ##/(se*np.sqrt(2*np.pi))
    # exp(|x-xc|/b + |y-yc|) est l'equation de base
    # on fait les transformations 2D vers 1D avec la relation: i = width.y + x
    M1 = A1*np.exp(-(  ((X%width-Y%width)**2)/(2*se1**2) + ((X/width-Y/width)**2)/(2*se2**2) ) )
    M2 = -A2*np.exp(-(  ((X%width-Y%width)**2)/(2*(K*se1)**2) + ((X/width-Y/width)**2)/(2*(K*se2)**2) ) )
    return(M1+M2)

def connection_map4s(width, se1, se2, K, inh): ### DoG:  width, Sx, Sy, aggrandiss inhib, amplit de inh
    x = np.arange(-width/2, width/2, 1)
    y = np.arange(-width/2, width/2, 1)
    x,y = np.meshgrid(x, y)
    A1 = 1.0 + inh
    A2 = inh
    M1 = A1*np.exp(-(  ((x)**2)/(2*se1**2) + ((y)**2)/(2*se2**2) ) )
    M2 = -A2*np.exp(-(  ((x)**2)/(2*(K*se1)**2) + ((y)**2)/(2*(K*se2)**2) ) )
    return(M1+M2)

def connection_map5(width, se1, se2, K, inh): ### DoG:  width, Sx, Sy, aggrandiss inhib, amplit de inh
    x = np.arange(0, width**2, 1)
    y = np.arange(0, width**2, 1)
    X,Y = np.meshgrid(x, y)
    A1 = 1.0 + inh
    A2 = inh
     ##/(se*np.sqrt(2*np.pi))
    # exp(|x-xc|/b + |y-yc|) est l'equation de base
    # on fait les transformations 2D vers 1D avec la relation: i = width.y + x
    M1 = A1*np.exp(-(  ((X%width-Y%width)**2)/(2*se1**2) + ((X/width-Y/width)**2)/(2*se2**2) ) )
    M2 = -A2*np.exp(-(  ((X%width-Y%width)**2)/(2*(K*se1)**2) + ((X/width-Y/width)**2)/(2*(K*se2)**2) ) )
    M4 = M1+M2
    B = 1/(2*se1**2) - 1/(2*(K*se1)**2)
    D = 1/(2*se2**2) - 1/(2*(K*se2)**2)
    M3 = -((X%width-Y%width)**2)*B - ((X/width-Y/width)**2)*D
    M4[M3< np.log(A2/(A1*K**2))] = np.min(M4)
    return(M4)

def connection_map6(width, se1, se2, K, a1, a2, c): ### DoG:  width, Sx, Sy, aggrandiss inhib, amplit de inh
    x = np.arange(0, width**2, 1)
    y = np.arange(0, width**2, 1)
    X,Y = np.meshgrid(x, y)
    A1 = a1
    A2 = a2
     ##/(se*np.sqrt(2*np.pi))
    # exp(|x-xc|/b + |y-yc|) est l'equation de base
    # on fait les transformations 2D vers 1D avec la relation: i = width.y + x
    M1 = A1*np.exp(-(  ((X%width-Y%width)**2)/(2*se1**2) + ((X/width-Y/width)**2)/(2*se2**2) ) )
    M2 = -A2*np.exp(-(  ((X%width-Y%width)**2)/(2*(K*se1)**2) + ((X/width-Y/width)**2)/(2*(K*se2)**2) ) )
    return(M1+M2-c)
 
def connection_map7(width, se1, se2, K, inh, c): ### DoG:  width, Sx, Sy, aggrandiss inhib, amplit de inh
    x = np.arange(0, width**2, 1)
    y = np.arange(0, width**2, 1)
    X,Y = np.meshgrid(x, y)
    A1 = 1.0 + inh
    A2 = inh
     ##/(se*np.sqrt(2*np.pi))
    # exp(|x-xc|/b + |y-yc|) est l'equation de base
    # on fait les transformations 2D vers 1D avec la relation: i = width.y + x
    M1 = A1*np.exp(-(  ((X%width-Y%width)**2)/(2*se1**2) + ((X/width-Y/width)**2)/(2*se2**2) ) )
    M2 = -A2*np.exp(-(  ((X%width-Y%width)**2)/(2*(K*se1)**2) + ((X/width-Y/width)**2)/(2*(K*se2)**2) ) )
    return(M1+M2-c)
    
def calculCentreGrav(M):
    xm, ym = 0, 0
    for x in xrange(0, M.shape[0]):
        xm += (x+1)*np.sum(M[x,:])
    xm = xm/np.sum(M) - 1
    
    for y in xrange(0, M.shape[1]):
        ym += (y+1)*np.sum(M[:,y])
    ym = ym/np.sum(M) - 1
    
    return(xm, ym  ) 

def r(x,y):
    return (np.sqrt(x**2+y**2))

def phi(x,y):
    return ( 2*np.arctan2(y,(x + np.sqrt(x**2+y**2))) )

def x(r, phi):
    return (r*np.cos(phi))

def y(r, phi):
    return (r*np.sin(phi))
    
def uColi(r, phi):
    return( 1.4* np.log(np.sqrt(r**2+2*3.0*r*np.cos(phi) + 3.0**2)/3.0) )

def vColi(r,phi):
    return( 1.8*np.arctan2((r*np.sin(phi)),(r*np.cos(phi)+3.0))  )

def rColi(u,v):
    return( 3.0*np.sqrt( np.exp(2*u/1.4) - 2*np.exp(u/1.4)*np.cos(v/1.8) + 1 ) )

def phiColi(u,v):
    return( np.arctan2( (np.exp(u/1.4)*np.sin(v/1.8)), ( np.exp(u/1.4)*np.cos(v/1.8)-1) ) )
    #return( np.arctan( (np.exp(u/1.4)*np.sin(v/1.8)) / ( np.exp(u/1.4)*np.cos(v/1.8)-1) ) )

def toColi(x,y):
    re, phie = r(x,y), phi(x,y)
    if type(x) != tuple:
        return(  ( uColi(re,phie), vColi(re,phie) )  )
    else:
        result = np.zeros((np.shape(x)[0], np.shape(x)[1],2))
        result[:,:,0] = uColi(re,phie)
        result[:,:,1] = vColi(re,phie)
        return(result)

def fromColi(u,v):
    r, phi = rColi(u,v), phiColi(u,v)
    #print np.shape(r), np.shape(phi), np.max(r), np.max(phi)
    if (len(np.shape(u))>1):
        result = np.zeros((np.shape(u)[0], np.shape(u)[1],2))
        result[:,:,0] = x(r,phi)
        result[:,:,1] = y(r,phi)
        return(result)
    else:
        return (x(r, phi), y(r, phi))

#e6 = 3
#e7 = 3
#a = toColi(6.0,0)
#b = toColi(7.5,0)
#c = toColi(6.0-e6, 0)
#d = toColi(7.5-e7, 0)    
#
#
#print "position de la cible:", str(a), str(b)
#print "landing position", str(c), str(d)
#print "erreur 6/7.5 pour barre 0.3:", str(a[0]-c[0]), str(b[0]-d[0])
