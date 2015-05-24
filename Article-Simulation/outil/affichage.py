'''
Created on 2 fevr. 2012

@author: Geoffrey
'''
import calcul
from brian import *
import os
import numpy as np
from scipy.interpolate import interp1d
from PIL import Image, ImageOps
from matplotlib import cm
import matplotlib.pyplot as plt

def surface3D(Z, width):
    x = np.arange(0, width, 1)
    y = np.arange(0, width, 1)
    X,Y = np.meshgrid(x, y)
    #imshow()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
            linewidth=0, antialiased=False)
    plt.show()

def statePlot(Mv, epi, elec_pos):
    plot(Mv.times/ms, Mv[epi]/mV)
    xlabel('Time in ms')
    ylabel('Menbrane potential in mV')
    title('Neurone '+str(epi) + ' loc: ' +str(elec_pos))


def rastPlot(Ms):
    raster_plot(Ms, showlast=200*ms)
    xlabel('Time in ms')
    ylabel('Neurons')
    title('Spike distribution over neurons and time')

def stimPlot(patterns, width):
    pattern = patterns[0]
    for i in xrange(1, len(patterns)):
        pattern += patterns[i]
    imshow(pattern.reshape((width,width)), interpolation='nearest', origin='lower')
    #imshow(CLi.W.todense(), interpolation='nearest', origin='lower')
    xlabel('Neurons')
    ylabel('Neurons')
    title('Carte du SCs/ Stimulation Visuelle')

def avSpikeplot(Cs, sim_time, width):
    Cs_reshape = Cs.count
    Cs_reshape = Cs_reshape.reshape((width,width))
    im = imshow(Cs_reshape/sim_time, interpolation='nearest', origin='lower',vmin = 0, vmax = 660)
    cbar2 = plt.colorbar(im, orientation='vertical', shrink=0.5, aspect=15)
    xlabel('X axis (Neurons)')
    ylabel('Y axis (Neurons)')
    title('Average of Spikes during Sim')

def avStatePlot(Mv, sim_time, width):
    Mv_reshape = Mv[:]/mV
    #print np.shape(Mv_reshape)
    Mv_reshape = np.sum(Mv_reshape, axis = 1)/(len(Mv[0])) ## 10000 = 1/dt dans la config normale de BRIAN 
    #print np.shape(Mv_reshape)
    Mv_reshape = Mv_reshape.reshape((width,width))
    im = imshow(Mv_reshape, interpolation='nearest', origin='lower',vmin = -80, vmax = -40)
    cbar2 = plt.colorbar(im, orientation='vertical', shrink=0.5, aspect=15)
    xlabel('X axis (Neurons)')
    ylabel('Y axis (Neurons)')
    title('Average of Membrane Potential during Sim')
    return(Mv_reshape)

def openImg(path):
    Img = Image.open(str(path))
    Img1 = ImageOps.grayscale(Img)
    largeur,hauteur = Img1.size
    imdata=Img1.getdata()
    tab=np.array(imdata)
    matrix = np.reshape(tab,(hauteur,largeur))
    return matrix

def saveImg(map, path):
    copie = Image.new("L",(map.shape[1],map.shape[0]))
    copie.putdata(list(map.flat))
    copie.save(path, format = "png")
 
def saveIEPSPDiff(Mge, Mgi, bin, width, adpath= ""):
    print "IEPSPDiff mapping ..."
    print "Enregistrement images..."
    Mgei_t = Mge[:] - Mgi[:]
    nb_frame = Mgei_t.shape[1]
    
    Mgei_t = Mgei_t.reshape(width, width, nb_frame)
    print "Max EIPSPDiff: " + str(np.max(Mgei_t))
    print np.shape(Mgei_t)
    dir = 'IEPSP/'+adpath
    if not os.path.isdir(dir):
            os.mkdir(dir)
    for i in xrange(0, nb_frame-bin, bin):
        image = sum(Mgei_t[:,:,i:i+bin],2)/bin
        print "\r image", str(i/bin), "sur 200",
        path = dir+'outfile'+str(i)+'.png'
        imsave(path, image, format='png', origin='lower', vmin = -40, vmax = 20)
    print "Enregistrement fini."
    
def saveMembPotCut(pos, on_x, Mv, bin, width,adpath= "", adname="",r=False):
    print "Membrane potential mapping ..."
    print "Enregistrement images..."
    Mv_time = Mv[:]/mV
    nb_frame = Mv_time.shape[1]
    Mv_time = Mv_time.reshape(width, width, nb_frame)
    print np.shape(Mv_time)
    dir = 'cut-pot/'+adpath
    if not os.path.isdir(dir):
            os.mkdir(dir)
    matrice = np.zeros(((nb_frame-bin)/bin, width))
    for i in xrange(0, nb_frame-bin, bin):
        image = sum(Mv_time[:,:,i:i+bin],2)/bin
        if on_x:
            matrice[i/bin,:] = image[pos,:]
        else:
            matrice[i/bin,:] = image[:,pos]
    path = dir + adname + "act3D.txt"
    np.savetxt(path, matrice, fmt = "%10.5f")
    print "Enregistrement fini."   

def saveCount(Ms, width, time):
    count = np.zeros((width, width))
    for i in xrange(0, width**2):
        if (len(Ms[i]) == 0):
            count[i%width, i/width] = 0
        else:
            Ci = Ms[i]
            Citime = Ci[Ci>time]
            if (len(Citime) > 0):
                count[i%width, i/width] = len(Citime)/time
    return(count)
            
def saveMembPotDyn(Mv, bin, width,adpath= "",r=False): ## bin is in frame
    print "Membrane potential mapping ..."
    print "Enregistrement images..."
    Mv_time = Mv[:]/mV
    nb_frame = Mv_time.shape[1]
    Mv_time = Mv_time.reshape(width, width, nb_frame)
    print np.shape(Mv_time)
    dir = 'potimg/'+adpath
    if not os.path.isdir(dir):
            os.mkdir(dir)
    for i in xrange(0, nb_frame-bin, bin):
        image = sum(Mv_time[:,:,i:i+bin],2)/bin
        #print "\r image", str(i/bin), "sur 200",
        path = dir+'outfile'+str(i)+'.png'
        imsave(path, image, format='png', origin='lower',vmin = -80, vmax = -40)
    print "Enregistrement fini."
    #print Mv_time[:,:,i]
    
def saveFiringRateDyn(Ms, bin, sim_time, width, adpath="", save = True): ## bin in seconde
    print "Firing rate mapping ..."
    print "Enregistrement images..."
    nb_frame = int(np.trunc(sim_time/bin))
    rate = np.zeros((width, width, nb_frame-1))
    for i in xrange(0,width**2):
        if (len(Ms[i]) == 0):
            rate[i%width, i/width, :] = 0
        else:
            #print len(Ms[i]) 
            spike_timei = np.zeros((len(Ms[i])+2))
            spike_timei[1:-1] = Ms[i]
            spike_timei[-1] = sim_time
            #print len(spike_timei)
            spike_count = np.arange(0, len(spike_timei))
            spike_count[-1] = len(spike_timei)-2
            #print len(spike_count)
            f = interp1d(spike_timei, spike_count, kind="linear")
            xnew = np.arange(0, sim_time, bin)
            freq_cum = f(xnew)
            rate[i%width, i/width, :] = np.diff(freq_cum)/bin
    if save:
        dir = 'fir2img/'+adpath
        if not os.path.isdir(dir):
                os.mkdir(dir)    
        for i in xrange(0, np.shape(rate)[2]):
            #print "image "+str(i)
            time = i*bin*1000
            path = dir + 'outfile'+str(time)+'.png'
            imsave(path, rate[:,:,i], format='png', origin='lower',vmin = 0, vmax = 666)
    print "Enregistrement fini."
    return(rate)
  

    