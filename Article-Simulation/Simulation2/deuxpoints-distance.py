'''
Created on 26 janv. 2012

@author: Geoffrey
'''


from brian import *
import numpy as np
import scipy as sp
import scipy.ndimage as snd
import time, sys, Tkinter, os
from outil import calcul
from outil import affichage as aff
import matplotlib.pyplot as plt



start = time.time() 
sim_time = 0.200   #duree de la simulation
width = 100   
elec_pos = (20, 23) ## electrode sur la barre
elec_pos2 = (30, 30) ## electrode sur la cible
distractor_potential = (1333*mV, 2000*mV, 3500*mV, 4000*mV)
target_potential = 4000 * mV  ## 4000mV poids Cible
frequence_distractor = 500 * Hz ## freq distract
frequence_target = 500 * Hz ## freq cible
sx, sy, K, B = 8.5, 8.5, 1.2, 6.0

control = (0,1) ## without or with the distractor
distances = np.arange(2,42,2)
pos_dist = (49,30) ## I write the position of the distractor here
path1 = "imgsource/distractor.png" ## the distractor position is then fixed!

Nb = width*width
taut = 10 * ms
taum = 10 * ms #100
taue = 3 * ms #30
taui = 10 * ms #70
Vt = -50 * mV
Vr = -80 * mV   ## resting value
El = -70 * mV   ## stable state
Ee = 0 * mV
Ei = -80 * mV
autocatalise = 1.0 #0.7
We = 200 * autocatalise * mV ## par defaut 200mV
Wi = 200 * autocatalise * mV ## par defaut 300mV
sigma = 0*mV
tau_n = 1*ms

eqs_elec = '''
dx/dt = freq : 1
freq : Hz
 '''
eqs = Equations('''
dv/dt = (-(v-El)-ge*(v-Ee)-gi*(v-Ei)+sigma*xi*tau_n**(0.5))/taum : volt
dge/dt = -ge/taue : 1
dgi/dt = -gi/taui : 1
''')
## v : menbrane potential
## w : spiking threshold
## ge : excitatory gate
## gi : inhibitory gate


epi = calcul.c2DTo1D(elec_pos, width) ## electrode enregistr 1D pour Brian
epi2 = calcul.c2DTo1D(elec_pos2, width) ## electrode enregistr 1D pour Brian
sim_clock = Clock(dt=0.01*ms)
state_clock = Clock(dt= 1*ms)

## Initialisation du Colliculus:
SC = NeuronGroup(N= Nb, model = eqs, threshold = Vt, reset= Vr, refractory =  1.5*ms, clock = sim_clock)
DoG = calcul.connection_map4(width,sx, sy, K, B) ## circular test
map_inhib  = DoG.copy()
map_inhib[map_inhib>0] = 0
map_excit = DoG.copy()
map_excit[map_excit<0] = 0
CLe = Connection(SC, SC, 'ge', clock = sim_clock)
CLi = Connection(SC, SC, 'gi', clock = sim_clock)
CLe.W = SparseConstructionMatrix(We*map_excit)
CLi.W = SparseConstructionMatrix(-Wi*map_inhib)

print "poids du neurone 10 sur le neurone 10: inhib: ", str(CLi.W[10,10]), "excit: ", str(CLe.W[10,10]) 

conn_time = time.time()-start
a = time.gmtime(conn_time)
print time.strftime("\ntask time for connection calculation: %H:%M:%S",a )        

print "Preparation structures and images..."
G1 = NeuronGroup(1, model = eqs_elec, threshold = 1., reset = 0, clock = sim_clock) ## N entree stim pour distract
G2 = NeuronGroup(1, model = eqs_elec, threshold = 1., reset = 0, clock = sim_clock) ## N entree stim pour cible

Mv = StateMonitor(SC, 'v', record=True, clock=state_clock)#[epi, epi2])
Ms = SpikeMonitor(SC)
Cs = SpikeCounter(SC)
net = [SC, CLe, CLi, G1, G2, Mv, Ms, Cs]

resume = np.zeros((len(distances), 4))
for dist_pot in distractor_potential:
    for c in control:
        xm, ym = np.zeros((2,4)), np.zeros((2,4)) ## def matrice pour centre de gravit
        header = "sx, sy, K, B, potential Target, potential Distractor, frequence Target, frequence Distractor, control"
        parameters = (sx, sy, K, B, target_potential, dist_pot, frequence_target, frequence_distractor, c)
        foldername = "twopoints-distance/exp1-noise/%d-%d-%d-%d-T%d-D%d-fT%d-fD%d-c%d/"%parameters
        foldername = calcul.test_folder(foldername, force=True)
        np.savetxt(foldername+"parameters.txt", np.array(parameters).T,header=header, delimiter = ',')
        stim_map1 = aff.openImg(path1).reshape((width**2)) ## decouper image pour remettre sur une ligne pour brian
        stim_map1[stim_map1<150] = 0
        stim_map1[stim_map1>=150] = c
        index_counter=0
        for j in distances:
            state_clock.reinit()
            sim_clock.reinit()
            for struct in net:
                struct.reinit()
            name = "distance-" + str(j) #+ str(index_counter)
            print name
            #reinit(True)
            SC.v = El + (rand(Nb) * (El - Vt))/2.0    # valeur initiale
            ## target map
            stim_map2 = np.zeros((width,width))
            stim_map2[pos_dist[0]:pos_dist[0]+2,pos_dist[1]+j:pos_dist[1]+2+j] = 1.0 ### cible
            stim_map2 = stim_map2.reshape((width**2))

            G1.freq = calcul.createBurst(frequence_distractor, 25, 80, 200) #par defaut: 25,50,200
            G2.freq = calcul.createBurst(frequence_target, 25, 80, 200)
            C1 = Connection(G1, SC, 'ge', clock = sim_clock)
            C2 = Connection(G2, SC, 'ge', clock = sim_clock)

            C1.W = SparseConstructionMatrix(stim_map1 * dist_pot)
            C2.W = SparseConstructionMatrix(stim_map2 * target_potential)

            run(sim_time * second, threads = 4, report = 'text', report_period = 0.5)


            fig = figure(figsize=(24,12.3), dpi=80)
            subplot(131)
            aff.stimPlot([dist_pot*stim_map1, target_potential*stim_map2], width)
            plt.subplot(132)
            pot = aff.avStatePlot(Mv, sim_time, width)
            plt.subplot(133)
            #aff.avSpikeplot(Cs, sim_time, width)
            aff.rastPlot(Ms)
            np.save(foldername+name+"-raw-mp.npy", Mv[:].reshape((width,width, len(Mv[0])))/mV)
            np.save(foldername+name+"-av-spikes.npy", Cs.count.reshape((width,width))/sim_time )
            np.save(foldername+name+"-spikes.npy", Ms.spikes)
            plt.savefig(foldername+name+".png",format='png')
            close()

            moyenne = aff.saveCount(Ms, width, 0.150)

            xm, ym = calcul.calculCentreGrav(moyenne)
            xd, xt = (pos_dist[1]+1), (pos_dist[1]+1+j)
            erreur_col= (xm - xt)
            header = "mm, md, mt, um, ud, ut, phim, phid, phit, col_error, vis_error"
            resume[index_counter, :] = (xm, xd, xt,  erreur_col)
            index_counter += 1

        np.savetxt(foldername+name+"-gravitycenter.txt", resume, fmt = "%10.5f", delimiter=',')
        comp_time = time.time()-start-conn_time
        a = time.gmtime (comp_time)
        print time.strftime("task time for simulation: %H:%M:%S",a )
