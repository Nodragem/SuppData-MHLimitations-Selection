'''
Created on 26 janv. 2012

@author: Geoffrey
'''

from brian import *
import numpy as np
import time, sys, os
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(1, "../")
from outil import calcul
from outil import affichage as aff
import matplotlib.pyplot as plt
import shutil

## stimulus' size and shape:
choose = "line" ## square or circle
beginning = 6 #2
end =  8 #43
step = 2 
PLOT_MH = False ## plot the Mexican hat before the simulation 
print "Preparation of the simulation ..."
start = time.time() 
sim_time = 0.200  
width = 100
 
## Inputs: 
eqs_elec = '''
dx/dt = freq : 1
freq : Hz
 '''
G1 = NeuronGroup(1, model = eqs_elec, threshold = 1., reset = 0) 

## the Network:
Nb = width*width
taut = 10 * ms
taum = 10 * ms 
taue = 3 * ms 
taui = 10 * ms 
Vt = -50*mV
Vr = -80 * mV   ## resting value
El = -70 * mV   ## stable state
Ee = 0 * mV
Ei = -80 * mV
autocatalise = 1.0 
We = 200 * autocatalise * mV 
Wi = 200 * autocatalise * mV 
sigma = 0*mV ## noise amplitude
tau_n = 1*ms

eqs = Equations('''
dv/dt = (-(v-El)-ge*(v-Ee)-gi*(v-Ei)+sigma*xi*tau_n**(0.5))/taum : volt
dge/dt = -ge/taue : 1
dgi/dt = -gi/taui : 1
''')
## v : menbrane potential
## w : spiking threshold
## ge : excitatory gate
## gi : inhibitory gate


print "SC connections initialisation..."
SC = NeuronGroup(N= Nb, model = eqs, threshold = Vt, reset= Vr, refractory =  1.5*ms)
#SC = NeuronGroup(N= Nb, model = eqs, threshold = Vt, reset= "V=Vr; ge=0; gi=0", refractory =  1.5*ms)

CLe = Connection(SC, SC, 'ge')
CLi = Connection(SC, SC, 'gi')

DoG = calcul.connection_map4(width, 5.0, 5.0, 1.2, 6.0) 
if PLOT_MH:
    Z = DoG[Nb/2+width/2].reshape((width, width))
    aff.surface3D(Z, width)

map_inhib = DoG.copy()
map_inhib[map_inhib>0] = 0
map_inhib = map_inhib
map_excit = DoG.copy()
map_excit[map_excit<0] = 0

CLe.W = SparseConstructionMatrix(We*map_excit)
CLi.W = SparseConstructionMatrix(-Wi*map_inhib)

conn_time = time.time()-start
a = time.gmtime(conn_time)
print time.strftime("\ntask time for connection calculation: %H:%M:%S",a )   
    
result_firing = np.zeros(((end-beginning)/step+1, 1))   
table_line = 0

date_simu =  time.strftime("%d%b%Y-%H.%M.%S", time.gmtime())
path = "./"+date_simu+choose+str(beginning)+"-"+str(end)
os.mkdir(path)
os.chdir(path)
shutil.copy2(os.path.realpath(__file__), "copy_source.py")
os.mkdir("./state-var")
os.mkdir("./maps")
os.mkdir("./img")

for i in xrange(beginning,end,step): 
    print "Preparation simulation taille", str(i)
    reinit_default_clock()
    reinit(True)
    SC.v = El# + (rand(Nb) * (El - Vr))#/2    

    G1.freq = calcul.createBurst(400*Hz, 25, 80, 200)
        
    elec_pos2 = (50-i/2, 50-i/2) ## electrode sur la barre
    elec_pos = (50, 50)
    w_elec1 = 4000*mV    ## 1000 mV
    epi = calcul.c2DTo1D(elec_pos, width) ## center of the bar
    epi2 = calcul.c2DTo1D(elec_pos2, width) ## bar extremities
    if choose == "line":
        pattern1 = calcul.creerRectangle(width, elec_pos[0]-i/2, elec_pos[1], i, 1)
    elif choose == "rectangle":
        pattern1 = calcul.creerRectangle(width, elec_pos[0]-i/2, elec_pos[1]-i/2, i, i)
    elif choose == "circle":
        pattern1 = calcul.createCircleStim(width, elec_pos[0], elec_pos[1], i/2, i/2)
    elif choose == "gaussian":
        pattern1 = calcul.createGaussianStim(elec_pos[0], elec_pos[1], i/2, width)
    else:
        print "Error in the stimuli's shape choice..."
        sys.exit()
        
    
    neuronsbarre = range((52*width+1),(52*width+101))
    
    C1 = Connection(G1, SC, 'ge')
    C1.W = SparseConstructionMatrix(pattern1 * w_elec1)
        
    Mv = StateMonitor(SC, 'v', record=True)
    Ms = SpikeMonitor(SC)
    #Ms.spikes = [(i, t) for (i, t) in Ms.spikes if i in neuronsbarre]
    MsG1 = SpikeMonitor(G1)
    Cs = SpikeCounter(SC)
    Mge1 = StateMonitor(SC, 'ge', record=neuronsbarre[60])
    Mge2 = StateMonitor(SC, 'gi', record=neuronsbarre[60])
    ## --> run <-- ##
    print "Simulation started..."
    run(sim_time * second, threads = 4, report = 'text', report_period = 0.5)

    comp_time = time.time()-start-conn_time
    a = time.gmtime(comp_time)
    print time.strftime("task time for simulation: %H:%M:%S",a )
    
    print "Save pictures..."
    fig = figure(figsize=(24,12.3), dpi=80)
    plt.subplot(131)
    aff.stimPlot([w_elec1*pattern1], width)
    plt.subplot(132)
    aff.avStatePlot(Mv, sim_time, width)
    plt.subplot(133)
    aff.avSpikeplot(Cs, sim_time, width)
    plt.savefig("img/picture"+str(i)+".png",format='png')
    plt.close()
    print "Save data..."
    bin = 0.0001
    nline = sim_time/bin
    np.savetxt("state-var/pm"+str(i)+".out", Mv.values, delimiter = ',') 
    np.savetxt("state-var/spikes_SC"+str(i)+".out", Ms.it, delimiter = ',')
    map_stim = pattern1.reshape(width, width) * w_elec1
    np.savetxt("maps/map-stim"+str(i)+".out", map_stim, delimiter = ',')
    map_mean_pm = np.array(Mv.mean)
    np.savetxt("maps/mean-map-pm"+str(i)+".out", map_mean_pm.reshape(width, width), delimiter = ',')
    map_mean_fr = np.array(Cs.count/sim_time)
    np.savetxt("maps/mean-map-fr"+str(i)+".out", map_mean_fr.reshape(width, width), delimiter = ',')
    table_line += 1
    #show()
np.savetxt("state-var/spikes_input.out", MsG1.it, delimiter = ',')
print "End of the Program."