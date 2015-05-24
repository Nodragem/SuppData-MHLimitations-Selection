import numpy as np
import sys, os, time
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib

def toIndex(xy, width):
	return((xy[0]-1)*width + xy[1])


def selectNeurons(orientation, bias = 0, start=0, end=100, map_size=100):
	index_neuron = []
	if (orientation == "diag"):
		
		for i in xrange(start,end):
			index_neuron.append(toIndex((i+bias,i), map_size))
		
	elif (orientation == "vert"):
		
		for i in xrange(start,end):
			index_neuron.append(index_neuron, toIndex((50+bias,i), map_size))
		
	elif (orientation == "hori"):
		
		for i in xrange(start,end):
			index_neuron.append(index_neuron, toIndex((i,50+bias), map_size))
		
	return(np.array(index_neuron))

def rasterPlot(ax, mat, start, index_neurons, color='k'):
	
	for line, i in enumerate(index_neurons):
		spikes = mat[1, mat[0,] == i]*1000
		if spikes != []:
			ax.vlines(spikes, start+line + .6, start+line + 1.4, color=color, linewidth=0.5)
	ax.set_ylim(start+.5, start+len(index_neurons) + .5)
	ax.set_xlim(0, 200)
	ax.set_xticks([0,100,200])
	ax.set_yticks([25,50,75])
	#return ax

def draw_a_line(size_circ, size_rect, ax, corner=False):
	## open the data
	index_neurons = selectNeurons("diag", start=25, end=75)
	os.chdir("D:/Users/Administrator/workspace/Modelisation/projet1/inhibitionsquare/")
	mat_rect = np.loadtxt("./Rectangle/maps/mean-map-fr"+ str(size_rect) +".out", delimiter=",") #18
	mat_circ = np.loadtxt("./Circle/maps/mean-map-fr"+ str(size_circ) +".out", delimiter=",") #22
	os.chdir("D:/Users/Administrator/workspace/Modelisation/projet1/inhibitionsquare/Circle/state-var")
	spikes_circle = np.loadtxt("spikes_SC"+str(size_circ)+".out", delimiter=",")
	spikes_circle = spikes_circle[:, np.in1d(spikes_circle[0,], index_neurons)]
	os.chdir("D:/Users/Administrator/workspace/Modelisation/projet1/inhibitionsquare/Rectangle/state-var")
	spikes_rect = np.loadtxt("spikes_SC"+str(size_rect)+".out", delimiter=",")
	spikes_rect = spikes_rect[:, np.in1d(spikes_rect[0,], index_neurons)]
	# plot the data
	
	ax[0].imshow(mat_rect[25:75,25:75], interpolation="nearest", vmin=0, vmax=600, extent=(25,75,25,75))
	ax[0].text(74,70, 'size = '+str(size_rect),horizontalalignment = "right", color="white", fontsize=8)
	ax[0].set_xticks([25,50,75])
	if corner:
		ax[0].set_yticks([50,75])
	else:
		ax[0].set_yticks([25,50,75])
	ax[0].tick_params(colors='white', labelcolor='black')
	rasterPlot(ax[1],spikes_rect, 25, index_neurons)
	im = ax[2].imshow(mat_circ[25:75,25:75], interpolation="nearest", vmin=0, vmax=600, extent=(25,75,25,75), label='size = '+str(size_circ))
	ax[2].text(74,70, 'size = '+str(size_circ),horizontalalignment = "right", color="white", fontsize=8)
	ax[2].set_xticks([25,50,75])
	ax[2].set_yticks([25,50,75])
	ax[2].tick_params(colors='white', labelcolor='black')
	rasterPlot(ax[3], spikes_circle, 25, index_neurons)
	#plt.colorbar(im)
	return (im)


os.chdir("D:/Users/Administrator/workspace/Modelisation/projet1/inhibitionsquare/")
#pp = PdfPages('figure_spikes_rect_circ.pdf')
fig = plt.figure(figsize=(7,5.25), dpi=300)	
all_ax = []
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(3, 4)
for i in xrange(0,12):
	all_ax.append(plt.subplot(gs[i/4, i%4]))

print "Affichage..."
draw_a_line(8, 8, all_ax[0:4])
draw_a_line(20, 18, all_ax[4:8]) # 20 18
im = draw_a_line(30, 30, all_ax[8:12], corner = True)
for i, ax in enumerate(all_ax):
	if not (i == 0 or i==4 or i==8):
		ax.get_yaxis().set_ticklabels([])
	else: 
		ax.set_ylabel("Space (cells)")
	if i < 8:
		ax.get_xaxis().set_ticklabels([])
	elif i%2==0:
		ax.set_xlabel("Space (cells)")
	else:
		ax.set_xlabel("Time (ms)")
#plt.colorbar()
fig.subplots_adjust(top=0.87)
cbar_ax = fig.add_axes([0.35, 0.93, 0.3, 0.01])
plt.figtext(0.275, 0.89, "Square" ,horizontalalignment="center",fontsize=12)
plt.figtext(0.74, 0.89, "Circle" ,horizontalalignment="center",fontsize=12)
cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal", ticks=[0,300,600])
cbar.ax.set_xticklabels(['0', '300', '600 Hz'])
for tick in cbar.ax.get_xaxis().iter_ticks():
    tick[0].label2On = True
    tick[0].label1On = False

fig.subplots_adjust(right=0.95)
fig.subplots_adjust(left=0.07)
matplotlib.rcParams.update({'font.size': 8})

#pp.savefig(fig)
#pp.close()
plt.show()

#print np.mean(mat, 1)
print "Fin"
