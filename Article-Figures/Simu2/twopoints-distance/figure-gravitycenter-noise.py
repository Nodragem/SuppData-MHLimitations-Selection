import numpy as np
import matplotlib.pyplot as plt
import os
import re
from outil import calcul

list_color = ['#FF6699', '#FF3300', '#99FF33', '#0099FF', '#CC0000', '#001F33', '#47008F', '#339933']
list_markers = np.array(["o", "v", "*"])
cm = plt.get_cmap('gist_rainbow')
print cm(0.5)
exp_folder = "./exp1-noise"
list_cond =  os.walk(exp_folder).next()[1]
print list_cond
cond_example = "8-8-1-6-T4-D4-fT500-fD500-c1"
img_example = []
tab_example = None
graphs = []

for cond in list_cond:
    folder = exp_folder + "/" + cond + "/"
    print folder
    tab_gc = np.loadtxt(folder+"distance-40-gravitycenter.txt", delimiter=",")
    target_files = np.array([file for file in os.listdir(folder) if re.match("^distance-\d+-av-spikes", file)])
    distances = np.array([int(re.findall("\d+", d)[0]) for d in target_files])
    target_files = target_files[np.argsort(distances)]
    distances = distances[np.argsort(distances)]
    #print target_files
    ref_oneblob = None
    col_newgc = []
    col_number_blob =[]
    for i in xrange(len(target_files)):
        m = np.load(folder+target_files[i])
        if cond_example == cond:
            img_example.append(m)
        d = distances[i]
        if i == 0:
            ref_oneblob = np.sum(m)
        if np.sum(m) < (0.6*ref_oneblob): ## zero blob
            col_newgc.append((np.nan, np.nan))
            col_number_blob.append(0)
        elif np.sum(m) < (1.2*ref_oneblob): ## one blobs
            gc = calcul.calculCentreGrav(m)
            col_newgc.append(gc)
            col_number_blob.append(1)
        else: ## two blobs
            xmidline = (tab_gc[i, 1] + tab_gc[i, 2])/ 2 ## position inbetween the target and the  distractor
            ## we will ignore  anything from the distractor side (erase the value on the left side)
            mcopy = m.copy()
            mcopy[:, 0:int(xmidline)] = 0 ## weirdly the x and y are inversed
            gc = calcul.calculCentreGrav(mcopy)
            col_newgc.append(gc)
            col_number_blob.append(2)
            ################################
            # fig_debug = plt.figure()
            # ax_dg1 = plt.subplot(211)
            # ax_dg1.imshow(m)
            # ax_dg2 = plt.subplot(212)
            # ax_dg2.imshow(mcopy)
            # ax_dg2.scatter(gc[1], gc[0], color="green")
            # ax_dg2.scatter(tab_gc[i, 2], 49.5, color="white")
            # plt.show()

    col_number_blob = np.array(col_number_blob).reshape((20,1))
    col_newgc = np.array(col_newgc).reshape((20,2))
    distances = distances.reshape((20,1))
    print col_newgc.shape, distances.shape, tab_gc.shape, col_number_blob.shape
    tab_gc = np.hstack((distances, tab_gc, col_newgc, col_number_blob))
    if cond == cond_example:
        tab_example = tab_gc
    np.savetxt(folder+"readytoGraph.txt", tab_gc, delimiter=",")
    graphs.append(tab_gc)

choice = [5, 9, 12]

fig = plt.figure(figsize=(7.4, 6.0), dpi = 80)
#plt.figtext(0.05, 1- 0.13, 'A', fontsize=30)
#plt.figtext(0.05, 1- 0.33-0.05, 'B', fontsize=30)
#plt.figtext(0.05, 1- 0.66-0.03, 'C', fontsize=30)
ax = fig.add_subplot(111)
for i, g in enumerate(graphs):
    if "c0" in list_cond[i]:
        continue
    print g
    deviation = g[:,-2]-g[:,3] + 0.5 ## let's put a minus one
    ax.hlines(0, 0, 45, linestyles="solid", colors="k", zorder=0, linewidth=0.5)
    ax.hlines(-5, 0, 45, linestyles="dashed", colors="gray", zorder=0)
    ax.hlines(5, 0, 45, linestyles="dashed", colors="gray", zorder=0)
    ax.plot(g[:,0], deviation , color=list_color[i], zorder=1)
    ax.scatter(g[:,0], deviation,s=50, lw = 2,
               edgecolor=list_color[i],
               facecolor= np.array(["gray", list_color[i], "white"])[g[g[:,-1]!=0,-1].astype("int")],
               alpha = 1, zorder=2)
    ax.set_xlim(0,45)
    ax.set_ylim(-26,12)
    # if (list_cond[i] == cond_example):
    #     for number, k in enumerate(choice):
    #         #ax.arrow(g[k-1,0]-3, deviation[k-1]-3, 2, 2, head_width=0.5, head_length=0.5, fc=list_color[i], ec=list_color[i])
    #         ax.annotate(str(number+1), xy=(g[k-1,0]-0.5, deviation[k-1]-0.5), xytext=(g[k-1,0]-3, deviation[k-1]-3),
    #         arrowprops=dict(headwidth=10,frac=0.40, shrink=0.03, fc=list_color[i], ec=list_color[i])            )
    # if (list_cond[i] == "8-8-1-6-T4-D2-fT500-fD500-c1"):
    #     for number, k in enumerate((36/2,38/2,40/2)):
    #         ##ax.arrow(g[k-1,0]+3, deviation[k-1]-3, -2, 2, head_width=0.5, head_length=0.5, fc=list_color[i], ec=list_color[i])
    #         ax.annotate("  ", xy=(g[k-1,0]+0.5, deviation[k-1]-0.5), xytext=(g[k-1,0]+2.5, deviation[k-1]-2.5),
    #         arrowprops=dict(headwidth=10,frac=0.40, shrink=0.03, fc=list_color[i], ec=list_color[i])            )
    #         #ax.annotate(str(j+1), xy=(g[k-1,0]+0.5, deviation[k-1]-0.5), xytext=(g[k-1,0]+2, deviation[k-1]-3),
    #         #arrowprops=dict(headwidth=10,frac=0.40, shrink=0.03, fc=list_color[i], ec=list_color[i])            )

ax.set_ylabel("Deviation from the \n stimulation B (cells)")
ax.set_xlabel("Distance between stimulations (cells)")
legends = ["4000-1333mV", "4000-2000mV","4000-3500mV","4000-4000mV"]
ax.legend(legends, loc="lower left", fontsize=12)
plt.tight_layout()

# axes = []
# for j, i in enumerate(choice):
#     img = img_example[i-1]
#     m = tab_example[i-1, :]
#     ax = fig.add_subplot(331 + j)
#     im = ax.imshow(img, vmin= 0, vmax=600, interpolation="none")
#     plt.scatter(m[-2], m[-3], facecolor = "white")
#     plt.scatter(m[3], 49.5, facecolor = "magenta", lw=1, edgecolor="k")
#     plt.scatter(30, 49.5, facecolor = "green", lw=1, edgecolor="k")
#     axes.append(ax)
#     ax.set_yticks(np.arange(0,101,25))
#     ax.set_xticks(np.arange(0,101,25))
#     ax.set_yticklabels(("0","","50","","100"))
#     ax.set_xticklabels(("0","","50","","100"))
#     ax.text(0,100, "distance %d"%m[0], verticalalignment = "bottom", color = "white")
#     #ax.text(0,0, str(j+1), verticalalignment = "top", color = "white")
#     if j == 2:
#         ax.set_title("3.Repulsion")
#
#         cbar = plt.colorbar(im, ax= axes, ticks = np.arange(0, 601,300),fraction=0.05, pad = 0.01, shrink=0.6 )
#         cbar.set_label('Average Firing Rate (Hz)', rotation=270, labelpad=15.0)
#     if j > 0:
#         ax.set_yticklabels([])
#     if j == 1:
#         ax.set_title("Condition 4000-4000mV (Green Curve) \n 2.Suppression")
#         ax.set_xlabel("x-axis (cells)")
#     if j == 0:
#         ax.set_title("1.Attraction")
#         ax.set_ylabel("y-axis (cells)")
#
# path = [
#     ".\exp1\8-8-1-6-T4-D2-fT500-fD500-c1\distance-36-raw-mp.npy",
#     ".\exp1\8-8-1-6-T4-D2-fT500-fD500-c1\distance-38-raw-mp.npy",
#     ".\exp1\8-8-1-6-T4-D2-fT500-fD500-c1\distance-40-raw-mp.npy"
#        ]
#
# axes = []
# for i, f in enumerate(path):
#     m = np.mean(np.load(f), axis=2)
#     print m.shape
#     ax = fig.add_subplot(337 + i)
#     im = ax.imshow(m, vmin= -80, vmax=-50, interpolation="none")
#     #plt.scatter(m[-2], m[-3], facecolor = "white")
#     plt.scatter(30 + 36 + i*2, 49.5, facecolor = "magenta", lw=1, edgecolor="k")
#     plt.scatter(30, 49.5, facecolor = "green", lw=1, edgecolor="k")
#     axes.append(ax)
#     ax.set_yticks(np.arange(0,101,25))
#     ax.set_xticks(np.arange(0,101,25))
#     ax.set_yticklabels(("0","","50","","100"))
#     ax.set_xticklabels(("0","","50","","100"))
#     if i == 2:
#         ax.text(0,100, "distance 40", verticalalignment = "bottom", color = "white")
#         cbar = plt.colorbar(im, ax= axes, ticks = np.arange(-80, -49, 15),fraction=0.05, pad = 0.01, shrink=0.6 )
#         cbar.set_label('Average Memb. Potential (mV)', rotation=270, labelpad=15.0)
#     if i > 0:
#         ax.set_yticklabels([])
#     if i == 1:
#         ax.set_title("Condition 4000-2000mV (Blue Curve)" )
#         ax.text(0,100, "distance 38", verticalalignment = "bottom", color = "white")
#         ax.set_xlabel("x-axis (cells)")
#     if i == 0:
#         ax.text(0,100, "distance 36", verticalalignment = "bottom", color = "white")
#         ax.set_ylabel("y-axis (cells)")


plt.show()





