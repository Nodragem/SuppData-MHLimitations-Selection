# Short Mexican Hat Limitations for Target Selection -- Code Source.

#### Full Title:
Limitations of Short Range Mexican Hat Connection for Driving Target Selection in a 2D Neural Field: Activity Suppression and Deviation from Input Stimuli.

#### Authors:
Geoffrey Mégardon (1,2); Christophe Tandonnet (2,3);  Petroc Sumner (1), Alain Guillaume (2,4)

1. School of Psychology, Cardiff University, Tower Building, 70 Park Place, Cardiff, UK.
2. Laboratoire de Neurobiologie de la Cognition, UMR CNRS 6155, Aix-Marseille Universitée, 3 place Victor-Hugo, 13331, Marseille, France.
3. Faculté de Psychologie et des Sciences de l'Education, Université de Genève, 40 bd du Pont d'Arve, CH-1205 Genève, Suisse.
4. Department of Psychology, New York University, 6 Washington Place, New York, NY, 10003, USA.

#### Abstract:
Dynamic Neural Field model (DNF) often uses a kernel of connection with short range excitation and long range inhibition. This organization has been suggested as a model for central nervous system structures or for artificial systems involved in winner-take-all processes such as perceptual decision or target selection. The superior colliculus (SC) is seen as a good example of neural structure containing such a DNF for target selection in the context of saccadic eye movements.

Recent results suggested that the superficial layers of the SC (SCs) is optimized for target selection, however it also exhibits a relatively short range of inhibition. The aims of the present study were 1) to further examine the properties of DNF when using relatively short range of inhibition, and 2) to highlight the consequences of those properties when applied to target selection.

We created a 2D spiking DNF with connections similar to those found in the SCs. It could be stimulated with input of different shapes and positions; consequential self-maintained cluster of firing neurons was considered the selected targets.

In a first protocol, we tested the target selection for a single stimulus with varying size and shape. For all tested shapes, a range of sizes led to the suppression of any activity on the network and hence to no target selection, or to a delayed selection of multiple loci. In a second protocol, we tested the target selection with two stimuli separated by a varying distance. Again no target selection and multi-selection could occur. In addition, attraction and repulsion influences between the two target candidates were observed.

Those properties and their functional consequences are discussed respectively, at the neurophysiological and behavioral level. Finally, as we used a spiking neuron model instead of the traditional population rate model, our work shed some light on properties they share.

## What do you need to run the programs:
*I recommend to install a scientific distribution of python as Anaconda (64-bit).*

You need python 2.7 - 64bit installed with the libraries:
- Matplotlib (64bit)
- Numpy (64bit)
- PIL (64bit)
- BRIAN 1.4.1

Note that BRIAN won't be provided by Anaconda. You will find it here:
http://brian.readthedocs.org/en/latest/installation.html

Note that the 64-bit versions of the above libraries are needed only for the Simulation1 and Simulation2 sources.
Indeed, those programs need 4-6 GB of RAM to run: only 64-bit programs can use more than 3 GB of RAM on your PC.
The figures and the SCtoVisual transformation work with the 32-bit versions.

Be careful, 32-bit libraries and python 32-bit do not work with 64-bit libraries.
Find 64-bits libraries for python here:
http://www.lfd.uci.edu/~gohlke/pythonlibs/
