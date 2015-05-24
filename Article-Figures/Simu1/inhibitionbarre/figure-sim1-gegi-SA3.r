main_folder <- "D:/Documents/Figures/Simu1/inhibitionbarre"
setwd(main_folder)

## work on SA3
MH_folder <- paste0(main_folder, "/cste-10-3-10-sy50sx50k1.2inh8.0/")

print(paste("Working on:", MH_folder))
sizes_plotted <- c(2, 28, 36)
tracked_neurons <- c(50, 59, 61)
stimulations_type <- c("transient", "constant")
states_var <- c("pm","ge", "gi")
list_matplot <- list() ## list of matrices, each matrice contains the curves for one plot.
## we want a 2x3 plots
nplot = 1
for (stim in stimulations_type) {
	print(paste("Stimulation type: ", stim))
	current_folder <- paste0(MH_folder,"gegi-",stim)
	for (i in (1:length(sizes_plotted))) {
		m <- c() ## sack of simp0licity, we now the number of rows of our files
		for (st in states_var) {
			f <- paste0(current_folder,"/" ,st, sizes_plotted[i], ".out")
			print(paste("Opening", basename(f), "..."))
			stat_var <- as.matrix(read.table(f, sep=','))
			if(st == 'pm'){
				stat_var <- (stat_var[,tracked_neurons[i]]+0.07) * 3000
			} else if (st=='ge') {
				#pm_last <- c(0, m[,1][-length(m)])/3-70
				stat_var <- stat_var[,tracked_neurons[i]]*65
			} else if (st=='gi') {
				#pm_last <- c(0, m[,1][-length(m)])/3-70
				#print(paste("moyenne: ",mean(tail(stat_var[,tracked_neurons[i]], 1000))))
				stat_var <- stat_var[,tracked_neurons[i]]*15 #+ 5
			} 			
			m <- cbind(m, stat_var)	
			print(dim(m))
		}
		#stat_var <- (70 + (0*m[,1] + (-80)*m[,2] + (-70))/(m[,1] + m[,2] + 1))*3
		#stat_var <- (70 + (0*m[,1] + (-80)*m[,2] + (-70))/(m[,1] + m[,2] + 1))*3
		#m <- cbind(m, stat_var)
		list_matplot[[nplot]] <- m
		nplot <- nplot + 1
	}
}

## margin setting:  c(bottom, left, top, right) == c(5, 4, 4, 2) + 0.1.
pdf("figure-sim1-gegi-SA3.pdf", width = 7.0, height=7.0/1.8, pointsize = 1/300)
## we know how the time is coded in our files:
time_line <- seq(0.1, 200, 0.1)
par(mfrow=c(2,3), cex=1.5, cex.main=1.0, cex.lab = 0.5, mgp=c(3,0.5,0))
par(oma = c(4.1,   1,   0,3.0))
par(mar=  c(0  , 1.6, 2.3,  0))
nplot = 1
col_order = gray.colors(3, start=0.0, end=0.8)[c(3,1,2)]#c(3,2,4)
typ_lines = c(1,1,1)
lwd_lines = c(0.5,1,1)
xvalues = 'n'
for (mat in list_matplot) {
	if (nplot == 1 | nplot == 4) {
		yvalues = 'l'
	} else {
		yvalues = 'n'
	}
	if (nplot > 3) {xvalues <- 'l'}
	print(dim(mat))
	print(length(time_line))
	
	matplot(time_line, mat, ## convert in mV
	xlab = NA, xaxt = xvalues, 
	ylab = NA, yaxt = yvalues, ylim = c(0, 80),
	col=col_order, type='l', lty = typ_lines, lwd = lwd_lines) 
	
	abline(h=(-51), lty=3)
	if(nplot == 3 | nplot == 6){Axis(side=4, at=seq(0,60, length.out=3),labels = c(-70,-60,-50), las = 0)}
	if(nplot != 1 & nplot != 4) { Axis(side=2, labels=FALSE)}
	if(nplot < 4) { Axis(side=1, labels=FALSE)}
	## legend creation:
	if (nplot <4){
	title_plot <- paste0("Stimulus size = ", sizes_plotted[nplot], 
	"\n", "Tracked neuron = 52x", tracked_neurons[nplot])
	#mtext(title_plot, at=c(0+nplot, 0),outer=TRUE, side =3)
	title(title_plot)
	}
	legend_text <- c("V", "ge", "gi")
	# for (st in states_var) {
		# legend_text <- c(legend_text, st)
	# }
	legend("topright", legend_text, col= col_order, lty = typ_lines, lwd = lwd_lines, cex=0.8, bg="white")
	nplot = nplot + 1
}
mtext("Channels Opening (Normalized in mV)", side = 2 , outer = T, cex=1.5)
mtext("Membrane Potential (mV)", side = 4 , outer = T, cex=1.5, line=1.6)
mtext("Simulation time (ms)", side=1, outer = T, cex=1.5, line=2.5)
dev.off()




