main_folder <- "D:/Documents/Figures/Simu1/inhibitionbarre"
setwd(main_folder)


namefolders <- dir('.', pattern = "^cste")
#"cste-10-3-10-sx50-sy50-k20-inh[1.42484]"
#"cste-10-3-10-sy50sx50k12inh60" 
#"cste-10-3-10-sy50sx50k12inh80"
namefolders <- paste0(namefolders[c(2,1,3)], "/pm-HD")
#"cste-10-3-10-sy50sx50k12inh60/pm" 
#"cste-10-3-10-sx50-sy50-k20-inh[1.42484]/pm"
#"cste-10-3-10-sy50sx50k12inh80/pm"
sizes_plotted <- matrix(c(2,18,20,42, 0,
						  2,22,24,42, 0,
						  2,14,16,30,38), 3, 5, byrow=TRUE)
#     [,1] [,2] [,3] [,4]
#[1,]    2   20   42    0
#[2,]    2   24   42    0
#[3,]    2   16   30   38
mat_speed_ini <- c()
mat_speed_second <- c()
list_matplot <- list() ## list of matrices, each matrice contains the curves for one plot.
list_parameters <- list()
#list_sizes <- list()
condition = 1

for (folder in namefolders) {
	mat_speed_ini <- cbind(mat_speed_ini, read.table(paste0('./',folder,"/speed_ini.txt"), sep=" ")[,2])
	mat_speed_second <- cbind(mat_speed_second, read.table(paste0('./',folder,"/speed_second.txt"), sep=" ")[,2])
	curves_to_plot <- c()
	print(paste("Working on:", folder))
	namefiles <- list.files(paste0('./',folder),pattern = "^pm\\d+\\.out")
	print(namefiles)
	sizefiles <- as.numeric(gsub("\\D", "", namefiles))
	for (size in sizes_plotted[condition,]) {
			if (size > 1){
				f <- paste(folder, namefiles[sizefiles==size], sep="/")
				m <- as.matrix(read.table(f, sep=','))
				average_space <- apply(m, 1, max)*1000 ## max of the line, converted in mV
				## suppress the curve's part after the spiking threshold is reached:
				reach_threshold <- which(average_space > -51.0 )[1]
				if (!is.na(reach_threshold)) { 
				print(paste("size:", size, "Delay to threshold:", reach_threshold))
				average_space[reach_threshold:length(average_space)] <- NA}
				## we assumed that all the files have the same number of rows 
				##(in a true program, you should not assume and should test it)
				curves_to_plot <- rbind(curves_to_plot, average_space) 
			}
	}
	list_matplot[[condition]] = curves_to_plot
	numbers_in_folder <- as.numeric(unlist(strsplit(folder, "[^0-9.]")))
	numbers_in_folder <- numbers_in_folder[!is.na(numbers_in_folder)]
	list_parameters[[condition]] <- tail(numbers_in_folder,2)
	condition = condition + 1
}

## margin setting:  c(bottom, left, top, right) == c(5, 4, 4, 2) + 0.1.
pdf("figure-sim1-essai2.pdf", width = 7.0, height=7.0/1.5, pointsize = 1/300)
## we know how the time is coded in our files:
time_line <- seq(0.1, 200, 0.1)
#par(mfrow=c(1,length(list_matplot)))
layout(matrix(c(1,2,3,4,4,5), 2, 3, byrow = TRUE))
par(cex=1.5)
par(oma = c(2.1,2,0,2))
condition = 1
col_order = 1:dim(sizes_plotted)[2]
for (mat in list_matplot) {
	if (condition == 1) {
		par(mar=c(2.1, 2.1, 2.1, 0))
		values = 'l'
	} else {
		par(mar=c(2.1, 2.1, 2.1, 0))
		values = 'n'
	}
	matplot(time_line, t(mat), ## convert in mV
	xlab = NA, 
	ylab = NA,
	yaxt = values,
	col=col_order, type='l', lty = 1, lwd = 1.0) 
	abline(h=(-51), lty=3)
	if(condition != 1) { Axis(side=2, labels=FALSE)}
	## legend creation:
	parameters <- list_parameters[[condition]]
	title_plot <- substitute(paste("SA", co, ": ",K == k, " , ", beta == b), list(co=condition, k=parameters[1], b=parameters[2]))
	title(title_plot)
	legend_text <- c()
	for (size in sizes_plotted[condition,]) {
		if (size > 1) {legend_text <- c(legend_text, paste("size = ", size))}
	}
	legend("topright", legend_text, col= col_order, lty=1, lwd=1.0, cex=0.8, bg="white")
	condition = condition + 1
}
mtext("Membrane potential speed (mV/s)     Max. of membrane potential (mV)", side = 2 , outer = T, cex=1.5)
mtext("Simulation time (ms)", side=1, outer = T, cex=1.5, line = -14.5)
mtext("Stimulus sizes (cells)", side=1, outer = T, cex=1.5)

matplot(seq(2,42,2), mat_speed_ini, col=c("black", "red", "green"), lty=1, type='l', lwd=1.5)
legend("topright","Initial Speed", cex=0.8)
text(30, mat_speed_ini[15,2]+0.05, "SA2", cex=0.8)
text(30, mat_speed_ini[15,1]+0.05, "SA1" , cex=0.8)
text(30, mat_speed_ini[15,3]+0.05, "SA3", cex=0.8)
##14 18 22
points(seq(2,18,2), mat_speed_ini[1:9,1])
points(seq(2,22,2), mat_speed_ini[1:11,2])
points(seq(2,14,2), mat_speed_ini[1:7,3])
k <- length(mat_speed_second[,1])
mat_speed_second[1:10,1] <- NA
mat_speed_second[1:12,2] <- NA
matplot(seq(16,42,2),
mat_speed_second[8:k,] , ylim=c(-0.05, 0.25), col=c("black", "red","green"), lty=1, type='l', lwd=1.5)

legend("topright", "Speed after the inhibition wave", cex=0.8)
points(seq(30, 36, 2), mat_speed_second[15:18,3])

dev.off()
