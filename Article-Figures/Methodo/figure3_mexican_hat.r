setwd("D:\\Documents\\Figures\\Methodo")

DOG <- function(x,y, sx, sy, K, inh){ ## similar than connectionmap4() in outil.calcul.py
    A1 <- 1.0 + inh
    A2 <- inh
    M1 <- A1*exp(-(  ((x)**2)/(2*sx**2) + ((y)**2)/(2*sy**2) ) )
    M2 <- -A2*exp(-(  ((x)**2)/(2*(K*sx)**2) + ((y)**2)/(2*(K*sy)**2) ) )
    return(M1+M2)
}

## chapeau used for auto-inhibition:
x <- seq(0,30,0.1)
SA1 <- DOG(x,0,5.0,5.0,1.2,6.0)
SA2 <- DOG(x,0,5.0,5.0,2.0,1.42868)
SA3 <- DOG(x,0,5.0,5.0,1.2,8.0)
pdf("pepin.pdf")
plot(x, SA3*200, col="green", type='l', lwd=3, lty=2)
points(x, SA2*200, col="red", type='l', lwd=3, lty=5)
points(x, SA1*200, type='l', col="black", lwd=3)
dev.off()