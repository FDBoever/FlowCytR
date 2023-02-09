#https://www.cyclismo.org/tutorial/R/s4Classes.html

#generate S4 object
setClass('predator',slots=list(name='character',age='numeric', GPA='numeric'))

pr1 <- new('predator',name='x',age=21,GPA=4.5)

pr1@name

#modifying
pr1@GPA <- 2.7

slot(pr1,'name') <- 'pred1'


#---

#install package simcol
install.packages('simecol',dependencies = TRUE)



#=====


library(deSolve)

## -----------------------------------------------------------------------------
## Define R-function
## -----------------------------------------------------------------------------

LV <- function(t, y, parms) {
  with(as.list(c(y, parms)), {

    dP <- rG * P * (1 - P/K) - rI * P * C
    dC <- rI * P * C * AE - rM * C

    return(list(c(dP, dC), sum = C+P))
  })
}

## -----------------------------------------------------------------------------
## Define parameters and variables
## -----------------------------------------------------------------------------

parms <- c(rI = 0.2, rG = 1.0, rM = 0.2, AE = 0.5, K = 10)
yini <- c(P = 1, C = 2)
times <- seq(from = 0, to = 200, by = 1)

## -----------------------------------------------------------------------------
## Solve the ODEs
## -----------------------------------------------------------------------------

out <- ode(y = yini, times = times, func = LV, parms = parms)

## -----------------------------------------------------------------------------
## Plot the results
## -----------------------------------------------------------------------------

matplot(out[ ,1], out[ ,2:4], type = "l", xlab = "time", ylab = "Conc",
        main = "Lotka-Volterra", lwd = 2)
legend("topright", c("prey", "predator", "sum"), col = 1:3, lty = 1:3)
