#Must run this program with the working directory that contains the Seed.csv file

seed<-read.csv("Seed.csv")

for (j in 1:length(seed[,1])){
  set.seed(seed$seed[j])

#Initial values for the modeled comparison PEM:
x1.act <- -.944957
x2.act <- .00301801
x3.act <- 7.401e-5
x4.act <- -1.88e-4
lamda.act <- 23
Rc.act <- .0001
B.act <- .02914489

param.act <- c(x1.act, x2.act, x3.act, x4.act, lamda.act, Rc.act, B.act)
amps <- c(1.1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)

#Ranges for each parameter:
x1.range <- c(-1.19969, -0.8532)
x2.range <- c(0.001, 0.005)
x3.range <- c(3.6e-5, 9.8e-5)
x4.range <- c(-2.6e-4, -9.54e-5)
lamda.range <- c(10, 24)
Rc.range <- c(0.0001, 0.0008)
B.range <- c(0.0136, 0.5)


#Single Cell Voltage Output Calculations:
Vout <- function (param) {
  v <- NULL

  for (i in amps) {
  Ph20 <- 10^(
    (2.95e-2*(353.15-273.15)) 
    - 9.19e-5*(353.15 - 273.15)^2 
    +1.44e-7*(353.15 - 273.15)^3
    -2.18)
  
  Ph2 <- 0.5*Ph20*((1/(Ph20/2.96077)^(1.635*(i/27)/353.15^1.334)) - 1)
  Po2 <- Ph20*((1/(Ph20/4.934616)^(4.192*(i/27)/353.15^1.334)) - 1)
  
  En <- 1.229 - 0.85*0.001*(353.15-298.15) + 4.3085e-5*353.15*log(Ph2*sqrt(Po2))
  
  co <- Po2/((5.08e6)^-(498/353.15))
  Nact <- -(param[1] + 353.15*param[2] + log(co)*353.15*param[3] + log(i)*353.15*param[4])
  
  Nconc <- -param[7]*log(1-((i/27)/0.860))
  
  pm <- (181.6*(1 + 0.03*(i/27) + 0.062*(353.15/303)*(i/27)^2.5))/(param[5] - 0.634 - 3*(i/27))^(4.18*((353.15-303)/353.15))
  Rm <- (pm*0.0127)/27
  Nohm <- i*(Rm + param[6])
  
  v <- append(v, (En - Nact - Nohm - Nconc))
  }
  return(v)
}

V.act <- Vout(param.act)

#SSE Evaluation (Fit Function)
SSE <- function (V.model) {
  temp2 <- data.frame(V.act, V.model)
  temp2$diff <- abs((temp2[,1]*24) - (temp2[,2]*24))  
  fit <- (sum(temp2$diff))^2
  return (fit)
}

#Setting the N of population:
n.pop <- 100

#Population of Randomly Generated numbers, within the parameters of each attribute respectively
pop <- data.frame("x1" = runif(n.pop, x1.range[1], x1.range[2]), "x2" = runif(n.pop, x2.range[1], x2.range[2]), "x3" = runif(n.pop, x3.range[1], x3.range[2]), 
        "x4" = runif(n.pop, x4.range[1], x4.range[2]), "Lamda" = runif(n.pop, lamda.range[1], lamda.range[2]), "Rc" = runif(n.pop, Rc.range[1], Rc.range[2]),
        "B" = runif(n.pop, B.range[1], B.range[2]))

pop

#Applying a fit SSE Score for each member of population
for (i in 1:n.pop){
  temp <- data.frame(Vout(c(pop[i,1], pop[i,2], pop[i,3], pop[i,4], pop[i,5], pop[i,6], pop[i,7])))
  pop$Fit[i] <- SSE(temp)
}

pop <- pop[order(pop$Fit),]

#Finding the best based on least amount of error
best.global <- pop[1, ]

count <- 1

###Loop to iniate stopping conditions here######
#while(best.global$Fit > .1){
for (i in 1:5000){

  if(best.global$Fit > pop[1,]$Fit) {
    best.global <- pop[1,]
    print(paste("New Best Found: Count:", count, best.global$Fit, sep = "     "))
  }
  
  #Selection--randomized between 1% and 90%
  selection.percent <- runif(1, .01, .9)
  cutoff <- floor(selection.percent*length(pop[,1])) 
  pop.new <- data.frame(pop[1:cutoff,])
  
  
  ##Crossover##
  #Number of replacement population needed:
  needed <- ceiling(n.pop*(1 - selection.percent))
  
  
  #This will split the needed pop refills into those that are from crossed parents, and those that get generated at random
  split <- floor(runif(1, 1, needed))
  if (split%%2 == 1) {split <- split + 1}
  
  
  parents <- sample(1:length(pop.new[,1]), size = split*2, replace = TRUE)
  crossPoint <- sample(1:7, size = split, replace= TRUE)
  
  start <- length(pop.new[,1])
  
  for (i in 1:split) {
    pop.new[start +i, ] <- (cbind(pop.new[parents[i*2-1], 1:crossPoint[i]], pop.new[parents[i*2], ((crossPoint[i] + 1):length(pop.new[1, ]))]))
  }
  
  #Now, create the rest by random generation:
  start <- length(pop.new[,1])
  
  for (i in 1:(needed-split)){
  pop.new[start + i, ] <- cbind(runif((needed-split), x1.range[1], x1.range[2]), runif((needed-split), x2.range[1], x2.range[2]), runif((needed-split), x3.range[1], x3.range[2]), 
                                runif((needed-split), x4.range[1], x4.range[2]), runif((needed-split), lamda.range[1], lamda.range[2]), runif((needed-split), Rc.range[1], Rc.range[2]),
                                runif((needed-split), B.range[1], B.range[2]), pop.new[cutoff,]$Fit)
  }
  
  ##Mutation--sigmoid function decreasing to min mutation
  mut_C <- 0.1
  mut_L <- 0.075
  count.max <- 5000
  
  
  for (i in 1:length(pop.new)){
    for (j in 1: length(pop.new[,1])) {
      mutation <- rbinom(1, size = 1, prob = mut_L + ((mut_C-mut_L)/(1+((2*count)/count.max)^2)))
      if (mutation == 1) {
        if(i == 1){
          pop.new[j,i] <- runif(1, x1.range[1], x1.range[2])
        }
        else if (i == 2) {
          pop.new[j,i] <- runif(1, x2.range[1], x2.range[2])
        }
        else if (i == 3) {
          pop.new[j,i] <- runif(1, x3.range[1], x3.range[2])
        }
        else if (i == 4) {
          pop.new[j,i] <- runif(1, x4.range[1], x4.range[2])
        }
        else if (i == 5) {
          pop.new[j,i] <- runif(1, lamda.range[1], lamda.range[2])
        }
        else if (i == 6) {
          pop.new[j,i] <- runif(1, Rc.range[1], Rc.range[2])
        }
        else if (i == 7) {
          pop.new[j,i] <- runif(1, B.range[1], B.range[2])
        }
      }
    }
  }
  
  #Recalculate and re-order before going back through the loop:
  for (i in 1:length(pop.new[,1])){
    temp3 <- data.frame(Vout(c(pop.new[i,1], pop.new[i,2], pop.new[i,3], pop.new[i,4], pop.new[i,5], pop.new[i,6], pop.new[i,7])))
    pop.new$Fit[i] <- SSE(temp3)
  }
  pop.new <- pop.new[order(pop.new$Fit),]
  
  pop <- pop.new
  count <- count + 1
  if (count == 5000) {break}
}
###End Loop here ####

print ("Optimal parameters found: ")
print(best.global)

print(c('Total iterations: ', count))

#table <- merge(V.act, Voltage(c(best.global[,1], best.global[,2], best.global[,3], best.global[,4], 
#                                best.global[,5], best.global[,6], best.global[,7])), by = "Amps")
#table
#plot(table)

output <- c("Fit: ", best.global$Fit, "Iterations: ", count)

##This was used to output the results for each seed run.  Commented out to prevent inadvertant information overwriting:
#write.table(output, file = "GA_Output_Sigmoid_.08_5000_Iter.csv", append = TRUE)
}



