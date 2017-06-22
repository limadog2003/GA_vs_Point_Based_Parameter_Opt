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


#Mutation function per heuristic guidelines:
mutate <- function(oldValue) {
    mutators <- sample(1:7, size = runif(1, 1, 7), replace = FALSE)
    
    for (i in mutators) {
      
      if (i == 1) {range <- x1.range}
      else if (i == 2) {range <- x2.range}
      else if (i == 3) {range <- x3.range}
      else if (i == 4) {range <- x4.range}
      else if (i == 5) {range <- lamda.range}
      else if (i == 6) {range <- Rc.range}
      else if (i == 7) {range <- B.range}
      
      x <- runif(1, 0, 1)
      v <- x*(range[1] - range[2])*runif(1,0,1)
      
      if (runif(1, 0, 1) >= .5) { v <- -v}
      
      new <- (oldValue)[i] + v
      
      lo <- range[1]
      if (new < lo) {
        if ((oldValue)[i] == lo) {
          new <- lo
        }
        else {
          new <- lo + runif(1, 0, 1)*((oldValue)[i]-lo)
        }
      }
      
      hi <- range[2]
      if (new > hi){
        if ((oldValue)[i] == hi){
          new <- hi
        }
        else{
          new <- hi - runif(1, 0, 1)*(hi - (oldValue)[i])
        }
      }
      oldValue[i] <- new
    }
    
    newValue <- oldValue
    
    return (newValue)
}


#Calculations for the average d:
n.pop <- 50

pop <- data.frame("x1" = runif(n.pop, x1.range[1], x1.range[2]), "x2" = runif(n.pop, x2.range[1], x2.range[2]), "x3" = runif(n.pop, x3.range[1], x3.range[2]), 
                  "x4" = runif(n.pop, x4.range[1], x4.range[2]), "Lamda" = runif(n.pop, lamda.range[1], lamda.range[2]), "Rc" = runif(n.pop, Rc.range[1], Rc.range[2]),
                  "B" = runif(n.pop, B.range[1], B.range[2]))

pop.mut <- as.data.frame(pop[0,])
for(i in 1:n.pop){
  pop.mut[i,] <- mutate(c(pop[i,1], pop[i,2], pop[i,3], pop[i,4], pop[i,5], pop[i,6], pop[i,7]))
  }

for (i in 1:n.pop){
  pop$Fit[i] <- SSE(Vout(c(pop[i,1], pop[i,2], pop[i,3], pop[i,4], pop[i,5], pop[i,6], pop[i,7])))
  
  pop.mut$Fit[i] <- SSE(Vout(c(pop.mut[i,1], pop.mut[i,2], pop.mut[i,3], pop.mut[i,4], pop.mut[i,5], pop.mut[i,6], pop.mut[i,7])))
}


d.avg <- abs(sum((pop.mut[,8] - pop[,8]))/n.pop)

pop <- pop[order(pop$Fit),]
pop.mut <- pop.mut[order(pop.mut$Fit), ]

if(pop$Fit[1] < pop.mut$Fit[1]){
  best <- pop[1, ]
  }
if (pop$Fit[1] >= pop.mut$Fit[1]){
  best <- pop.mut[1, ]
  }

previous <- best
print(best)
count <- 1

###Looping conditions###
#while(count < 500000){
while (best$Fit > 0.1){
  maxCount <- 1000000
  current <- mutate(previous[1:7])
  current$Fit <- SSE(Vout(c(current[,1], current[,2], current[,3], current[,4], current[,5], current[,6], current[,7])))
  
  if (current$Fit < best$Fit){
    best <- current
    print("Found new best")
    print(paste(count, best$Fit, sep = "     "))
    next
  }
  
  if (current$Fit < previous$Fit){
    previous <- current
    next
  }
  else {
    d <- current$Fit - previous$Fit
    evalFrac <- 0.05
    alpha <- .999
    
    if (count <= evalFrac*maxCount){
      prob <- alpha^count
      }
    else{
      prob <- (alpha^count)/(1+ (d/d.avg)^2)
      }
    change <- sample(c(0,1), size = 1, prob=c(1-prob, prob))
    if (change == 1) {
      previous <- current
    } 
  }
  
  count <- count +1
  if (count == 1000000) {break}
}
###End of Loop###

print(paste0("Count: ", count))
print (best) 

output <- c("Fit: ", best$Fit, "Iterations: ", count)
write.table(output, file = "JPS_Output.csv", append = TRUE)
}










