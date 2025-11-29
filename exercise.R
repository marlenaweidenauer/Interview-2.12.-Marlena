#-------------------------------------------------------------
#------------------ Definition variables ---------------------
#-------------------------------------------------------------

beta <- 24
delta_E <- 0.6
mu_E <- 0.15
delta_J <- 0.08
mu_J <- 0.05
alpha <- 0.003
omega <- 0.5
mu_A <- 0.1
t_0 <- 0
t_end <- 365
h <- 0.01
y_0 <- c(10,0,0)

#-------------------------------------------------------------
#----------------------- Algorithms --------------------------
#-------------------------------------------------------------

# Algorithm 1 - function RK4
RK4 <- function (f, y_n, t_n, h){ 
  k_1 <- f(t_n,y_n)
  k_2 <- f(t_n + 1/2*h, y_n + 1/2*h*k_1)
  k_3 <- f(t_n + 1/2*h, y_n + 1/2*h*k_2)
  k_4 <- f(t_n + h, y_n + h*k_3)
  y_new <- y_n + 1/6 *h*(k_1 + 2*k_2 + 2*k_3 + k_4)
  return(y_new)
}

# Algorithm 2 - function f with predefined variables
f <- function(t_n, y_n){
  E <- y_n[1]
  J <- y_n[2]
  A <- y_n[3]
  dE <- beta*A - delta_E*E - mu_E*E
  dJ <- delta_E*E - delta_J*J- alpha*J^2 - mu_J*J
  dA <- omega*delta_J*J - mu_A*A
  dy <- c(dE,dJ,dA)
  return(dy)
}


#-------------------------------------------------------------
#------------------  Vector T --------------------------------
#-------------------------------------------------------------

T <- seq(t_0, t_end, by = h) # vector of length 36501

#-------------------------------------------------------------
#------------------ for-loop ---------------------------------
#-------------------------------------------------------------

n <- t_end/h # number of iterations

# matrix of y_0 and simulated observations
Y <- matrix(0, nrow = 3, ncol = n + 1)
Y[,1] <- y_0

for (i in 1:n){
  Y[,i+1] <- RK4(f,Y[,i],T[i],h)
}

#-------------------------------------------------------------
#------------------ Plot -------------------------------------
#-------------------------------------------------------------

par(mfrow = c(3, 1))   # three subplots

plot(T, Y[1, ], type = "l", xlab = "T", ylab = "E", col = "darkgreen")
plot(T, Y[2, ], type = "l", xlab = "T", ylab = "J", col = "darkorange")
plot(T, Y[3, ], type = "l", xlab = "T", ylab = "A", col = "darkblue")

