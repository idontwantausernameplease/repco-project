size <- 99

c <- toeplitz(0:size) #creates a matrix with size+1 rows and columns where first row goes from 0-size - original had size = 9

# from excel file 
# h_share <- c(0.05,0.075,0.1,0.15,0.15,0.15,0.15,0.1,0.05,0.025)
# e_share <- c(0.6,0.15,0.025,0.025,0.025,0.025,0.025,0.05,0.05,0.025)

h_share <- diff(c(0, sort(runif(size)), 1)) #creates an array of [size] random numbers in a random sequence that all add up to 1
e_share <- diff(c(0, sort(runif(size)), 1)) #creates an array of [size] random numbers in a random sequence that all add up to 1

gamma <- -0.07

# alpha <- numeric(10)
# beta <- numeric(10)

alpha <- numeric(size+1)
beta <- numeric(size+1)

get_beta <- function(alpha) {
  a_matrix <- exp(alpha+t(c)*gamma)
  #print(a_matrix)
  beta <- log(e_share)-log(colSums(a_matrix))
  beta_1 <<- log(e_share)-log(colSums(a_matrix))
  beta <<- beta-beta[[1]]
}

get_alpha <- function(beta) {
  b_matrix <- t(exp(beta+t(c)*gamma))
  #print(b_matrix)
  alpha <- log(h_share)-log(rowSums(b_matrix))
  alpha <<- alpha-alpha[[1]]
}

get_pi_ij <- function(alpha,beta) {
  alphas <- matrix(alpha, nrow=length(alpha), ncol=length(alpha), byrow=FALSE)
  #print(alphas)
  betas <- matrix(beta_1, nrow=length(beta_1), ncol=length(beta_1), byrow=TRUE)
  #print(betas)
  pi_ij <<- exp(alphas+betas+gamma*c)
}

get_commutes <- function(pi_ij) {
  commute <- pi_ij*c
  commute <<- pi_ij*c
  avgcommute <<- sum(commute)
}

iteration <- 0
precision <- 0.0001
repeat { #iterates through the above functions until the right alphas and betas have been found
  get_beta(alpha)
  #print(beta)
  get_alpha(beta)
  #print(alpha)
  get_pi_ij(alpha,beta)
  #print(pi_ij)
  get_commutes(pi_ij)
  #print(commute)
  print(paste0("average commute: ", avgcommute))
  iteration <- iteration + 1
  print(paste0("iterations: ", iteration))
  if(abs(h_share-rowSums(pi_ij))<precision && abs(e_share-colSums(pi_ij))<precision) 
    break
} 

