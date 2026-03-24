library(MASS)
library(ggplot2)
library(dplyr)

load(file = "data/speaker1_tongue_shapes.RData")
load(file = "data/speaker2_tongue_shapes.RData")

dim(speaker1)
dim(speaker2)

n_shapes1 <- 10
n_shapes2 <- 10

n_points1 <- 65
n_points2 <- 45

##############
# Data preparation (hopefully has everything)
# in theory, set seed too for simulations!

# first speaker
df1 <- data.frame(
  x = as.vector(speaker1[,1,]),
  y = as.vector(speaker1[,2,]),
  speaker = 1,
  shape = rep(1:n_shapes1, each = n_points1),  # number
  point = rep(1:n_points1, times = n_shapes1)  # point along the curve
)

# second speaker
df2 <- data.frame(
  x = as.vector(speaker2[,1,]),
  y = as.vector(speaker2[,2,]),
  speaker = 2,                                 # Updated to 2
  shape = rep(1:n_shapes2, each = n_points2),  # Updated to n_shapes2 and n_points2
  point = rep(1:n_points2, times = n_shapes2)  # Updated to n_points2 and n_shapes2
)

# check them
ggplot(df1, aes(x = x, y = y, group = point)) +
  geom_point(size = 1) +                       
  facet_wrap(~shape, ncol = 5) +      
  coord_fixed() +                    
  theme_minimal()+
  labs(title = "Speaker 1 Tongue Shapes (10 images, first two have missing points)")

ggplot(df2, aes(x = x, y = y, group = point)) +
  geom_point(size = 1) +                       
  facet_wrap(~shape, ncol = 5) +      
  coord_fixed() +                    
  theme_minimal()+
  labs(title = "Speaker 2 Tongue Shapes (10 images, first two have missing points)")

# MARGINAL CURVES
# x curve

ggplot(df1, aes(x = point, y = x, group = point)) +
  geom_point(size = 1) +                       
  facet_wrap(~shape, ncol = 5) +      
  coord_fixed() +                    
  theme_minimal()+
  labs(title = "Speaker 1 Tongue Shapes (first 10 images, first two have missing points)")
# y curve
ggplot(df1, aes(x = point, y = y, group = point)) +
  geom_point(size = 1) +                       
  facet_wrap(~shape, ncol = 5) +      
  coord_fixed() +                    
  theme_minimal()+
  labs(title = "Speaker 1 Tongue Shapes (10 images, first two have missing points)")

# y, curve 2 has a slightly weird point. may keep, may remove (currently not removed)


# speaker 2
# x curve
ggplot(df2, aes(x = point, y = x, group = point)) +
  geom_point(size = 1) +                       
  facet_wrap(~shape, ncol = 5) +      
  coord_fixed() +                    
  theme_minimal()+
  labs(title = "Speaker 2 Tongue Shapes (10 images, first two have missing points)")
# y
ggplot(df2, aes(x = point, y = y, group = point)) +
  geom_point(size = 1) +                       
  facet_wrap(~shape, ncol = 5) +      
  coord_fixed() +                    
  theme_minimal()+
  labs(title = "Speaker 2 Tongue Shapes (10 images, first two have missing points)")


####
# reorder the second speaker curves
# (speaker 2 has 45 points)
df2$point2 <- 45 - df2$point + 1

# but the 5th curve is ok, it just has a weird point (so reorder 5th again and remove that point)
df2 <- df2 %>%
  mutate(point2 = ifelse(shape == 5,
                         45 - point2 + 1,
                         point2)) %>% 
  filter(!(shape == 5 & point2 == 45)) # remove the odd point

# compare
ggplot(df2, aes(x = point2, y = x, group = point2)) +
  geom_point(size = 1) +                       
  facet_wrap(~shape, ncol = 5) +      
  coord_fixed() +                    
  theme_minimal()+
  labs(title = "Speaker 2 Tongue Shapes (10 images, first two have missing points)")
# should be s-shaped
# y
ggplot(df2, aes(x = point2, y = y, group = point2)) +
  geom_point(size = 1) +                       
  facet_wrap(~shape, ncol = 5) +      
  coord_fixed() +                    
  theme_minimal()+
  labs(title = "Speaker 2 Tongue Shapes (10 images, first two have missing points)")


##############


# function that takes a shape and calculate the "betas" (coefficients for the polynomial basis)
coefs <- function(x, y, degree){
  
  # design matrix
  #X_poly <- cbind(1, x, x^2, x^3, x^4, x^5)
  X_poly <- sapply(0:degree, function(i) x^i)
  
  m <- lm(y ~ X_poly - 1)
  beta <- coef(m)
  
  return(beta)
}



# function that receives coefficients and a curve and will try to find the 6 best betas
curve_selection <- function(X_observed, y_observed, number = 6){
  
  mat <- matrix(NA, nrow = ncol(X_observed), ncol = 2)
  for (i in 1:ncol(X_observed)){
    
    similarity <- sum((X_observed[,i] - y_observed)^2)
    mat[i,1] <- i
    mat[i,2] <- similarity
  }
  
  # sort by the second column (similarity) in ascending order
  mat <- mat[order(mat[,2]), ]
  
  # only return indices
  return(mat[1:number,1])
  
}

# function that takes a curve and the coefficients and fits it

s <- function(x, y, coefs, observed, simulations){
  
  # this was used to test the internal parts
  #x <- 1:65
  #coefs <- betas2
  #observed <- 21:65
  #y <- subset(df1, shape == 1)$x
  #simulations <- 5000
  
  # build design matrix
  #X_poly <- cbind(1, x,x^2,x^3,x^4,x^5)
  X_poly <- sapply(0:(ncol(coefs)-1), function(i) x^i)
  X_full <- t((coefs) %*% t(X_poly))
  
  # choose the observed parts
  y_observed <- y[observed]
  X_observed <- X_full[observed, ]
  n <- length(y_observed)
  
  # curve selection
  selected <- curve_selection(X_observed, y_observed, ncol(coefs))
  X_observed <- X_observed[, selected]
  X_full <- X_full[,selected]
  
  
  # initialise parameters for priors
  p <- ncol(X_observed)
  beta0 <- rep(0.5, p)
  Sigma0_inv <- 1/0.5* diag(1,p)
  #Sigma0 <- diag(1, p)
  
  # choose 
  nu0 <- 2
  s20 <- 5
  
  beta_samples <- matrix(NA, nrow=simulations, ncol=p)
  sigma2_samples <- numeric(simulations)
  
  # initial values
  beta_samples[1, ] <- rep(1/p, p)
  sigma2_samples[1] <- var(y_observed)
  
  for(iter in 2:simulations){
    
    # BETA
    # find V
    V_beta <- solve(t(X_observed) %*% X_observed / sigma2_samples[iter-1] + Sigma0_inv)
    # find m
    m_beta <- V_beta %*% (t(X_observed) %*% y_observed / sigma2_samples[iter-1] + Sigma0_inv %*% beta0)
    # sample beta
    #beta_samples[iter, ] <- m_beta + t(chol(V_beta)) %*% rnorm(p)
    # should be the same as 
    beta_samples[iter, ] <- mvrnorm(1, mu = m_beta, Sigma = V_beta)
    
    # SIGMA
    
    # first calculate this
    SSR <- sum((y_observed - X_observed %*% beta_samples[iter, ])^2)
    
    # parameters
    shape_sig <- (nu0 + n) / 2
    rate_sig  <- (nu0 * s20 + SSR) / 2
    
    sigma2_samples[iter] <- 1 / rgamma(1, shape = shape_sig, rate = rate_sig)
  }
  
  predictions <- X_full %*% t(beta_samples)
  
  return(list(
    predictions = predictions,
    beta_samples = beta_samples,
    selected = selected
  ))
  
}

############################################
# important:
# for speaker 1, curve 1 observed points are (obs_id <- 21:65)
# for speaker 1, curve 2 observed points are (obs_id <- c(1:19, 41:65))

# trying it out
# SPEAKER 1, curve 1

# putting together the coefficients matrix
betas_list <- lapply(3:10, function(i){
  
  curve_i <- subset(df1, shape == i)
  y <- curve_i$x
  
  coefs(1:65, y, 5) # up to degree 5!
})

betas2 <- do.call(rbind, betas_list)
betas2

# simply call this function
b1 <- s(1:65, subset(df1, shape == 1)$x, betas2, 21:65, 5000)

y_pred_samples1 <- b1$predictions

results1 <- data.frame(
  point = 1:65,
  best_estimate = rowMeans(y_pred_samples1),
  lower_95 = apply(y_pred_samples1, 1, quantile, probs = 0.025),
  upper_95 = apply(y_pred_samples1, 1, quantile, probs = 0.975))

head(results1)

# first curve
plot(1:65, results1$best_estimate)
points(1:65, subset(df1, df1$shape == 1)$x, col = "red")
lines(1:65, results1$lower_95)
lines(1:65, results1$upper_95)

# y curve
betas_list <- lapply(3:10, function(i){
  
  curve_i <- subset(df1, shape == i)
  y <- curve_i$y
  
  coefs(1:65, y,4) # up to degree 4! (5 is too flexible!)
})

betas2 <- do.call(rbind, betas_list)
betas2

b2 <- s(1:65, subset(df1, shape == 1)$y, betas2, 21:65, 5000)

y_pred_samples2 <- b2$predictions
results2 <- data.frame(
  point = 1:65,
  best_estimate = rowMeans(y_pred_samples2),
  lower_95 = apply(y_pred_samples2, 1, quantile, probs = 0.025),
  upper_95 = apply(y_pred_samples2, 1, quantile, probs = 0.975)
)

head(results2)

# first curve y
plot(1:65, results2$best_estimate)
points(1:65, subset(df1, df1$shape == 1)$y, col = "red")
lines(1:65, results2$lower_95)
lines(1:65, results2$upper_95)

###################
# looking at (x,y) together
plot(results1$best_estimate, results2$best_estimate)
# comparing
points(subset(df1, df1$shape == 1)$x, subset(df1, df1$shape == 1)$y, col = "red")
points(subset(df1, df1$shape == 5)$x, subset(df1, df1$shape == 5)$y, col = "blue")
points(subset(df1, df1$shape == 2)$x, subset(df1, df1$shape == 2)$y, col = "green")
# they are all slightly different curves; perhaps the 2nd matches the part where the 1st has missing points the best



