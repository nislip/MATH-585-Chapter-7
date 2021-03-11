### Packages ==================================================

install.packages("matlib")
library(matlib)
install.packages("stargazer")
library(stargazer)
install.packages("ggplot2")
library(ggplot2)
install.packages("latex2exp")
library(latex2exp)
install.packages("ellipse")
library(ellipse)
install.packages("MASS")
library(MASS)
install.packages("mvtnorm")

### Problem 7.9 =====================================================
{
# construct data matrix Z
df1 <- Problem_7_9
n1 <- nrow(df1)
Z1 <- matrix(c(rep(1,n1), df1$z1), ncol=2)
ZZinv1 <- inv(t(Z1) %*% Z1)
# Compute least squares estimates
y1 <- matrix(c(df1$Y1), ncol = 1)
y2 <- matrix(c(df1$Y2), ncol = 1)

# Z'y1 and Z'y2

Zy1 <- t(Z1) %*% y1
Zy2 <- t(Z1) %*% y2

# Least squares estimators (Beta = (Z'Z)-1 Z'y(1,2))

Beta1 <- ZZinv1 %*% Zy1 # Individual least squares estimators
Beta2 <- ZZinv1 %*% Zy2

BETA <- cbind(Beta1, Beta2) # Combined least squares estimators 
Yhat <- Z1 %*% BETA # Yhat = Z * Beta hat 

# Now for epsilon hat (residuals)

Y <- cbind(df1$Y1, df1$Y2)
ehat <- Y - Yhat
YY <- t(Y) %*% Y
YYhat <- t(Yhat) %*% Yhat
eehat <- t(ehat) %*% ehat 
test <- YYhat + eehat # same as YY therefore verified
}

# Verify the sum of squares and cross products decomposition

# LMz1 <- lm(df1$Y1 ~ df1$z1, data = df1)
# LMz2 <- lm(df1$Y2 ~ df1$z1, data = df1)

### Problem 7.15 =====================================================

df <- table_71
linear_model <- lm(df$y ~ df$z1 + df$z2, data=df)
summary(linear_model)
n <- nrow(df)
stargazer(linear_model)

Z <- matrix(c(rep(1,n), df$z1, df$z2), ncol = 3)
Y <- matrix(c(df$y), ncol = 1)

ZZinv <- inv(t(Z) %*% Z)
Beta <- ZZinv %*% t(Z) %*% Y

# Plotting the residuals 

model_z1 <- lm(df$y ~ df$z1)
model_z2 <- lm(df$y ~ df$z2)
res <- residuals(model_z1)
res2 <- residuals(model_z2)
plot(res, main = "Residuals with Z1", ylab = "Residuals", xlab = "z1", xlim = c(0,25), ylim = c(-7,7))
abline(h = 0, col = "blue")
plot(res2, main = "Residuals with Z2", ylab = "Residuals", xlab = "z2")
abline(h = 0, col = "blue")

# Plotting the residuals for different scenarios 
df_test <- c(df$z1 * df$z2)
df_test <- data.frame(df_test)

# Plotting the residuals as a histogram

ggplot(data = df, aes(x = linear_model$residuals)) +
  geom_histogram(fill = 'steelblue', color = 'black') +
  labs(title = 'Histogram of Residuals', x = 'Residuals', y = 'Frequency')


### Problem 8.2 =====================================================

E <- matrix(c(5,2,2,2), nrow = 2, ncol = 2)
Rho <- cov2cor(E)
eigen(Rho)


### Problem 8.6 ======================================================

### Part (A)

# Construct S

Xbar <- c(155.60, 14,70)
S <- matrix(c(7476.45, 303.62, 303.62, 26.19), nrow = 2, ncol=2 )
# Compute eigenvalues and eigenvectors
E <- eigen(S)

(-0.99917337 * sqrt(13.83))/(S[2,2])

cor(sales_versus_profits_chapter_8)

# Ellipse 
p <- ellipse(Xbar, S, alpha = 0.05, npoints = 250)
ellipse(p)


### Problem 8.10 =====================================================

### Part A) 

# Construct Sample Covariance matrix 

WRR <- weekly_rates_of_return_T_84
S_WRR <- cov(WRR) # Sample Covariance Matrix
col_mean_WRR <- colMeans(WRR) # Column means

stargazer(S_WRR, digits = 10) # Output to LaTex table

E_S <- eigen(S_WRR) # Construct Eigen values & vectors
E_values <- E_S[["values"]] # Eigen Values
E_vectors <- E_S[["vectors"]] # Eigen Vectors
stargazer(E_values, digits = 7) # Output LaTex
stargazer(E_vectors, digits = 7) # Output Latex

pcal <- princomp(WRR, scores = TRUE, cor = FALSE)

# Part C) Construct Confidence Intervals

z <- qnorm(p=0.05, lower.tail = FALSE)
val <- sqrt(2/103)

B_Interval_1 <- E_values[1] / (1 + z*val)
B_Interval_12 <- E_values[1] / (1 - z*val)

B_Interval_2 <- E_values[2] / (1 + z*val)
B_Interval_22 <- E_values[2] / (1 - z*val)

B_Interval_3 <- E_values[3] / (1 + z*val)
B_Interval_32 <- E_values[3] / (1 - z*val)

### Problem 8.26 ====================================================

# Part A) 

df <- t46
X <- cbind(df$indep, df$supp, df$Benev, df$conform, df$leader) 
R <- cor(X)
EE <- eigen(R)
E_values_1 <- EE[["values"]]
E_vectors_1 <- EE[["vectors"]]

pcal <- princomp(X)
PCA_test <- summary(pcal)
plot(pcal)
screeplot(pcal, type="line", main = "Scree Plot")

# Part C) 

# Test with Table 5.8 data 

PP <- t58
pca_police <- princomp(PP, cor = TRUE)

pcal_fixed <- princomp(X)
df4 <- PCA_df_scores
data.frame(df4)
plot(y = df4$Comp.1, x = df4$Comp.2, ylab = TeX('$\\alpha x^\\alpha'))
abline(h = 0)

ggplot(df4, aes(x = df4$Comp.1, y = df4$Comp.2, colour = Comp.1, color = Comp.2)) + geom_point() + xlab(TeX('$\\hat{y}_{2}')) + ylab(TeX('$\\hat{y}_{1}')) + geom_hline(yintercept = 0, linetype = "dashed") + ggtitle(TeX('Plotting the scores of $\\hat{y}_{1}$ and $\\hat{y}_{2}')) + stat_ellipse()  

# Part D) 

z <- qnorm(p=0.025, lower.tail = FALSE)
CI_1 <- E_values_1[1]/(1 + z * sqrt(2/130))
CI_2 <- E_values_1[1]/(1 - z * sqrt(2/130))










sigma <- matrix(c(1, .3, .3, 1.5), 2, 2)
mu <- c(1, 3)

p <- ellipse(mu, sigma, npoints = 200, draw = TRUE)
plot(p)
