# -----------------------
# GET SIMULATION RESULTS
# -----------------------

library("ProSGPV")
library("MASS")
library("BeSS") # USE VERSION 2.0.0 WHICH SUPPORTS POISSON
library("SIS")
library("ggplot2")
library("grid")
library("ggpubr")
library("wesanderson")
library("survival")
library("pROC")
library("doParallel")

source("utils.R")

# -----------------------
# FIGURE 1
# -----------------------


# TO DELETE LATER

i = 45
while(T){
  
  set.seed(i)
  data <- gen.data.fig.1()
  
  # load data
  x <- data[[1]]
  y <- data[[2]]
  
  n <- length(y)
  
  # find the lambda cutoffs
  lasso.poisson <- glmnet(x,y,family = "poisson")
  
  yhat.m <- predict(lasso.poisson,newx=x,type="response")
  p.k <- lasso.poisson$df + 1
  
  lasso.lambda.index <- which.min(sapply(1:ncol(yhat.m),function(z)
    poisson()$aic(y, rep(1,n), yhat.m[,z], rep(1,n)) +
      p.k[z] * log(log(n))*log(p.k[z]) ) )
  lambda.gic <- lasso.poisson$lambda[lasso.lambda.index]
  candidate.index <- which(
    coef(lasso.poisson, s = lambda.gic)[-1] != 0)
  
  if(length(candidate.index) == 2){
    
    sgpv.index <- pro.sgpv(x,y,family="poisson")$var.index
    
    if(setequal(sgpv.index,3)) break
  }
  
  i = i + 1
}

i

set.seed(72)
fig.1.data <- gen.data.fig.1()

fig.1 <- get.plot.fig.1(fig.1.data, lambda.max = 0.25, x.offset = -0.01, y.p = -0.7)

# -----------------------
# OTHER MAIN FIGURES
# -----------------------

set.seed(1)

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
registerDoParallel(no_cores)

# logistic low-d sparse
out.main.logit.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = n, p = 20, s = 4, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# logistic low-d dense
out.main.logit.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = n, p = 20, s = 14, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# logistic high-d sparse
out.main.logit.hd.s <- foreach(
  p = seq(200, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = 200, p = p, s = 4, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# poisson low-d sparse
out.main.poisson.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = n, p = 20, s = 4, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# poisson low-d dense
out.main.poisson.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = n, p = 20, s = 14, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# poisson high-d sparse
out.main.poisson.hd.s <- foreach(
  p = seq(120, 480, 24),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = 120, p = p, s = 10, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# cox low-d sparse
out.main.cox.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = n, p = 20, s = 4, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}

# cox low-d dense
out.main.cox.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = n, p = 20, s = 14, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}

# cox high-d sparse
out.main.cox.hd.s <- foreach(
  p = seq(80, 320, 16),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  main.many(
    n = 80, p = p, s = 4, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}

# ----------------------------------------------
# SUPPLEMENTARY FIGURE 1: SENSITIVITY OF LAMBDA
# ----------------------------------------------

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
registerDoParallel(no_cores)

# logistic low-d sparse
out.supp.fig.1.logit.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = n, p = 20, s = 4, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# logistic low-d dense
out.supp.fig.1.logit.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = n, p = 20, s = 14, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# logistic high-d sparse
out.supp.fig.1.logit.hd.s <- foreach(
  p = seq(200, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = 200, p = p, s = 4, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# poisson low-d sparse
out.supp.fig.1.poisson.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = n, p = 20, s = 4, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# poisson low-d dense
out.supp.fig.1.poisson.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = n, p = 20, s = 14, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# poisson high-d sparse
out.supp.fig.1.poisson.hd.s <- foreach(
  p = seq(120, 480, 24),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = 120, p = p, s = 10, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# cox low-d sparse
out.supp.fig.1.cox.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = n, p = 20, s = 4, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}

# cox low-d dense
out.supp.fig.1.cox.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = n, p = 20, s = 14, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}

# cox high-d sparse
out.supp.fig.1.cox.hd.s <- foreach(
  p = seq(80, 320, 16),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.1(
    n = 80, p = p, s = 4, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}


# -----------------------------------
# SUPPLEMENTARY FIGURE 2: NULL BOUND
# -----------------------------------

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
registerDoParallel(no_cores)

supp.fig.2.ind <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.2(n = n, rho = 0)
}

supp.fig.2.medium <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.2(n = n, rho = 0.35)
}

supp.fig.2.high <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.2(n = n, rho = 0.7)
}

# ----------------------------------------
# SUPPLEMENTARY FIGURE 3: ML AND JEFFREYS
# ----------------------------------------

set.seed(1)

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
registerDoParallel(no_cores)

# logistic low-d sparse
out.supp.fig.3.low.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  supp.fig.3.many(
    n = n, p = 20, s = 4, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# logistic low-d dense
out.supp.fig.3.low.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  supp.fig.3.many(
    n = n, p = 20, s = 14, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# logistic high-d sparse
out.supp.fig.3.high.s <- foreach(
  p = seq(200, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  supp.fig.3.many(
    n = 200, p = p, s = 4, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# ---------------------------------------------
# SUPPLEMENTARY FIGURE 6: LOW-D DENSE LOGISTIC
# ---------------------------------------------

set.seed(1)

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
registerDoParallel(no_cores)

# rho = 0. 35
out.supp.fig.6.1 <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.6(
    n = n, rho = 0.35
  )
}

# rho = 0.7
out.supp.fig.6.2 <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.6(
    n = n, rho = 0.7
  )
}

# ----------------------------------------
# SUPPLEMENTARY FIGURE 7: COMPUTATION TIME
# ----------------------------------------

set.seed(1)

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
registerDoParallel(no_cores)

# logistic low-d sparse
out.supp.fig.7.logit.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = n, p = 20, s = 4, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# logistic low-d dense
out.supp.fig.7.logit.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = n, p = 20, s = 14, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# logistic high-d sparse
out.supp.fig.7.logit.hd.s <- foreach(
  p = seq(200, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = 200, p = p, s = 4, beta.min = 0.5, beta.max = 1.5,
    intercept = 0, family = "binomial"
  )
}

# poisson low-d sparse
out.supp.fig.7.poisson.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = n, p = 20, s = 4, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# poisson low-d dense
out.supp.fig.7.poisson.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = n, p = 20, s = 14, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# poisson high-d sparse
out.supp.fig.7.poisson.hd.s <- foreach(
  p = seq(120, 480, 24),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = 80, p = p, s = 10, beta.min = 0.1, beta.max = 0.4,
    intercept = 2, family = "poisson"
  )
}

# cox low-d sparse
out.supp.fig.7.cox.ld.s <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = n, p = 20, s = 4, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}

# cox low-d dense
out.supp.fig.7.cox.ld.d <- foreach(
  n = seq(40, 800, 40),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = n, p = 20, s = 14, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}

# cox high-d sparse
out.supp.fig.7.cox.hd.s <- foreach(
  p = seq(80, 320, 16),
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  many.sim.supp.fig.7(
    n = 80, p = p, s = 4, beta.min = 0.2, beta.max = 0.8,
    intercept = 0, family = "cox"
  )
}

# --------------------
# REAL WORLD ANALYSIS
# --------------------

set.seed(1)

spine.temp <- replicate(1e3, suppressMessages(
  spine.one.time(spine[,1:12],spine[,13]) ))

# most frequent model
# sgpv: 4 5 6
sort(table(paste(spine.temp[1, ])), decreasing = T)[1]

# lasso: 2:12 but the frequency is only 80 out of 1000
sort(table(paste(spine.temp[4, ])), decreasing = T)[1]

# bess: 1 2 5 6
sort(table(paste(spine.temp[7, ])), decreasing = T)[1]

# sis: 3
sort(table(paste(spine.temp[11, ])), decreasing = T)[1]

# get histogram

hist.d <- data.frame(
  Size = unlist(spine.temp[c(2,5,8,11), ]),
  Algorithm = rep(c("ProSGPV", "Lasso", "BeSS", "ISIS"), 1e3)
)
hist.d$Algorithm <- factor(hist.d$Algorithm,
                           levels = c("ProSGPV", "Lasso", "BeSS", "ISIS")
)

cols <- c("black", "springgreen3", "blue", "red")

# histogram
spine.hist <- ggplot(hist.d, aes(x = Size, fill = Algorithm)) +
  geom_density(alpha = 0.6, color = NA) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  scale_fill_manual(values = cols) +
  labs(x = "Model size", y = "Density", title = "Spine training set (70%)") +
  annotate("text",
           x = 0, y = 5,
           label = paste(c(
             "Lasso most often selects", paste("V", 2:12, sep = "")
           ),
           collapse = " "
           ), hjust = 0, col = "springgreen3"
  ) +
  annotate("text",
           x = 0, y = 4,
           label = paste(c(
             "BeSS most often selects", 
             paste("V", c(1, 2, 5, 6) , sep = "")
           ),
           collapse = " "
           ), hjust = 0, col = "blue"
  ) +
  annotate("text",
           x = 0, y = 3, col = "red",
           label =  "ISIS most often selects V3" , hjust = 0
  ) +
  annotate("text",
           x = 0, y = 6, col = "black",
           label = paste(c(
             "ProSGPV most often selects", paste("V", 4:6, sep = "")
           ),
           collapse = " "
           ), hjust = 0
  )


# get boxplot of prediction accuracy

box.d <- data.frame(
  error = unlist(spine.temp[c(3,6,9,12), ]),
  Method = rep(c("ProSGPV", "Lasso", "BeSS", "ISIS"
  ), 1e3)
)
box.d$Method <- factor(box.d$Method,
                       levels = c("ProSGPV", "Lasso", "BeSS", "ISIS")
)

cols <- c("black", "springgreen3", "blue", "red")

# boxplot 
spine.box <- ggplot(box.d, aes(x = Method, fill = Method , y = error, color = Method)) +
  geom_boxplot(alpha = 0.6) +
  theme_classic() +  theme(legend.position = "none") +
  scale_fill_manual(values = cols) + 
  scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
  scale_color_manual(values = cols) +
  labs(x = "Method", y = "Prediction AUC", title = "Spine test set (30%)") 



# -----------------
# SAVE ALL RESULTS
# -----------------
load("manuscript_data.RData")

save(fig.1,

  out.main.logit.ld.s,
  out.main.logit.ld.d,
  out.main.logit.hd.s,
  
  out.main.poisson.ld.s,
  out.main.poisson.ld.d,
  out.main.poisson.hd.s,

  out.main.cox.ld.s,
  out.main.cox.ld.d,
  out.main.cox.hd.s,

  out.supp.fig.1.logit.ld.s,
  out.supp.fig.1.logit.ld.d,
  out.supp.fig.1.logit.hd.s,

  out.supp.fig.1.poisson.ld.s,
  out.supp.fig.1.poisson.ld.d,
  out.supp.fig.1.poisson.hd.s,

  out.supp.fig.1.cox.ld.s,
  out.supp.fig.1.cox.ld.d,
  out.supp.fig.1.cox.hd.s,
  
  supp.fig.2.ind,
  supp.fig.2.medium,
  supp.fig.2.high,

  out.supp.fig.3.low.s,
  out.supp.fig.3.low.d,
  out.supp.fig.3.high.s,
  
  out.supp.fig.6.1,
  out.supp.fig.6.2,
  
  out.supp.fig.7.logit.ld.s,
  out.supp.fig.7.logit.ld.d,
  out.supp.fig.7.logit.hd.s,
  
  out.supp.fig.7.poisson.ld.s,
  out.supp.fig.7.poisson.ld.d,
  out.supp.fig.7.poisson.hd.s,
  
  out.supp.fig.7.cox.ld.s,
  out.supp.fig.7.cox.ld.d,
  out.supp.fig.7.cox.hd.s,
  
  spine.hist,
  spine.box,
  
  file = "manuscript_data.RData"
)


