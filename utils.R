# --------------------
# UTILITY FUNCTIONS
# --------------------

# --------------------
# FIGURE 1
# --------------------

# function to generate data in figure 1
gen.data.fig.1 <- function() {
  beta <- c(0, 0, 0.25, 0, 0)
  intercept <- 0

  cov.structure <- matrix(0, 5, 5)
  for (i in 1:5) {
    for (j in 1:5) {
      cov.structure[i, j] <- 0.5^(abs(i - j))
    }
  }

  # generate X
  X <- mvrnorm(n = 100, rep(0, 5), Sigma = cov.structure)

  # generate z
  z <- X %*% beta + intercept

  # generate Y
  Y <- rpois(n = 100, exp(z))

  return(list(X, Y))
}


# function to extract coefficients in figure 1
get.coef.fig.1 <- function(x, y, lambda, lasso) {
  p <- ncol(x)

  # evaluate lasso at lambda
  pe.lasso <- coef(lasso, s = lambda)[-1]
  index <- which(pe.lasso != 0)

  # define the output
  out.coef <- numeric(p)
  out.lb <- numeric(p)
  out.ub <- numeric(p)

  if (lambda == 0) {

    # full glm model
    full.glm <- glm(y ~ x, family = "poisson")

    pe <- summary(full.glm)$coef[-1, 1]
    se <- summary(full.glm)$coef[-1, 2]
    lb <- pe - 1.96 * se
    ub <- pe + 1.96 * se

    null.bound.lasso <- mean(summary(full.glm)$coef[-1, 2])

    out.coef <- as.numeric(pe)
    out.lb <- as.numeric(lb)
    out.ub <- as.numeric(ub)
  } else if (length(index) != 0) {

    # run fully relaxed LASSO
    f.l <- glm(y ~ x[, index], family = "poisson")

    # get confidence bands
    pe <- summary(f.l)$coef[-1, 1]
    se <- summary(f.l)$coef[-1, 2]
    lb <- pe - 1.96 * se
    ub <- pe + 1.96 * se

    null.bound.lasso <- mean(summary(f.l)$coef[-1, 2])

    out.coef[index] <- as.numeric(pe)
    out.lb[index] <- as.numeric(lb)
    out.ub[index] <- as.numeric(ub)
  } else if (length(index) == 0) {

    # intercept only model
    null.bound.lasso <- 0
  }

  return(c(pe.lasso, out.coef, out.lb, out.ub, null.bound.lasso))
}


# output figure 1
get.plot.fig.1 <- function(data, p = 5, lambda.max = 0.2, 
                           x.offset = -0.01, y.p = -0.5) {

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
  
  # lasso selection results
  q.index <- which( coef(lasso.poisson, s = lambda.gic)[-1] != 0)

  # get a range of lambda
  lambda.grid <- seq(0, lambda.max, length.out = 100)

  # input lambda and output coefficients
  lasso <- glmnet(x, y, family = "poisson")
  results <- sapply(
    lambda.grid,
    function(z) get.coef.fig.1(x, y, lambda = z, lasso)
  )

  # create a data set for plot
  to.plot <- data.frame(
    lambda = rep(lambda.grid, each = p),
    v = rep(paste("V", 1:p, sep = ""), length(lambda.grid)),
    pe.lasso = c(results[1:p, ]),
    pe.ols = c(results[(p + 1):(2 * p), ]),
    lb = c(results[(2 * p + 1):(3 * p), ]),
    ub = c(results[(3 * p + 1):(4 * p), ])
  )

  to.plot$closer <- sapply(
    1:nrow(to.plot),
    function(x) {
      ifelse(abs(to.plot$lb[x]) <
        abs(to.plot$ub[x]), "lb", "ub")
    }
  )

  to.plot$curve <- sapply(
    1:nrow(to.plot),
    function(x) to.plot[x, to.plot$closer[x]]
  )

  # find the limit of the canvas
  ylim <- c(
    min(c(to.plot$pe.ols, to.plot$lb, to.plot$ub)) * 1.1,
    max(c(to.plot$pe.ols, to.plot$lb, to.plot$ub)) * 1.1
  )

  # panel 1: lasso

  # change color
  color.lines <- c("red", "chocolate4", "blue", "springgreen3", "hotpink")
  color.use <- rep("black", 5)
  color.use[3] <- "blue"

  # find the coefficient in full glm
  location.beta <- to.plot$pe.lasso[to.plot$lambda == 0]

  p1 <- ggplot(data = to.plot) +
    geom_line(aes(x = lambda, y = pe.lasso, col = v)) +
    guides(col = FALSE) +
    scale_color_manual(values = color.lines) +
    scale_x_continuous(
      breaks = seq(0, ceiling(lambda.max * 100) / 100, length.out = 5),
      limits = c(x.offset, lambda.max * 1.01)
    ) +
    labs(x = expression(lambda), y = "Log rate ratio estimates", title = "(1) Lasso") +
    scale_y_continuous(
      breaks = seq(-1, 1, 0.2),
      limits = ylim
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    annotate("text",
      x = x.offset, y = location.beta, label = c(
        as.expression(bquote(V[1])),
        as.expression(bquote(V[2])),
        as.expression(bquote(V[3])),
        as.expression(bquote(V[4])),
        as.expression(bquote(V[5]))
      ),
      col = color.use
    ) +
    geom_vline(xintercept = lambda.gic, linetype = "dashed", col = "black") +
    annotation_custom(textGrob(bquote(lambda["gic"])),
      xmin = lambda.gic, xmax = lambda.gic,
      ymin = y.p, ymax = y.p
    ) + coord_cartesian(clip = "off")

  # panel 2: fully relaxed lasso
  location.beta <- to.plot$pe.ols[to.plot$lambda == 0]

  p2 <- ggplot(data = to.plot) +
    geom_line(aes(x = lambda, y = pe.ols, col = v)) +
    guides(col = FALSE) +
    scale_color_manual(values = color.lines) +
    scale_x_continuous(
      breaks = seq(0, ceiling(lambda.max * 100) / 100, length.out = 5),
      limits = c(x.offset, lambda.max * 1.01)
    ) +
    labs(
      x = expression(lambda), y = "Log rate ratio estimates",
      title = "(2) Fully-relaxed lasso"
    ) +
    scale_y_continuous(
      breaks = seq(-1, 1, 0.2),
      limits = ylim
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    annotate("text",
      x = x.offset, y = location.beta, label = c(
        as.expression(bquote(V[1])),
        as.expression(bquote(V[2])),
        as.expression(bquote(V[3])),
        as.expression(bquote(V[4])),
        as.expression(bquote(V[5]))
      ),
      col = color.use
    ) +
    geom_vline(xintercept = lambda.gic, linetype = "dashed", col = "black") +
    annotation_custom(textGrob(bquote(lambda["gic"])),
      xmin = lambda.gic, xmax = lambda.gic,
      ymin = y.p, ymax = y.p
    ) +
    coord_cartesian(clip = "off")

  # panel 3: fully relaxed lasso with confidence interval

  p3 <- ggplot(data = to.plot) +
    geom_line(aes(x = lambda, y = pe.ols, col = v)) +
    guides(col = FALSE) +
    scale_color_manual(values = color.lines) +
    geom_line(aes(x = lambda, y = lb, col = v), alpha = 0.5) +
    geom_line(aes(x = lambda, y = ub, col = v), alpha = 0.5) +
    scale_x_continuous(
      breaks = seq(0, ceiling(lambda.max * 100) / 100, length.out = 5),
      limits = c(x.offset, lambda.max * 1.01)
    ) +
    labs(
      x = expression(lambda), y = "Log rate ratio estimates",
      title = "(3) Fully-relaxed lasso w/ confidence interval"
    ) +
    scale_y_continuous(
      breaks = seq(-1, 1, 0.2),
      limits = ylim
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    annotate("text",
      x = x.offset, y = location.beta, label = c(
        as.expression(bquote(V[1])),
        as.expression(bquote(V[2])),
        as.expression(bquote(V[3])),
        as.expression(bquote(V[4])),
        as.expression(bquote(V[5]))
      ),
      col = color.use
    ) +
    geom_vline(xintercept = lambda.gic, linetype = "dashed", col = "black") +
    annotation_custom(textGrob(bquote(lambda["gic"])),
      xmin = lambda.gic, xmax = lambda.gic,
      ymin = y.p, ymax = y.p
    ) +
    coord_cartesian(clip = "off")

  # panel 4: one bound from fully relaxed lasso + null region
  location.beta <- to.plot$curve[to.plot$lambda == 0]

  # create a data set for null bound
  n.bound <- results[(4 * p + 1), ]
  bound.d <- data.frame(lambda = lambda.grid, lb = -n.bound, ub = n.bound)

  to.plot$v <- factor(to.plot$v, levels = paste("V", 1:p, sep = ""))

  p4 <- ggplot(data = to.plot) +
    geom_line(aes(x = lambda, y = curve, col = v), alpha = 0.5) +
    guides(col = FALSE) +
    scale_color_manual(values = color.lines) +
    geom_ribbon(
      data = bound.d, aes(x = lambda, ymin = lb, ymax = ub),
      alpha = .2
    ) +
    scale_x_continuous(
      breaks = seq(0, ceiling(lambda.max * 100) / 100, length.out = 5),
      limits = c(x.offset, lambda.max * 1.01)
    ) +
    labs(
      x = expression(lambda), y = "Log rate ratio estimates",
      title = "(4) Proposed algorithm"
    ) +
    scale_y_continuous(
      breaks = seq(-1, 1, 0.2),
      limits = ylim
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    annotate("text",
      x = x.offset, y = location.beta, label = c(
        as.expression(bquote(V[1])),
        as.expression(bquote(V[2])),
        as.expression(bquote(V[3])),
        as.expression(bquote(V[4])),
        as.expression(bquote(V[5]))
      ),
      col = color.use
    ) +
    geom_vline(xintercept = lambda.gic, linetype = "dashed", col = "black") +
    annotation_custom(textGrob(bquote(lambda["gic"])),
      xmin = lambda.gic, xmax = lambda.gic,
      ymin = y.p, ymax = y.p
    ) +
    coord_cartesian(clip = "off")

  ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
}




# --------------------
# REST OF THE FIGURES
# --------------------

main.one.time <- function(n, p, s, rho = 0.35, sig = 2,
                          beta.min = 0.1, beta.max = 0.4,
                          intercept = 2,
                          family = c("binomial", "poisson", "cox"),
                          scale = 2, shape = 1, rateC = 0.2) {
  family <- match.arg(family)

  # simulate training data
  sim.data <- gen.sim.data(
    n = round(n * 5 / 3), p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max,
    intercept = intercept, family = family,
    scale = 2, shape = 1, rateC = 0.2
  )
  x <- sim.data[[1]]
  y <- sim.data[[2]]
  true.index <- sim.data[[3]]
  true.beta <- sim.data[[4]]
  noise.index <- setdiff(1:p, true.index)

  # generate training and testing indices
  train.index <- sample(1:round(n * 5 / 3), n, replace = F)
  test.index <- setdiff(1:round(n * 5 / 3), train.index)

  # -------------------
  # lasso
  # -------------------

  if (family != "cox") {
    cv.m <- cv.glmnet(x[train.index, ], y[train.index], family = family)
    lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)[-1]
  } else {
    cv.m <- cv.glmnet(x[train.index, ],
      Surv(y[train.index, 1], y[train.index, 2]),
      family = "cox"
    )
    lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)
  }
  lasso.min.index <- which(lasso.min.coef != 0)

  # out 1: support recovery
  out.lasso.1 <- setequal(true.index, lasso.min.index)

  # out 2: power
  out.lasso.2 <- sum(true.index %in% lasso.min.index) / (length(true.index))

  # out 3: type i error rate
  out.lasso.3 <- sum(lasso.min.index %in% noise.index) / (p - s)

  # out 4: false discovery rate
  out.lasso.4 <- ifelse(length(lasso.min.index) > 0,
    sum(!lasso.min.index %in% true.index) / (length(lasso.min.index)),
    0
  )

  # out 5: false non discovery rate
  out.lasso.5 <- ifelse(length(lasso.min.index) < p,
    sum(setdiff(1:p, lasso.min.index) %in% true.index) /
      (p - length(lasso.min.index)), 0
  )

  # out 6: parameter estimation
  out.lasso.6 <- mean(abs(as.numeric(lasso.min.coef) - true.beta))

  # out 7: prediction performance
  if (family == "binomial") {
    lasso.pred <- predict(cv.m,
      s = cv.m$lambda.min, newx = x[test.index, ],
      type = "response"
    )

    out.lasso.7 <- roc(y[test.index] ~ as.vector(lasso.pred), plot = F, print.auc = F)$auc
  } else if (family == "poisson") {
    lasso.pred <- predict(cv.m,
      s = cv.m$lambda.min, newx = x[test.index, ],
      type = "response"
    )

    out.lasso.7 <- sqrt(mean((lasso.pred - y[test.index])^2))
  } else {
    out.lasso.7 <- 0
  }

  # -------------------
  # pro.sgpv
  # -------------------

  if (family != "cox") {
    out.sgpv <- pro.sgpv(x[train.index, ], y[train.index], family = family)
    sgpv.index <- out.sgpv$var.index
  } else {
    out.sgpv <- try(pro.sgpv(x[train.index, ], y[train.index, ], family = family) , 
        silent=T)
    if(inherits(out.sgpv,"try-error")){
      sgpv.index <- numeric(0)
      
    }else{
      sgpv.index <- out.sgpv$var.index
    }
  }
  

  # out 1: support recovery
  out.sgpv.1 <- setequal(true.index, sgpv.index)

  # out 2: power
  out.sgpv.2 <- sum(true.index %in% sgpv.index) / (length(true.index))

  # out 3: type i error rate
  out.sgpv.3 <- sum(sgpv.index %in% noise.index) / (p - s)

  # out 4: false discovery rate
  out.sgpv.4 <- ifelse(length(sgpv.index) > 0,
    sum(!sgpv.index %in% true.index) / (length(sgpv.index)),
    0
  )

  # out 5: false non discovery rate
  out.sgpv.5 <- ifelse(length(sgpv.index) < p,
    sum(setdiff(1:p, sgpv.index) %in% true.index) /
      (p - length(sgpv.index)), 0
  )

  # out 6: parameter estimation
  if(inherits(out.sgpv,"try-error")){
    out.sgpv.6 <- mean(abs(numeric(p) - true.beta))
  }else{
    out.sgpv.6 <- mean(abs(coef(out.sgpv) - true.beta))
  }

  # out 7: prediction performance
  if (family == "binomial") {
    sgpv.pred <- predict(out.sgpv,
      newdata = x[test.index, ],
      type = "response"
    )

    out.sgpv.7 <- roc(y[test.index] ~ sgpv.pred, plot = F, print.auc = F)$auc
  } else if (family == "poisson") {
    sgpv.pred <- predict(out.sgpv,
      newdata = x[test.index, ]
    )

    out.sgpv.7 <- sqrt(mean((sgpv.pred - y[test.index])^2))
  } else {
    out.sgpv.7 <- 0
  }

  # -------------------
  # BeSS
  # -------------------
  if (family != "cox") {
    bess.fit <- bess(x[train.index, ], y[train.index], family = family)
    bess.index <- which(bess.fit$beta != 0)
    
  } else {
    bess.fit <- bess(x[train.index, ], y[train.index, ], family = family)
    bess.index <- which(bess.fit$beta != 0)
  }

  # out 1: support recovery
  out.bess.1 <- setequal(true.index, bess.index)

  # out 2: power
  out.bess.2 <- sum(true.index %in% bess.index) / (length(true.index))

  # out 3: type i error rate
  out.bess.3 <- sum(bess.index %in% noise.index) / (p - s)

  # out 4: false discovery rate
  out.bess.4 <- ifelse(length(bess.index) > 0,
    sum(!bess.index %in% true.index) / (length(bess.index)),
    0
  )

  # out 5: false non discovery rate
  out.bess.5 <- ifelse(length(bess.index) < p,
    sum(setdiff(1:p, bess.index) %in% true.index) /
      (p - length(bess.index)), 0
  )

  # out 6: parameter estimation
  out.bess.6 <- mean(abs(bess.fit$beta - true.beta))

  # out 7: prediction performance
  if (family == "binomial") {
    bess.pred <- predict(bess.fit, newx = x[test.index, ], type = "response")
    
    out.bess.7 <- roc(y[test.index] ~ bess.pred, plot = F, print.auc = F)$auc
  } else if (family == "poisson") {
    if (length(bess.index) > 0) {
      d.bess <- data.frame(y = y, x[, bess.index])
      
      poisson.bess <- glm(y ~ ., data = d.bess[train.index, ], family = "poisson", maxit = 3e2)
      
      bess.pred <- predict(poisson.bess, newdata = d.bess[test.index, ], type = "response")
      rmse.bess <- sqrt(mean((bess.pred - y[test.index])^2))
    } else {
      d.train <- data.frame(y = y[train.index], x[train.index, ])
      poisson.bess <- glm(y ~ 1, data = d.train, family = "poisson", maxit = 3e2)
      d.test <- data.frame(y = y[test.index], x[test.index, ])
      bess.pred <- predict(poisson.bess, newdata = d.test, type = "response")
      rmse.bess <- sqrt(mean((bess.pred - y[test.index])^2))
    }
    
    out.bess.7 <- rmse.bess
  } else {
    out.bess.7 <- 0
  }

  # -------------------
  # SIS
  # -------------------
  if (family != "cox") {
    invisible(capture.output(sis.m <- SIS(x[train.index, ],
      y[train.index],
      family = family, tune = "ebic"
    )))
  } else {
    invisible(capture.output(sis.m <- SIS(x[train.index, ],
      Surv(
        y[train.index, 1],
        y[train.index, 2]
      ),
      family = "cox",
      tune = "ebic", penalty = "lasso"
    )))
  }
  sis.index <- sis.m$ix

  # out 1: support recovery
  out.sis.1 <- setequal(true.index, sis.index)

  # out 2: power
  out.sis.2 <- sum(true.index %in% sis.index) / (length(true.index))

  # out 3: type i error rate
  out.sis.3 <- sum(sis.index %in% noise.index) / (p - s)

  # out 4: false discovery rate
  out.sis.4 <- ifelse(length(sis.index) > 0,
    sum(!sis.index %in% true.index) / (length(sis.index)),
    0
  )

  # out 5: false non discovery rate
  out.sis.5 <- ifelse(length(sis.index) < p,
    sum(setdiff(1:p, sis.index) %in% true.index) /
      (p - length(sis.index)), 0
  )

  # out 6: parameter estimation
  sis.coef <- integer(p)

  if (length(sis.index) > 0) {
    if (family != "cox") {
      out.sis.coef <- as.numeric(sis.m$coef.est[-1])
    } else {
      out.sis.coef <- as.numeric(sis.m$coef.est)
    }


    for (i in 1:length(sis.index)) {
      sis.coef[sis.index[i]] <- out.sis.coef[i]
    }
  }

  out.sis.6 <- mean(abs(sis.coef - true.beta))

  # out 7: prediction performance
  if (family == "binomial") {
    sis.pred <- predict(sis.m, newx = x[test.index, ], type = "response")

    out.sis.7 <- roc(y[test.index] ~ sis.pred, plot = F, print.auc = F)$auc
  } else if (family == "poisson") {
    sis.pred <- predict(sis.m, newx = x[test.index, ], type = "response")

    out.sis.7 <- sqrt(mean((sis.pred - y[test.index])^2))
  } else {
    out.sis.7 <- 0
  }


  return(c(
    out.sgpv.1, out.sgpv.2, out.sgpv.3, out.sgpv.4,
    out.sgpv.5, out.sgpv.6, out.sgpv.7,

    out.lasso.1, out.lasso.2, out.lasso.3, out.lasso.4,
    out.lasso.5, out.lasso.6, out.lasso.7,

    out.bess.1, out.bess.2, out.bess.3, out.bess.4,
    out.bess.5, out.bess.6, out.bess.7,

    out.sis.1, out.sis.2, out.sis.3, out.sis.4,
    out.sis.5, out.sis.6, out.sis.7
  ))
}


main.many <- function(num.sim = 1e3, n, p, s, rho = 0.35, sig = 2,
                      beta.min = 0.1, beta.max = 0.4,
                      intercept = 0,
                      family = c("binomial", "poisson", "cox"),
                      scale = 2, shape = 1, rateC = 0.2) {
  suppressMessages(
    out <- replicate(num.sim, main.one.time(
      n = n, p = p, s = s, rho = rho,
      sig = sig, beta.min = beta.min,
      beta.max = beta.max,
      intercept = intercept,
      family = family, scale = scale,
      shape = shape, rateC = rateC
    ))
  )

  return(c(
    mean(out[1, ]), # sgpv capture rate
    mean(out[2, ]), # sgpv power
    mean(out[3, ]), # sgpv type i error rate
    mean(out[4, ]), # sgpv fdr
    mean(out[5, ]), # sgpv fndr
    median(out[6, ]), # sgpv median parameter estimation
    quantile(out[6, ], 0.25), # sgpv parameter estimation first quartile
    quantile(out[6, ], 0.75), # sgpv parameter estimation third quartile
    median(out[7, ]), # sgpv median prediction
    quantile(out[7, ], 0.25), # sgpv prediction first quartile
    quantile(out[7, ], 0.75), # sgpv prediction third quartile

    mean(out[8, ]), # lasso capture rate
    mean(out[9, ]), # lasso power
    mean(out[10, ]), # lasso type i error rate
    mean(out[11, ]), # lasso fdr
    mean(out[12, ]), # lasso fndr
    median(out[13, ]), # lasso median parameter estimation
    quantile(out[13, ], 0.25), # lasso parameter estimation first quartile
    quantile(out[13, ], 0.75), # lasso parameter estimation third quartile
    median(out[14, ]), # lasso median prediction
    quantile(out[14, ], 0.25), # lasso prediction first quartile
    quantile(out[14, ], 0.75), # lasso prediction third quartile

    mean(out[15, ]), # bess capture rate
    mean(out[16, ]), # bess power
    mean(out[17, ]), # bess type i error rate
    mean(out[18, ]), # bess fdr
    mean(out[19, ]), # bess fndr
    median(out[20, ]), # bess median parameter estimation
    quantile(out[20, ], 0.25), # bess parameter estimation first quartile
    quantile(out[20, ], 0.75), # bess parameter estimation third quartile
    median(out[21, ]), # bess median prediction
    quantile(out[21, ], 0.25), # bess prediction first quartile
    quantile(out[21, ], 0.75), # bess prediction third quartile

    mean(out[22, ]), # sis capture rate
    mean(out[23, ]), # sis power
    mean(out[24, ]), # sis type i error rate
    mean(out[25, ]), # sis fdr
    mean(out[26, ]), # sis fndr
    median(out[27, ]), # sis median parameter estimation
    quantile(out[27, ], 0.25), # sis parameter estimation first quartile
    quantile(out[27, ], 0.75), # sis parameter estimation third quartile
    median(out[28, ]), # sis median prediction
    quantile(out[28, ], 0.25), # sis prediction first quartile
    quantile(out[28, ], 0.75) # sis prediction third quartile
  ))
}

# -----------------------
# FIGURE 2
# -----------------------

get.plot.fig.2 <- function(data, p = 200,
                           type = c("lds", "ldd", "hds"), num.sim = 1e3) {
  title.p <- ifelse(type == "lds", "low-d sparse",
    ifelse(type == "ldd", "low-d dense",
      ifelse(type == "hds", "high-d sparse")
    )
  )

  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")

  if (type != "hds") {
    plot.d <- data.frame(
      x = rep(seq(2, 40, 2), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 20),
      rate = c(data[1, ], data[12, ], data[23, ], data[34, ])
    )

    xlim <- c(1, 40)
    xbreaks <- seq(2, 40, 4)
    xlab <- "n/p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(2, 40, 2), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 20)
    )
    
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, p * 4, p / 5), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 16),
      rate = c(data[1, ], data[12, ], data[23, ], data[34, ])
    )

    xlim <- c(p, p * 4)
    xbreaks <- seq(p, p * 4, p / 2)
    xlab <- "p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(p, p * 4, p / 5), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 16)
    )
    
  }


  plot.d$Method <- factor(plot.d$Method,
    levels = c(
      "ProSGPV", "Lasso", "BeSS", "ISIS"
    )
  )

  lb.out <- NULL
  ub.out <- NULL

  for (i in c(1, 12, 23, 34)) {
    pe <- as.numeric(data[i, ])
    lb.out <- c(lb.out, pe - 1.96 * sqrt(pe * (1 - pe) / num.sim))
    ub.out <- c(ub.out, pe + 1.96 * sqrt(pe * (1 - pe) / num.sim))
  }

  ci.d$lb <- lb.out
  ci.d$ub <- ub.out

  ci.d$lb[ci.d$lb < 0] <- 0
  ci.d$ub[ci.d$ub > 1] <- 1

  ci.d$Method <- factor(ci.d$Method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS"
  ))

  ggplot()+
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method))  + 
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = xlab, y = "Support recovery rate", col = "Method",
      title = title.p
    ) + geom_ribbon(data = ci.d, aes(x = x, ymin = lb, ymax = ub, fill = Method), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

# -----------------------
# FIGURE 3
# -----------------------

get.plot.fig.3 <- function(data, p = 200, type = c("lds", "ldd", "hds"), 
                           ycap, ybreaks) {
  title.p <- ifelse(type == "lds", "low-d sparse",
    ifelse(type == "ldd", "low-d dense",
      ifelse(type == "hds", "high-d sparse")
    )
  )

  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")

  if (type != "hds") {
    plot.d <- data.frame(
      x = rep(seq(2, 40, 2), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 20),
      rate = c(data[6, ], data[17, ], data[28, ], data[39, ])
    )

    xlim <- c(1, 40)
    xbreaks <- seq(2, 40, 4)
    xlab <- "n/p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(2, 40, 2), 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 20)
    )
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, p * 4, p / 5), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 16),
      rate = c(data[6, ], data[17, ], data[28, ], data[39, ])
    )

    xlim <- c(p, p * 4)
    xbreaks <- seq(p, p * 4, p / 2)
    xlab <- "p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(p, p * 4, p / 5), 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 16)
    )
  }

  plot.d$rate[plot.d$rate > ycap] <- ycap

  plot.d$Method <- factor(plot.d$Method,
    levels = c("ProSGPV", "Lasso", "BeSS", "ISIS")
  )

  lb.out <- NULL
  ub.out <- NULL

  for (i in c(7, 18, 29, 40)) {
    lb.out <- c(lb.out, as.numeric(data[i, ]))
    ub.out <- c(ub.out, as.numeric(data[i + 1, ]))
  }

  ci.d$lb <- lb.out
  ci.d$ub <- ub.out

  ci.d$method <- factor(ci.d$method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS"
  ))

  ci.d$ub[ci.d$ub > ycap] <- ycap

  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method)) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(
      limits = c(0, ycap),
      breaks = ybreaks
    ) +
    labs(
      x = xlab, y = "Estimation MAE", col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = ci.d, aes(x = x, ymin = lb, ymax = ub, fill = method), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

# -----------------------
# FIGURE 4
# -----------------------

get.plot.fig.4 <- function(data, p = 200, model = c("logistic", "poisson"),
                           ycap, ybreaks,
                           type = c("lds", "ldd", "hds")) {
  title.p <- ifelse(type == "lds", "low-d sparse",
    ifelse(type == "ldd", "low-d dense",
      ifelse(type == "hds", "high-d sparse")
    )
  )

  if (model == "logistic") {
    ylab <- "Prediction AUC"
    ylim <- c(0.5, 1)
    ybreaks <- seq(0.5, 1, 0.1)
  } else {
    ylab <- "Prediction RMSE"
    ylim <- c(0, ycap)
  }

  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")

  if (type != "hds") {
    plot.d <- data.frame(
      x = rep(seq(2, 40, 2), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 20),
      rate = c(data[9, ], data[20, ], data[31, ], data[42, ])
    )

    xlim <- c(1, 40)
    xbreaks <- seq(2, 40, 4)
    xlab <- "n/p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(2, 40, 2), 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 20)
    )
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, p * 4, p / 5), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 16),
      rate = c(data[9, ], data[20, ], data[31, ], data[42, ])
    )

    xlim <- c(p, p * 4)
    xbreaks <- seq(p, p * 4, p / 2)
    xlab <- "p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(p, p * 4, p / 5), 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 16)
    )
  }


  plot.d$Method <- factor(plot.d$Method,
    levels = c("ProSGPV", "Lasso", "BeSS", "ISIS")
  )

  lb.out <- NULL
  ub.out <- NULL

  for (i in c(10, 21, 32, 43)) {
    lb.out <- c(lb.out, as.numeric(data[i, ]))
    ub.out <- c(ub.out, as.numeric(data[i + 1, ]))
  }

  ci.d$lb <- lb.out
  ci.d$ub <- ub.out
  
  if(model=="poisson"){
    plot.d$rate[plot.d$rate > ycap] <- ycap
    ci.d$ub[ci.d$ub > ycap] <- ycap
  }
 

  ci.d$method <- factor(ci.d$method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS"
  ))

  ggplot() +
    geom_line(
      data = plot.d,
      aes(x = x, y = rate, col = Method)
    ) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = ylim, breaks = ybreaks) +
    labs(
      x = xlab, y = ylab, col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = ci.d, aes(
      x = x, ymin = lb, ymax = ub, group = method,
      fill = method
    ), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

# -----------------------
# SUPPLEMENTARY FIGURE 1
# -----------------------

# prosgpv - glm
sgpv.supp.fig.1 <- function(x, y, family = c("binomial", "poisson", "cox")) {
  family <- match.arg(family)

  if (is.null(colnames(x))) colnames(x) <- paste("X", 1:ncol(x), sep = "")

  if (family != "cox") {
    cv.lasso <- cv.glmnet(x, y, family = family)
    lambda.use <- runif(1, 0.8 * cv.lasso$lambda.min, 1.2 * cv.lasso$lambda.1se)
    candidate.index <- which(coef(cv.lasso, s = lambda.use)[-1] != 0)
  } else {
    cv.lasso <- cv.glmnet(x, Surv(y[, 1], y[, 2]), family = "cox")
    lambda.use <- runif(1, 0.8 * cv.lasso$lambda.min, 1.2 * cv.lasso$lambda.1se)
    candidate.index <- which(coef(cv.lasso, s = lambda.use) != 0)
  }

  if (length(candidate.index) > 0) {
    if (family == "binomial") {
      glm.m <- glm(y ~ x[, candidate.index],
        family = family, method = "brglmFit", type = "MPL_Jeffreys"
      )
      pe <- coef(glm.m)[-1]
      se <- summary(glm.m)$coef[-1, 2]
    } else if (family == "poisson") {
      glm.m <- glm(y ~ x[, candidate.index],
        family = family
      )
      pe <- coef(glm.m)[-1]
      se <- summary(glm.m)$coef[-1, 2]
    }

    else {
      cox.train <- coxph(Surv(y[, 1], y[, 2]) ~ x[, candidate.index])
      pe <- coef(cox.train)
      se <- summary(cox.train)$coef[, 3]
    }

    sgpv.index <- candidate.index[as.numeric(which(abs(pe) > 1.96 * se + mean(se)))]
  } else {
    sgpv.index <- NULL
  }

  out <- list(
    var.index = sgpv.index,
    var.label = colnames(x)[sgpv.index],
    x = x,
    y = y,
    family = family
  )

  class(out) <- "sgpv.glm"

  return(out)
}

one.time.supp.fig.1 <- function(n, p, s, rho = 0.35, sig = 2,
                                beta.min = 0.1, beta.max = 0.4,
                                intercept = 2,
                                family = c("binomial", "poisson", "cox"),
                                scale = 2, shape = 1, rateC = 0.2) {
  family <- match.arg(family)

  # simulate training data
  sim.data <- gen.sim.data(
    n = n, p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max,
    intercept = intercept, family = family,
    scale = 2, shape = 1, rateC = 0.2
  )
  x <- sim.data[[1]]
  y <- sim.data[[2]]
  true.index <- sim.data[[3]]
  
  if(family == "cox"){
    try(out.sgpv <- pro.sgpv(x, y, family = family), 
        silent=T)
    try(out.sgpv.unif <- sgpv.supp.fig.1(x, y, family = family), silent=T)
    
    if(inherits(out.sgpv,"try-error")){
      sgpv.index.fix <- numeric(0)
    }else{
      sgpv.index.fix <- out.sgpv$var.index
    }
    
    if(inherits(out.sgpv.unif,"try-error")){
      sgpv.index.unif <- numeric(0)
    }else{
      sgpv.index.unif <- out.sgpv.unif$var.index
    }
    
  }else{
    out.sgpv.fix <- pro.sgpv(x, y, family = family)
    sgpv.index.fix <- out.sgpv.fix$var.index
    
    sgpv.index.unif <- sgpv.supp.fig.1(x, y, family = family)
    sgpv.index.fix <- sgpv.index.unif$var.index
  }  
  
  
  return(c(
    setequal(true.index, sgpv.index.fix),
    setequal(true.index, sgpv.index.unif)
  ))
}


many.sim.supp.fig.1 <- function(num.sim = 1e3, n, p, s, rho = 0.35, sig = 2,
                                beta.min = 0.1, beta.max = 0.4,
                                intercept = 0,
                                family = c("binomial", "poisson", "cox"),
                                scale = 2, shape = 1, rateC = 0.2) {
  out <- replicate(num.sim, one.time.supp.fig.1(
    n = n, p = p, s = s, rho = rho,
    sig = sig, beta.min = beta.min,
    beta.max = beta.max,
    intercept = intercept,
    family = family, scale = scale,
    shape = shape, rateC = rateC
  ))


  return(apply(out, 1, mean))
}


get.plot.supp.fig.1 <- function(data, p = 200,
                                type = c("lds", "ldd", "hds"), num.sim = 1e3) {
  title.p <- ifelse(type == "lds", "low-d sparse",
    ifelse(type == "ldd", "low-d dense",
      ifelse(type == "hds", "high-d sparse")
    )
  )

  # color scheme
  dcols <- c("black", "springgreen3")

  if (type != "hds") {
    plot.d <- data.frame(
      x = rep(seq(2, 40, 2), 2),
      Method = rep(c(
        "Fixed", "Uniform"
      ), each = 20),
      rate = c(data[1, ], data[2, ])
    )

    xlim <- c(1, 40)
    xbreaks <- seq(2, 40, 4)
    xlab <- "n/p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(2, 40, 2), 2),
      method = rep(c(
        "Fixed", "Uniform"
      ), each = 20)
    )
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, p * 4, p / 5), 2),
      Method = rep(c(
        "Fixed", "Uniform"
      ), each = 16),
      rate = c(data[1, ], data[2, ])
    )

    xlim <- c(p, p * 4)
    xbreaks <- seq(p, p * 4, p / 2)
    xlab <- "p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(p, p * 4, p / 5), 2),
      method = rep(c(
        "Fixed", "Uniform"
      ), each = 16)
    )
  }


  plot.d$Method <- factor(plot.d$Method,
    levels = c(
      "Fixed", "Uniform"
    )
  )

  lb.out <- NULL
  ub.out <- NULL

  for (i in c(1, 2)) {
    pe <- as.numeric(data[i, ])
    lb.out <- c(lb.out, pe - 1.96 * sqrt(pe * (1 - pe) / num.sim))
    ub.out <- c(ub.out, pe + 1.96 * sqrt(pe * (1 - pe) / num.sim))
  }

  ci.d$lb <- lb.out
  ci.d$ub <- ub.out

  ci.d$lb[ci.d$lb < 0] <- 0
  ci.d$ub[ci.d$ub > 1] <- 1

  ci.d$method <- factor(ci.d$method, levels = c(
    "Fixed", "Uniform"
  ))

  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method)) +
    scale_color_manual(
      values = dcols,
      labels = c(
        bquote(lambda["gic"]),
        bquote(lambda["cv"])
      )
    ) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = xlab, y = "Support recovery rate", col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = ci.d, aes(x = x, ymin = lb, ymax = ub, fill = method), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


# -----------------------
# SUPPLEMENTARY FIGURE 2
# -----------------------

# get the index of selected variables
get.out.supp.fig.2 <- function(candidate.index, x, y, n, p = 20) {
  out.original <- out.0 <- out.sigma <- out.i.sqrt.log.n <-
    out.sqrt.log.n <- integer(0)

  if (length(candidate.index) > 0) {
    # run fully relaxed LASSO
    f.l <- glm(y ~ x[, candidate.index],
      family = "binomial",
      method = "brglmFit", type = "MPL_Jeffreys"
    )

    # get confidence bands
    pe <- summary(f.l)$coef[-1, 1]
    se <- summary(f.l)$coef[-1, 2]
    null.bound.p <- mean(se)

    null.bound.original <- null.bound.p
    null.bound.sqrt.log.n <- null.bound.p * sqrt(log(n / p))
    null.bound.i.sqrt.log.n <- null.bound.p / sqrt(log(n / p))
    null.bound.sigma <- null.bound.p * sqrt(n / p) /2
    null.bound.0 <- 0

    out.original <- candidate.index[which(abs(pe) > 1.96 * se + null.bound.original)]
    out.sqrt.log.n <- candidate.index[which(abs(pe) > 1.96 * se + null.bound.sqrt.log.n)]
    out.i.sqrt.log.n <- candidate.index[which(abs(pe) > 1.96 * se + null.bound.i.sqrt.log.n)]
    out.sigma <- candidate.index[which(abs(pe) > 1.96 * se + null.bound.sigma)]
    out.0 <- candidate.index[which(abs(pe) > 1.96 * se)]
  } else {
    null.bound.original <- null.bound.0 <-
      null.bound.sqrt.log.n <- null.bound.sigma <- null.bound.i.sqrt.log.n <- NA
  }


  return(list(
    null.bound.original, null.bound.sqrt.log.n,
    null.bound.sigma, null.bound.0, null.bound.i.sqrt.log.n,
    out.original, out.sqrt.log.n, out.sigma, out.0, out.i.sqrt.log.n
  ))
}

# main function
one.time.supp.fig.2 <- function(n, p, s, rho) {

  # generate data
  sim.data <- gen.sim.data(
    n = n, p = p, s = s, rho = rho,
    beta.min = 0.5, beta.max = 1.5,
    sig = 2, intercept = 0, family = "binomial"
  )
  x <- sim.data[[1]]
  y <- sim.data[[2]]
  true.index <- sim.data[[3]]

  lasso.log <- glmnet(x,y,family = "binomial")
  
  yhat.m <- predict(lasso.log,newx=x,type="response")
  p.k <- lasso.log$df + 1
  
  lasso.lambda.index <- which.min(sapply(1:ncol(yhat.m),function(z)
    binomial()$aic(y, rep(1,n), yhat.m[,z], rep(1,n)) +
      p.k[z] * log(log(n))*log(p.k[z]) ) )
  lambda <- lasso.log$lambda[lasso.lambda.index]
  candidate.index <- which(coef(lasso.log,s = lambda)[-1] != 0)

  # original
  temp <- get.out.supp.fig.2(candidate.index, x, y, n)
  bound.original <- temp[[1]]
  index.original <- temp[[6]]

  # root log n
  bound.root.log.n <- temp[[2]]
  index.root.log.n <- temp[[7]]

  # sigma
  bound.sigma <- temp[[3]]
  index.sigma <- temp[[8]]

  # 0
  bound.0 <- temp[[4]]
  index.0 <- temp[[9]]

  # inverse root log n
  bound.i.root.log.n <- temp[[5]]
  index.i.root.log.n <- temp[[10]]


  return(c(
    bound.original,
    bound.root.log.n,
    bound.sigma,
    bound.0,
    bound.i.root.log.n,
    setequal(true.index, index.original),
    setequal(true.index, index.root.log.n),
    setequal(true.index, index.sigma),
    setequal(true.index, index.0),
    setequal(true.index, index.i.root.log.n)
  ))
}

many.sim.supp.fig.2 <- function(num.sim = 1e3, n, p = 20, s = 4,
                                rho) {
  out <- NULL
  out <- replicate(num.sim, one.time.supp.fig.2(
    n = n, p = p, s = s, rho = rho
  ))

  return(c(
    median(out[1, ], na.rm = T), # median bound original
    quantile(out[1, ], 0.25, na.rm = T), # first quartile bound original
    quantile(out[1, ], 0.75, na.rm = T), # third quartile bound original

    median(out[2, ], na.rm = T), # median bound root log n
    quantile(out[2, ], 0.25, na.rm = T), # first quartile bound root log n
    quantile(out[2, ], 0.75, na.rm = T), # third quartile bound root log n

    median(out[3, ], na.rm = T), # median bound sigma
    quantile(out[3, ], 0.25, na.rm = T), # first quartile bound sigma
    quantile(out[3, ], 0.75, na.rm = T), # third quartile bound sigma

    median(out[4, ], na.rm = T), # median bound 0
    quantile(out[4, ], 0.25, na.rm = T), # first quartile bound 0
    quantile(out[4, ], 0.75, na.rm = T), # third quartile bound 0

    median(out[5, ], na.rm = T), # median bound inverse root log n
    quantile(out[5, ], 0.25, na.rm = T), # first quartile bound i. root log n
    quantile(out[5, ], 0.75, na.rm = T), # third quartile bound i. root log n

    mean(out[6, ]), # capture rate original
    mean(out[7, ]), # capture rate root log n
    mean(out[8, ]), # capture rate sigma
    mean(out[9, ]), # capture rate 0
    mean(out[10, ]) # capture rate inverse root log n
  ))
}


get.plot.supp.fig.2.bound <- function(data, cor) {

  # get keywords from data names
  title.p <- ifelse(cor == "ind", "Independent",
    ifelse(cor == "medium", "Medium correlation",
      ifelse(cor == "high", "High correlation", "Weird")
    )
  )


  # color scheme
  dcols <- c("black", "springgreen3", "blue", "yellow4", "red")

  plot.d <- data.frame(
    n = rep(seq(2,40,2), each = 5),
    Method = rep(c(
      "Original null bound", "Null bound root log n", "Null bound sigma",
      "0 null bound", "Null bound inverse"
    ), 20),
    rate = c(data[c(1, 4, 7, 10, 13), ])
  )

  plot.d$Method <- factor(plot.d$Method,
    levels = c(
      "Original null bound", "Null bound root log n", "Null bound inverse",
      "Null bound sigma",
      "0 null bound"
    )
  )

  # add confidence interval
  ci.d <- data.frame(
    n = rep(seq(2,40,2), each = 5),
    method = rep(c(
      "Original null bound", "Null bound root log n", "Null bound sigma",
      "0 null bound", "Null bound inverse"
    ), 20),
    lb = c(data[c(2, 5, 8, 11, 14), ]),
    ub = c(data[c(3, 6, 9, 12, 15), ])
  )

  ci.d$method <- factor(ci.d$method,
    levels = c(
      "Original null bound",  "Null bound root log n", "Null bound inverse",
      "Null bound sigma",
      "0 null bound"
    )
  )

  ci.d$ub[ci.d$ub > 0.5] <- 0.5
  plot.d$rate[plot.d$rate > 0.5] <- 0.5
  
  ggplot() +
    geom_line(data = plot.d, aes(x = n, y = rate, col = Method)) +
    scale_color_manual(
      values = dcols,
      labels = c(
        "Original null bound",
        bquote("Null bound *" ~ sqrt(log("n/p"))),
        bquote("Null bound /" ~ sqrt(log("n/p"))),
        bquote("Null bound *" ~ sqrt("n/p") ~ "/2"),
        "0 null bound"
      )
    ) +
    geom_ribbon(data = ci.d, aes(
      x = n, ymin = lb, ymax = ub,
      fill = method
    ), alpha = .2) +
    scale_fill_manual(values = dcols) +
    scale_x_continuous(limits = c(1, 40), breaks = seq(2, 40, 4)) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    labs(
      x = "n/p", y = "Null bound size", col = "Method",
      title = title.p
    ) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


get.plot.supp.fig.2.capture <- function(data, cor, num.sim = 1e3) {

  # get keywords from data names
  title.p <- ifelse(cor == "ind", "Independent",
    ifelse(cor == "medium", "Medium correlation",
      ifelse(cor == "high", "High correlation", "Weird")
    )
  )


  # color scheme
  dcols <- c("black", "springgreen3", "blue", "yellow4", "red")

  plot.d <- data.frame(
    n = rep(seq(2,40,2), each = 5),
    Method = rep(c(
      "Original null bound", "Null bound root log n", "Null bound sigma",
      "0 null bound", "Null bound inverse"
    ), 20),
    rate = c(data[16:20, ])
  )

  plot.d$Method <- factor(plot.d$Method,
    levels = c(
      "Original null bound",  "Null bound root log n", "Null bound inverse",
      "Null bound sigma",
      "0 null bound"
    )
  )

  # add confidence interval
  ci.d <- data.frame(
    n = rep(seq(2,40,2), 5),
    method = rep(c(
      "Original null bound", "Null bound root log n", "Null bound sigma",
      "0 null bound", "Null bound inverse"
    ), each = 20)
  )

  lb.out <- NULL
  ub.out <- NULL

  for (i in 16:20) {
    pe <- as.numeric(data[i, ])
    lb.out <- c(lb.out, pe - 1.96 * sqrt(pe * (1 - pe) / num.sim))
    ub.out <- c(ub.out, pe + 1.96 * sqrt(pe * (1 - pe) / num.sim))
  }

  ci.d$lb <- lb.out
  ci.d$ub <- ub.out

  ci.d[ci.d < 0] <- 0

  ci.d$method <- factor(ci.d$method,
    levels = c(
      "Original null bound",  "Null bound root log n", "Null bound inverse",
      "Null bound sigma",
      "0 null bound"
    )
  )

  ggplot() +
    geom_line(data = plot.d, aes(x = n, y = rate, col = Method)) +
    scale_color_manual(
      values = dcols,
      labels = c(
        "Original null bound",
        bquote("Null bound *" ~ sqrt(log("n/p"))),
        bquote("Null bound /" ~ sqrt(log("n/p"))),
        bquote("Null bound *" ~ sqrt("n/p")),
        "0 null bound"
      )
    ) +
    geom_ribbon(data = ci.d, aes(
      x = n, ymin = lb, ymax = ub,
      fill = method
    ), alpha = .2) +
    scale_fill_manual(values = dcols) +
    scale_x_continuous(limits = c(1, 40), breaks = seq(2, 40, 4)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = "n/p", y = "Support recovery rate", col = "Method",
      title = title.p
    ) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


# ----------------------
# SUPPLEMENTARY FIGURE 3
# ----------------------

# prosgpv - glm
pro.sgpv.glm <- function(x, y, family = "binomial") {
  
  n <- length(y)
  
  lasso.log <- glmnet(x,y,family = "binomial")
  
  yhat.m <- predict(lasso.log,newx=x,type="response")
  p.k <- lasso.log$df + 1
  
  lasso.lambda.index <- which.min(sapply(1:ncol(yhat.m),function(z)
    binomial()$aic(y, rep(1,n), yhat.m[,z], rep(1,n)) +
      p.k[z] * log(log(n))*log(p.k[z]) ) )
  
  lambda <- lasso.log$lambda[lasso.lambda.index]
  candidate.index <- which(coef(lasso.log,s = lambda)[-1] != 0)
  
  if (is.null(colnames(x))) colnames(x) <- paste("X", 1:ncol(x), sep = "")
  
  if (length(candidate.index) > 0) {
    glm.m <- glm(y ~ x[, candidate.index],
                 family = family, maxit = 3e2
    )
    pe <- coef(glm.m)[-1]
    se <- summary(glm.m)$coef[-1, 2]
    
    sgpv.index <- candidate.index[as.numeric(which(abs(pe) > 1.96 * se + mean(se)))]
  } else {
    sgpv.index <- NULL
  }
  
  out <- list(
    var.index = sgpv.index,
    var.label = colnames(x)[sgpv.index],
    x = data.frame(x),
    y = y,
    family = family
  )
  
  class(out) <- "sgpv.glm"
  
  return(out)
}

coef.sgpv.glm <- function(object, ...) {
  out.coef <- numeric(ncol(object$x))
  if (length(object$var.index) > 0) {
    glm.d <- data.frame(yy = object$y, xx = object$x[, object$var.index])
    
    out.sgpv.coef <- coef(glm(yy ~ ., family = object$family, data = glm.d))[-1]
    
    for (i in 1:length(object$var.index)) {
      out.coef[object$var.index[i]] <- out.sgpv.coef[i]
    }
  }
  
  out.coef
}

predict.sgpv.glm <- function(object, newdata, ...) {
  if (length(object$var.index) > 0) {
    glm.d <- data.frame(yy = object$y, xx = object$x[, object$var.index])
    colnames(glm.d)[-1] <- object$var.label
    glm.m <- glm(yy ~ ., family = object$family, data = glm.d)
    if (missing(newdata)) newdata <- object$x
    
    predict(glm.m, data.frame(newdata), type = "response")
  } else {
    glm.m <- glm(yy ~ 1,
                 data = data.frame(
                   yy = object$y,
                   xx = object$x
                 ),
                 family = object$family
    )
    if (missing(newdata)) newdata <- object$x
    predict(glm.m, data.frame(newdata), type = "response")
  }
}

supp.fig.3.one.time <- function(n, p, s, rho = 0.35, sig = 2,
                                beta.min = 0.5, beta.max = 1.5,
                                intercept = 0, family = "binomial") {
  
  # simulate training data
  sim.data <- gen.sim.data(
    n = round(n * 5 / 3), p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max,
    intercept = intercept, family = family,
    scale = 2, shape = 1, rateC = 0.2
  )
  x <- sim.data[[1]]
  y <- sim.data[[2]]
  true.index <- sim.data[[3]]
  true.beta <- sim.data[[4]]
  
  train.index <- sample(1:round(n * 5 / 3), n, replace = F)
  test.index <- setdiff(1:round(n * 5 / 3), train.index)
  
  # prosgpv with ML
  out.original <- pro.sgpv.glm(x[train.index, ], y[train.index], family = "binomial"
  )
  original.index <- out.original$var.index
  original.pe <- mean(abs(coef(out.original) - true.beta))
  original.pred <- roc(y[test.index] ~ predict(out.original, newdata = x[test.index, ]),
                       plot = F, print.auc = F
  )$auc
  
  # jeffreys prior
  out.jeffreys <- pro.sgpv(x[train.index, ], y[train.index], family = "binomial"
  )
  jeffreys.index <- out.jeffreys$var.index
  jeffreys.pe <- mean(abs(coef(out.jeffreys) - true.beta))
  jeffreys.pred <- roc(y[test.index] ~ predict(out.jeffreys, newdata = x[test.index, ]),
                       plot = F, print.auc = F
  )$auc
  
  return(c(
    setequal(true.index, original.index),
    original.pe,
    original.pred,
    
    setequal(true.index, jeffreys.index),
    jeffreys.pe,
    jeffreys.pred
  ))
}


supp.fig.3.many <- function(num.sim = 1e3, n = 100, p = 50, s = 4,
                            rho = 0.35, sig = 2,
                            beta.min = 0.5, beta.max = 1.5,
                            intercept = 0,
                            family = "binomial",
                            scale = 2, shape = 1, rateC = 0.2) {
  suppressMessages(out <- replicate(num.sim, supp.fig.3.one.time(
    n = n, p = p, s = s, rho = rho,
    sig = sig, beta.min = beta.min,
    beta.max = beta.max,
    intercept = intercept,
    family = family
  )))
  
  return(c(
    mean(out[1, ]), # capture rate of the true model
    median(out[2, ]), # parameter estimation
    as.numeric(quantile(out[2, ], 0.25)),
    as.numeric(quantile(out[2, ], 0.75)),
    median(out[3, ]), # prediction
    as.numeric(quantile(out[3, ], 0.25)),
    as.numeric(quantile(out[3, ], 0.75)),
    
    mean(out[4, ]), # capture rate of the true model
    median(out[5, ]), # parameter estimation
    as.numeric(quantile(out[5, ], 0.25)),
    as.numeric(quantile(out[5, ], 0.75)),
    median(out[6, ]), # prediction
    as.numeric(quantile(out[6, ], 0.25)),
    as.numeric(quantile(out[6, ], 0.75))
  ))
}


get.plot.supp.fig.3 <- function(data, outcome = c("rate", "pe", "pr"), 
                                x.n = c(seq(2,40,2), seq(200,800,40)), 
                                xbreaks = c(seq(2,40,4), seq(200,800,100)),
                                pe.flag = c("ld.s", "ld.d", "hd.s"),
                                num.sim = 1e3) {
  
  # get keywords from data names
  title.p <- ifelse(outcome == "rate", "Support recovery",
                    ifelse(outcome == "pe", "Parameter estimation",
                           ifelse(outcome == "pr", "Prediction", "Weird")
                    )
  )
  
  dcols <- c("black", "blue")
  
  if (outcome == "rate") {
    plot.d <- data.frame(
      n = rep(x.n , each = 2),
      Method = rep(c(
        "Maximum likelihood", "Jeffreys prior"
      ), length(x.n)),
      rate = c(data[c(1, 8), ])
    )
    
    ylim <- c(0, 1)
    ybreaks <- seq(0, 1, 0.2)
    
    ci.d <- data.frame(
      x = rep(x.n, 2),
      method = rep(c(
        "Maximum likelihood", "Jeffreys prior"
      ), each = length(x.n))
    )
    
    lb.out <- NULL
    ub.out <- NULL
    
    for (i in c(1, 8)) {
      pe <- as.numeric(data[i, ])
      lb.out <- c(lb.out, pe - 1.96 * sqrt(pe * (1 - pe) / num.sim))
      ub.out <- c(ub.out, pe + 1.96 * sqrt(pe * (1 - pe) / num.sim))
    }
    
    ci.d$lb <- lb.out
    ci.d$ub <- ub.out
    
    ci.d$lb[ci.d$lb < 0] <- 0
    ci.d$ub[ci.d$ub > 1] <- 1
    
  } else if (outcome == "pe") {
    plot.d <- data.frame(
      x = rep(x.n, each = 2),
      Method = rep(c(
        "Maximum likelihood", "Jeffreys prior"
      ), length(x.n)),
      pe = c(data[c(2, 9), ]),
      lb = c(data[c(3, 10), ]),
      ub = c(data[c(4, 11), ])
    )
    
    ylab <- "MAE"
    
    if(pe.flag == "ld.s"){
      ylim <- c(0, 0.2)
      ybreaks <- seq(0, 0.2, 0.05)
    }else if(pe.flag == "ld.d"){
      ylim <- c(0, 0.7)
      ybreaks <- seq(0, 0.7, 0.1)
    }else{
      ylim <- c(0, 0.006)
      ybreaks <- seq(0, 0.006, 0.001)
    }
    
  } else if (outcome == "pr") {
    plot.d <- data.frame(
      x = rep(x.n, each = 2),
      Method = rep(c(
        "Maximum likelihood", "Jeffreys prior"
      ), length(x.n)),
      pe = c(data[c(5, 12), ]),
      lb = c(data[c(6, 13), ]),
      ub = c(data[c(7, 14), ])
    )
    
    ylab <- "AUC"
    ylim <- c(0, 1)
    ybreaks <- seq(0, 1, 0.2)
  }
  
  if (outcome == "rate") {
    ggplot() +
      geom_line(data = plot.d, aes(x = n, y = rate, col = Method)) +
      scale_color_manual(values = dcols) +
      scale_x_continuous(limits = range(x.n), breaks = xbreaks) +
      scale_y_continuous(limits = ylim, breaks = ybreaks) +
      labs(
        x = "p", y = "Average capture rate", col = "Method",
        title = title.p
      ) + 
      geom_ribbon(data = ci.d, aes(
        x = x, ymin = lb, ymax = ub,
        fill = method
      ), alpha = .3) +
      scale_fill_manual(values = dcols) +
      guides(fill = F) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      )
  } else if (outcome == "pe") {
    ggplot() +
      geom_line(data = plot.d, aes(x = x, y = pe, col = Method)) +
      scale_color_manual(values = dcols) +
      scale_x_continuous(limits = range(x.n), breaks = xbreaks) +
      scale_y_continuous(limits = ylim, breaks = ybreaks) +
      labs(
        x = "p", y = ylab, col = "Method",
        title = title.p
      ) +
      geom_ribbon(data = plot.d, aes(
        x = x, ymin = lb,
        ymax = ub, fill = Method
      ), alpha = .3) +
      scale_fill_manual(values = dcols) +
      guides(fill = F) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      )
  } else if (outcome == "pr") {
    ggplot() +
      geom_line(data = plot.d, aes(x = x, y = pe, col = Method)) +
      scale_color_manual(values = dcols) +
      scale_x_continuous(limits = range(x.n), breaks = xbreaks) +
      scale_y_continuous(limits = ylim, breaks = ybreaks) +
      labs(
        x = "p", y = ylab, col = "Method",
        title = title.p
      ) +
      geom_ribbon(data = plot.d, aes(
        x = x, ymin = lb,
        ymax = ub, fill = Method
      ), alpha = .3) +
      scale_fill_manual(values = dcols) +
      guides(fill = F) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      )
  }
}


# -----------------------
# SUPPLEMENTARY FIGURE 4
# -----------------------


get.plot.supp.fig.4 <- function(t1.data, power.data, p = 200,
                                type = c("lds", "ldd", "hds")) {
  title.p <- ifelse(type == "lds", "low-d sparse",
    ifelse(type == "ldd", "low-d dense",
      ifelse(type == "hds", "high-d sparse")
    )
  )

  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")

  if (type != "hds") {
    plot.d <- data.frame(
      x = rep(rep(seq(2, 40, 2), each = 4), 2),
      Method = rep(rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), 20), 2),
      rate = c(c(t1.data), c(power.data)),
      Type = rep(c("Type I Error rate", "Power"), each = 80)
    )
    xlab <- "n/p"
    xlim <- c(1, 40)
    xbreaks <- seq(2, 40, 4)
  } else {
    plot.d <- data.frame(
      x = rep(rep(seq(p, 4 * p, p / 5), each = 4), 2),
      Method = rep(rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), 16), 2),
      rate = c(c(t1.data), c(power.data)),
      Type = rep(c("Type I Error rate", "Power"), each = 64)
    )
    xlab <- "p"
    xlim <- c(p, 4 * p)
    xbreaks <- seq(p, 4 * p, p / 2)
  }


  plot.d$Method <- factor(plot.d$Method,
    levels = c(
      "ProSGPV", "Lasso", "BeSS", "ISIS"
    )
  )

  plot.d$Type <- factor(plot.d$Type, levels = c("Power", "Type I Error rate"))

  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method, linetype = Type)) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(
      limits = xlim,
      breaks = xbreaks
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = xlab, y = "Average rate", col = "Method",
      title = title.p
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

# -----------------------
# SUPPLEMENTARY FIGURE 5
# -----------------------

get.plot.supp.fig.5 <- function(fdr.data, fndr.data, p = 200,
                                type = c("lds", "ldd", "hds")) {
  title.p <- ifelse(type == "lds", "low-d sparse",
    ifelse(type == "ldd", "low-d dense",
      ifelse(type == "hds", "high-d sparse")
    )
  )

  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")

  if (type != "hds") {
    plot.d <- data.frame(
      x = rep(rep(seq(2, 40, 2), each = 4), 2),
      Method = rep(rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), 20), 2),
      rate = c(c(fdr.data), c(fndr.data)),
      Type = rep(c("FDR", "FNDR"), each = 80)
    )
    xlab <- "n/p"
    xlim <- c(1, 40)
    xbreaks <- seq(2, 40, 4)
  } else {
    plot.d <- data.frame(
      x = rep(rep(seq(p, 4 * p, p / 5), each = 4), 2),
      Method = rep(rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), 16), 2),
      rate = c(c(fdr.data), c(fndr.data)),
      Type = rep(c("FDR", "FNDR"), each = 64)
    )
    xlab <- "p"
    xlim <- c(p, 4 * p)
    xbreaks <- seq(p, 4 * p, p / 2)
  }


  plot.d$Method <- factor(plot.d$Method,
    levels = c(
      "ProSGPV", "Lasso", "BeSS", "ISIS"
    )
  )

  plot.d$Type <- factor(plot.d$Type, levels = c("FDR", "FNDR"))

  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method, linetype = Type)) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(
      limits = xlim,
      breaks = xbreaks
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = xlab, y = "Average rate", col = "Method",
      title = title.p
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

# ---------------------------------------------
# SUPPLEMENTARY FIGURE 6: LOW-D DENSE LOGISTIC
# ---------------------------------------------
one.time.supp.fig.6 <- function(n = 40, p = 20, s = 14, rho = 0.35, 
                                beta.min = 0.5, beta.max = 1.5, family = "binomial",
                                sig = 2, intercept = 0){
  
  # simulate training data
  sim.data <- gen.sim.data(
    n = round(n * 5 / 3), p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max,
    intercept = intercept, family = family
  )
  
  x <- sim.data[[1]]
  y <- sim.data[[2]]
  true.index <- sim.data[[3]]
  true.beta <- sim.data[[4]]
  
  # generate training and testing indices
  train.index <- sample(1:round(n * 5 / 3), n, replace = F)
  test.index <- setdiff(1:round(n * 5 / 3), train.index)
  
  # -------------------
  # lasso
  # -------------------
  
  cv.m <- cv.glmnet(x[train.index, ], y[train.index], family = family)
  lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)[-1]
  
  lasso.min.index <- which(lasso.min.coef != 0)
  
  # out 1: support recovery
  out.lasso.1 <- setequal(true.index, lasso.min.index)
  
  # out 6: parameter estimation
  out.lasso.6 <- mean(abs(as.numeric(lasso.min.coef) - true.beta))
  
  # out 7: prediction performance
  lasso.pred <- predict(cv.m, s = cv.m$lambda.min, newx = x[test.index, ],
                        type = "response")
  
  out.lasso.7 <- roc(y[test.index] ~ as.vector(lasso.pred), plot = F, print.auc = F)$auc
  
  # -------------------
  # pro.sgpv
  # -------------------
  
  out.sgpv <- pro.sgpv(x[train.index, ], y[train.index], family = family)
  sgpv.index <- out.sgpv$var.index
  
  # out 1: support recovery
  out.sgpv.1 <- setequal(true.index, sgpv.index)
  
  # out 6: parameter estimation
  out.sgpv.6 <- mean(abs(coef(out.sgpv) - true.beta))
  
  
  # out 7: prediction performance
  sgpv.pred <- predict(out.sgpv,
                       newdata = x[test.index, ],
                       type = "response"
  )
  out.sgpv.7 <- roc(y[test.index] ~ sgpv.pred, plot = F, print.auc = F)$auc
  
  # -------------------
  # BeSS
  # -------------------
  
  bess.fit <- bess(x[train.index, ], y[train.index], family = family)
  bess.index <- which(bess.fit$beta != 0)
  
  # out 1: support recovery
  out.bess.1 <- setequal(true.index, bess.index)
  
  # out 6: parameter estimation
  out.bess.6 <- mean(abs(bess.fit$beta - true.beta))
  
  # out 7: prediction performance
  bess.pred <- predict(bess.fit, newx = x[test.index, ], type = "response")
  out.bess.7 <- roc(y[test.index] ~ bess.pred, plot = F, print.auc = F)$auc
  
  # -------------------
  # SIS
  # -------------------
  invisible(capture.output(sis.m <- SIS(x[train.index, ],
                                        y[train.index],
                                        family = family, tune = "ebic"
  )))
  
  sis.index <- sis.m$ix
  
  # out 1: support recovery
  out.sis.1 <- setequal(true.index, sis.index)
  
  # out 6: parameter estimation
  sis.coef <- integer(p)
  
  if (length(sis.index) > 0) {
    out.sis.coef <- as.numeric(sis.m$coef.est[-1])
    
    for (i in 1:length(sis.index)) {
      sis.coef[sis.index[i]] <- out.sis.coef[i]
    }
  }
  
  out.sis.6 <- mean(abs(sis.coef - true.beta))
  
  # out 7: prediction performance
  sis.pred <- predict(sis.m, newx = x[test.index, ], type = "response")
  
  out.sis.7 <- roc(y[test.index] ~ sis.pred, plot = F, print.auc = F)$auc
  
  
  # -------------------
  # pro.sgpv.gvif
  # -------------------
  
  out.sgpv.gvif <- pro.sgpv(x[train.index, ], y[train.index], family = family,
                            gvif = T)
  sgpv.gvif.index <- out.sgpv.gvif$var.index
  
  # out 1: support recovery
  out.sgpv.gvif.1 <- setequal(true.index, sgpv.gvif.index)
  
  # out 6: parameter estimation
  out.sgpv.gvif.6 <- mean(abs(coef(out.sgpv.gvif) - true.beta))
  
  # out 7: prediction performance
  sgpv.gvif.pred <- predict(out.sgpv.gvif,
                            newdata = x[test.index, ],
                            type = "response"
  )
  out.sgpv.gvif.7 <- roc(y[test.index] ~ sgpv.gvif.pred, plot = F, print.auc = F)$auc
  
  
  return(c(out.sgpv.1, 
           out.sgpv.6, 
           out.sgpv.7, 
           
           out.lasso.1, 
           out.lasso.6, 
           out.lasso.7, 
           
           out.bess.1, 
           out.bess.6, 
           out.bess.7, 
           
           out.sis.1,
           out.sis.6,
           out.sis.7,
           
           out.sgpv.gvif.1, 
           out.sgpv.gvif.6, 
           out.sgpv.gvif.7
           
  ))
  
}

many.sim.supp.fig.6 <- function(num.sim = 1e3, n = 40, rho = 0.35, 
                                beta.min = 0.5, beta.max = 1.5, 
                                sig = 2, intercept = 0, p = 20, s = 14){
  
  suppressMessages(
    out <- replicate(num.sim, one.time.supp.fig.6(
      n = n, rho = rho, sig = sig, p = p , s = s,
      beta.min = beta.min, beta.max = beta.max,
      intercept = intercept
    ))
  )
  
  return(c(
    mean(out[1, ]), # sgpv capture rate
    median(out[2, ]), # sgpv median parameter estimation
    quantile(out[2, ], 0.25), # sgpv parameter estimation first quartile
    quantile(out[2, ], 0.75), # sgpv parameter estimation third quartile
    median(out[3, ]), # sgpv median prediction
    quantile(out[3, ], 0.25), # sgpv prediction first quartile
    quantile(out[3, ], 0.75), # sgpv prediction third quartile
    
    mean(out[4, ]), # lasso capture rate
    median(out[5, ]), # lasso median parameter estimation
    quantile(out[5, ], 0.25), # lasso parameter estimation first quartile
    quantile(out[5, ], 0.75), # lasso parameter estimation third quartile
    median(out[6, ]), # lasso median prediction
    quantile(out[6, ], 0.25), # lasso prediction first quartile
    quantile(out[6, ], 0.75), # lasso prediction third quartile
    
    mean(out[7, ]), # bess capture rate
    median(out[8, ]), # bess median parameter estimation
    quantile(out[8, ], 0.25), # bess parameter estimation first quartile
    quantile(out[8, ], 0.75), # bess parameter estimation third quartile
    median(out[9, ]), # bess median prediction
    quantile(out[9, ], 0.25), # bess prediction first quartile
    quantile(out[9, ], 0.75), # bess prediction third quartile
    
    mean(out[10, ]), # sis capture rate
    median(out[11, ]), # sis median parameter estimation
    quantile(out[11, ], 0.25), # sis parameter estimation first quartile
    quantile(out[11, ], 0.75), # sis parameter estimation third quartile
    median(out[12, ]), # sis median prediction
    quantile(out[12, ], 0.25), # sis prediction first quartile
    quantile(out[12, ], 0.75), # sis prediction third quartile
    
    mean(out[13, ]), # sgpv.gvif capture rate
    median(out[14, ]), # sgpv.gvif median parameter estimation
    quantile(out[14, ], 0.25), # sgpv.gvif parameter estimation first quartile
    quantile(out[14, ], 0.75), # sgpv.gvif parameter estimation third quartile
    median(out[15, ]), # sgpv.gvif median prediction
    quantile(out[15, ], 0.25), # sgpv.gvif prediction first quartile
    quantile(out[15, ], 0.75) # sgpv.gvif prediction third quartile
    
  ))
}

get.plot.supp.fig.6.cap <- function(data, num.sim = 1e3, title.p = "Support recovery") {
  
  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red", "darkgoldenrod")
  
  plot.d <- data.frame(
    x = rep(seq(2, 40, 2), 5),
    Method = rep(c(
      "ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
    ), each = 20),
    rate = c(data[1, ], data[8, ], data[15, ], data[22, ], data[29,])
  )
  
  xlim <- c(1, 40)
  xbreaks <- seq(2, 40, 4)
  xlab <- "n/p"
  
  # add confidence interval
  ci.d <- data.frame(
    x = rep(seq(2, 40, 2), 5),
    Method = rep(c(
      "ProSGPV","Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
    ), each = 20)
  )
  
  plot.d$Method <- factor(plot.d$Method,
                          levels = c(
                            "ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
                          )
  )
  
  lb.out <- NULL
  ub.out <- NULL
  
  for (i in c(1,8,15,22,29)) {
    pe <- as.numeric(data[i, ])
    lb.out <- c(lb.out, pe - 1.96 * sqrt(pe * (1 - pe) / num.sim))
    ub.out <- c(ub.out, pe + 1.96 * sqrt(pe * (1 - pe) / num.sim))
  }
  
  ci.d$lb <- lb.out
  ci.d$ub <- ub.out
  
  ci.d$lb[ci.d$lb < 0] <- 0
  ci.d$ub[ci.d$ub > 1] <- 1
  
  ci.d$Method <- factor(ci.d$Method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
  ))
  
  ggplot()+
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method))  + 
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = xlab, y = "Support recovery rate", col = "Method",
      title = title.p
    ) + geom_ribbon(data = ci.d, aes(x = x, ymin = lb, ymax = ub, fill = Method), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

get.plot.supp.fig.6.pe <- function(data, title.p = "Parameter estimation", 
                                   ycap, ybreaks){
  
  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red", "darkgoldenrod")
  
  plot.d <- data.frame(
    x = rep(seq(2, 40, 2), 5),
    Method = rep(c(
      "ProSGPV", "Lasso", "BeSS", "ISIS",  "ProSGPV-GVIF"
    ), each = 20),
    rate = c(data[2, ], data[9, ], data[16, ], data[23, ], data[30,])
  )
  
  xlim <- c(1, 40)
  xbreaks <- seq(2, 40, 4)
  xlab <- "n/p"
  
  # add confidence interval
  ci.d <- data.frame(
    x = rep(seq(2, 40, 2), 5),
    method = rep(c(
      "ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
    ), each = 20)
  )
  
  
  plot.d$rate[plot.d$rate > ycap] <- ycap
  
  plot.d$Method <- factor(plot.d$Method,
                          levels = c("ProSGPV", "Lasso", "BeSS", "ISIS","ProSGPV-GVIF")
  )
  
  lb.out <- NULL
  ub.out <- NULL
  
  for (i in c(3, 10, 17, 24, 31)) {
    lb.out <- c(lb.out, as.numeric(data[i, ]))
    ub.out <- c(ub.out, as.numeric(data[i + 1, ]))
  }
  
  ci.d$lb <- lb.out
  ci.d$ub <- ub.out
  
  ci.d$method <- factor(ci.d$method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
  ))
  
  ci.d$ub[ci.d$ub > ycap] <- ycap
  
  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method)) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(
      limits = c(0, ycap),
      breaks = ybreaks
    ) +
    labs(
      x = xlab, y = "Estimation MAE", col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = ci.d, aes(x = x, ymin = lb, ymax = ub, fill = method), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


get.plot.supp.fig.6.pr <- function(data,title.p = "Prediction") {
  
  
  ylab <- "Prediction AUC"
  ylim <- c(0.5, 1)
  ybreaks <- seq(0.5, 1, 0.1)
  
  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red", "darkgoldenrod")
  
  plot.d <- data.frame(
    x = rep(seq(2, 40, 2), 5),
    Method = rep(c(
      "ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
    ), each = 20),
    rate = c(data[5, ], data[12, ], data[19, ], data[26, ], data[33,])
  )
  
  xlim <- c(1, 40)
  xbreaks <- seq(2, 40, 4)
  xlab <- "n/p"
  
  # add confidence interval
  ci.d <- data.frame(
    x = rep(seq(2, 40, 2), 5),
    method = rep(c(
      "ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
    ), each = 20)
  )
  
  plot.d$Method <- factor(plot.d$Method,
                          levels = c("ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF")
  )
  
  lb.out <- NULL
  ub.out <- NULL
  
  for (i in c(6, 13, 20, 27, 34)) {
    lb.out <- c(lb.out, as.numeric(data[i, ]))
    ub.out <- c(ub.out, as.numeric(data[i + 1, ]))
  }
  
  ci.d$lb <- lb.out
  ci.d$ub <- ub.out
  
  
  ci.d$method <- factor(ci.d$method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS", "ProSGPV-GVIF"
  ))
  
  ggplot() +
    geom_line(
      data = plot.d,
      aes(x = x, y = rate, col = Method)
    ) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = ylim, breaks = ybreaks) +
    labs(
      x = xlab, y = ylab, col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = ci.d, aes(
      x = x, ymin = lb, ymax = ub, group = method,
      fill = method
    ), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


# -----------------------------------------
# SUPPLEMENTARY FIGURE 7: COMPUTATION TIME
# -----------------------------------------

run.1.sgpv <- function(x, y, family) {
  
  sgpv.index <- try(pro.sgpv(x, y, family = family)$var.index, silent=T)
  
}

run.1.lasso <- function(x, y, family) {
  
  if (family != "cox") {
    cv.m <- cv.glmnet(x, y, family = family)
    lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)[-1]
  } else {
    cv.m <- cv.glmnet(x,
                      Surv(y[, 1], y[, 2]),
                      family = "cox"
    )
    lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)
  }
  lasso.min.index <- which(lasso.min.coef != 0)
  
}

run.1.bess <- function(x, y, family) {
  
    bess.fit <- bess(x, y, family = family)
    bess.index <- which(bess.fit$beta != 0)

}

run.1.isis <- function(x, y, family ) {
  
  if (family != "cox") {
    invisible(capture.output(
      sis.m <- SIS(x,  y, family = family, tune = "ebic"
      )))
  } else {
    invisible(capture.output(
      sis.m <- SIS(x, Surv( y[, 1],  y[, 2] ),  
                   family = "cox", tune = "ebic", penalty = "lasso"
      )))
  }
  sis.index <- sis.m$ix
  
}

# main function
one.time.supp.fig.7 <- function(n, p, s, rho,  sig = 2,
                                beta.min = 0.1, beta.max = 0.4,
                                intercept = 2,
                                family = c("binomial", "poisson", "cox"),
                                scale = 2, shape = 1, rateC = 0.2 ) {
  
  # generate data
  sim.data <- gen.sim.data(
    n = n, p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max,
    intercept = intercept, family = family,
    scale = 2, shape = 1, rateC = 0.2
  )
  x <- sim.data[[1]]
  y <- sim.data[[2]]
  true.index <- sim.data[[3]]
  
  # method 1: pro.sgpv
  time.sgpv <- as.numeric(system.time(run.1.sgpv(x, y, family))[3])
  
  # method 2: lasso
  time.lasso <- as.numeric(system.time(run.1.lasso(x, y, family))[3])
  
  # method 3: bess
  time.bess <- as.numeric(system.time(run.1.bess(x, y, family))[3])

  # method 4: isis
  time.isis <- as.numeric(system.time(run.1.isis(x, y, family))[3])
  
  return(c(time.sgpv, time.lasso, time.bess, time.isis))
}


many.sim.supp.fig.7 <- function(num.sim = 300, n, p, s, rho = 0.35, sig = 2,
                                beta.min = 0.1, beta.max = 0.4,
                                intercept = 2,
                                family = c("binomial", "poisson", "cox"),
                                scale = 2, shape = 1, rateC = 0.2) {
  out <- NULL
  out <- replicate(num.sim, one.time.supp.fig.7(
    n = n, p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max,
    intercept = intercept, family = family,
    scale = 2, shape = 1, rateC = 0.2
  )) 
  
  median.sgpv <- median(out[1, ])
  l.sgpv <- as.numeric(quantile(out[1, ], 0.25))
  u.sgpv <- as.numeric(quantile(out[1, ], 0.75))
  
  median.lasso <- median(out[2, ])
  l.lasso <- as.numeric(quantile(out[2, ], 0.25))
  u.lasso <- as.numeric(quantile(out[2, ], 0.75))
  
  median.bess <- median(out[3, ])
  l.bess <- as.numeric(quantile(out[3, ], 0.25))
  u.bess <- as.numeric(quantile(out[3, ], 0.75))
  
  median.isis <- median(out[4, ])
  l.isis <- as.numeric(quantile(out[4, ], 0.25))
  u.isis <- as.numeric(quantile(out[4, ], 0.75))
  
  return(c(
    median.sgpv, l.sgpv, u.sgpv,
    median.lasso, l.lasso, u.lasso,
    median.bess, l.bess, u.bess,
    median.isis, l.isis, u.isis
  ))
}

# function to plot the supp fig 6

get.plot.supp.fig.7 <- function(data, type, p = 200, ybreaks,
                                cap=0.5, np = c("n", "p")) {
  
  data[data>cap] <- cap
  
  # get keywords from data names
  title.p <- ifelse(type == "lds", "low-d sparse",
                    ifelse(type == "ldd", "low-d dense",
                           ifelse(type == "hds", "high-d sparse")
                    )
  )
  
  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")
  
  if (np == "n") {
    plot.d <- data.frame(
      x = rep(seq(2,40,2), each = 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), 20),
      pe = c(data[c(1, 4, 7, 10), ]),
      lb = c(data[c(2, 5, 8, 11), ]),
      ub = c(data[c(3, 6, 9, 12), ])
    )
    
    xlim <- c(1, 40)
    xbreaks <- seq(2, 40, 4)
    xlab <- "n/p"
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, 4 * p, p / 5), each = 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), 16),
      pe = c(data[c(1, 4, 7, 10), ]),
      lb = c(data[c(2, 5, 8, 11), ]),
      ub = c(data[c(3, 6, 9, 12), ])
    )
    
    xlim <- c(p, 4 * p)
    xbreaks <- seq(p, 4 * p, p / 2)
    xlab <- "p"
  }
  
  
  plot.d$Method <- factor(plot.d$Method,
                          levels = c(
                            "ProSGPV", "Lasso", "BeSS", "ISIS"
                          )
  )
  
  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = pe, col = Method)) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = c(0, cap), breaks = ybreaks) +
    labs(
      x = xlab, y = "Running time in seconds", col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = plot.d, aes(
      x = x, ymin = lb,
      ymax = ub, fill = Method
    ), alpha = .3) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


# ------------------------
# RUN ANALYSIS ONE TIME
# ------------------------

spine.one.time <- function(X, Y){
  
  n <- nrow(X)
  
  # generate training and test indices
  train.index <- sample(1:n, n * 0.7, replace = F)
  test.index <- setdiff(1:n, train.index)
  
  # method 1: prosgpv
  sgpv.out <- pro.sgpv(X[train.index, ], Y[train.index], family = "binomial")
  sgpv.index <- sgpv.out$var.index
  sgpv.pred <- predict(sgpv.out, newdata = X[test.index,])
  sgpv.roc <- roc(Y[test.index] ~ sgpv.pred, plot = F, print.auc = F)$auc 
  
  # method 2: lasso
  cv.lasso <- cv.glmnet(as.matrix(X[train.index,]), 
                        Y[train.index], family = "binomial"
  )
  
  lasso.index <- which(coef(cv.lasso, s = cv.lasso$lambda.min)[-1] != 0)
  
  
  lasso.pred <- predict(cv.lasso, s = cv.lasso$lambda.min, 
                        newx = as.matrix(X[test.index, ]),
                        type = "response"
  ) 
  lasso.roc <- roc(Y[test.index] ~ as.vector(lasso.pred), plot = F, print.auc = F)$auc 
  
  # method 3: bess
  bess.fit <- bess(X[train.index, ], Y[train.index], family = "binomial")
  bess.index <- which(bess.fit$beta != 0)
  bess.pred <- predict(bess.fit, newx = X[test.index, ], type = "response")
  bess.roc <- roc(Y[test.index] ~ bess.pred, plot = F, print.auc = F)$auc
 
  # method 4: isis
  
  invisible(capture.output(sis.m <- SIS(as.matrix(X[train.index, ]),
                                        Y[train.index],
                                        family = "binomial", tune = "ebic"
  )))
  sis.index <- sis.m$ix
  
  sis.pred <- predict(sis.m, newx = as.matrix(X[test.index, ]), type = "response")
  
  sis.roc <- roc(Y[test.index] ~ sis.pred, plot = F, print.auc = F)$auc
  
  return(list(
    sgpv.index,
    length(sgpv.index),
    sgpv.roc,
    
    lasso.index,
    length(lasso.index),
    lasso.roc,
    
    bess.index,
    length(bess.index),
    bess.roc,
    
    sis.index,
    length(sis.index),
    sis.roc
    
  ))
  
}

