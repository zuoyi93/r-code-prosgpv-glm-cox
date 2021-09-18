# ------------------------------------
# REPLICATE RESULTS IN THE MANUSCRIPT
# ------------------------------------

library("ggplot2")
library("ggpubr")
library("corrr")
library("dplyr")

load("manuscript_data.RData")
source("utils.R")

# -----------------------
# FIGURE 1
# -----------------------

png("Figures/Figure 1.png", units = "in", width = 9, height = 6, res = 300)
fig.1
dev.off()

# -----------------------
# FIGURE 2
# -----------------------

fig.2.logit.ld.s <- get.plot.fig.2(out.main.logit.ld.s, type = "lds")
fig.2.logit.ld.d <- get.plot.fig.2(out.main.logit.ld.d, type = "ldd")
fig.2.logit.hd.s <- get.plot.fig.2(out.main.logit.hd.s, type = "hds", p = 200)

fig.2.logit <- ggarrange(fig.2.logit.ld.s, fig.2.logit.ld.d, fig.2.logit.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

fig.2.logit.3 <- annotate_figure(fig.2.logit,
  top = text_grob("Logistic regression",
    size = 14, face = "bold"
  )
)

fig.2.poisson.ld.s <- get.plot.fig.2(out.main.poisson.ld.s, type = "lds")
fig.2.poisson.ld.d <- get.plot.fig.2(out.main.poisson.ld.d, type = "ldd")
fig.2.poisson.hd.s <- get.plot.fig.2(out.main.poisson.hd.s, type = "hds", p = 80)

fig.2.poisson <- ggarrange(fig.2.poisson.ld.s, fig.2.poisson.ld.d, fig.2.poisson.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

fig.2.poisson.3 <- annotate_figure(fig.2.poisson,
  top = text_grob("Poisson regression",
    size = 14, face = "bold"
  )
)

fig.2.cox.ld.s <- get.plot.fig.2(out.main.cox.ld.s, type = "lds")
fig.2.cox.ld.d <- get.plot.fig.2(out.main.cox.ld.d, type = "ldd")
fig.2.cox.hd.s <- get.plot.fig.2(out.main.cox.hd.s, type = "hds", p = 80)

fig.2.cox <- ggarrange(fig.2.cox.ld.s, fig.2.cox.ld.d, fig.2.cox.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

fig.2.cox.3 <- annotate_figure(fig.2.cox,
  top = text_grob("Cox regression",
    size = 14, face = "bold"
  )
)

png("Figures/Figure 2.png", width = 9, height = 9, units = "in", res = 300)
ggarrange(fig.2.logit.3, fig.2.poisson.3, fig.2.cox.3, ncol = 1)
dev.off()

# -----------------------
# FIGURE 3
# -----------------------

fig.3.logit.ld.s <- get.plot.fig.3(out.main.logit.ld.s, type = "lds", ycap = 0.2, 
                                   ybreaks = seq(0,0.2,0.05))
fig.3.logit.ld.d <- get.plot.fig.3(out.main.logit.ld.d, type = "ldd", ycap = 1, 
                                   ybreaks = seq(0,1,0.2))
fig.3.logit.hd.s <- get.plot.fig.3(out.main.logit.hd.s, type = "hds", p = 200,
                                   ycap = 0.05, ybreaks = seq(0,0.05,0.01))

fig.3.logit <- ggarrange(fig.3.logit.ld.s, fig.3.logit.ld.d, fig.3.logit.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

fig.3.logit.3 <- annotate_figure(fig.3.logit,
  top = text_grob("Logistic regression",
    size = 14, face = "bold"
  )
)

fig.3.poisson.ld.s <- get.plot.fig.3(out.main.poisson.ld.s, type = "lds", ycap = 0.02, 
                                     ybreaks = seq(0,0.02,0.005))
fig.3.poisson.ld.d <- get.plot.fig.3(out.main.poisson.ld.d, type = "ldd", ycap = 0.1, 
                                     ybreaks = seq(0,0.1,0.02))
fig.3.poisson.hd.s <- get.plot.fig.3(out.main.poisson.hd.s, type = "hds", p = 120, 
                                     ycap = 0.01, ybreaks = seq(0,0.01,0.002))

fig.3.poisson <- ggarrange(fig.3.poisson.ld.s, fig.3.poisson.ld.d, fig.3.poisson.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

fig.3.poisson.3 <- annotate_figure(fig.3.poisson,
  top = text_grob("Poisson regression",
    size = 14, face = "bold"
  )
)

fig.3.cox.ld.s <- get.plot.fig.3(out.main.cox.ld.s, type = "lds", ycap = 0.2, 
                                 ybreaks = seq(0,0.2,0.05))
fig.3.cox.ld.d <- get.plot.fig.3(out.main.cox.ld.d, type = "ldd", ycap = 0.2, 
                                 ybreaks = seq(0,0.2,0.05))
fig.3.cox.hd.s <- get.plot.fig.3(out.main.cox.hd.s, type = "hds", p = 80, 
                                 ycap = 0.02, ybreaks = seq(0,0.02,0.005))

fig.3.cox <- ggarrange(fig.3.cox.ld.s, fig.3.cox.ld.d, fig.3.cox.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

fig.3.cox.3 <- annotate_figure(fig.3.cox,
  top = text_grob("Cox regression",
    size = 14, face = "bold"
  )
)

png("Figures/Figure 3.png", width = 9, height = 9, units = "in", res = 300)
ggarrange(fig.3.logit.3, fig.3.poisson.3, fig.3.cox.3, ncol = 1)
dev.off()

# -----------------------
# FIGURE 4
# -----------------------

fig.4.logit.ld.s <- get.plot.fig.4(out.main.logit.ld.s, type = "lds", model = "logistic")
fig.4.logit.ld.d <- get.plot.fig.4(out.main.logit.ld.d, type = "ldd", model = "logistic")
fig.4.logit.hd.s <- get.plot.fig.4(out.main.logit.hd.s, type = "hds", p = 200, model = "logistic")

fig.4.logit <- ggarrange(fig.4.logit.ld.s, fig.4.logit.ld.d, fig.4.logit.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

fig.4.logit.3 <- annotate_figure(fig.4.logit,
  top = text_grob("Logistic regression",
    size = 14, face = "bold"
  )
)

fig.4.poisson.ld.s <- get.plot.fig.4(out.main.poisson.ld.s, type = "lds", model = "poisson",
                                     ycap = 8, ybreaks = seq(0,8,2))
fig.4.poisson.ld.d <- get.plot.fig.4(out.main.poisson.ld.d, type = "ldd", model = "poisson",
                                     ycap = 40, ybreaks = seq(0,40,10))
fig.4.poisson.hd.s <- get.plot.fig.4(out.main.poisson.hd.s, type = "hds", p = 120, 
                                     model = "poisson", ycap = 50, ybreaks = seq(0,50,10))

fig.4.poisson <- ggarrange(fig.4.poisson.ld.s, fig.4.poisson.ld.d, fig.4.poisson.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

fig.4.poisson.3 <- annotate_figure(fig.4.poisson,
  top = text_grob("Poisson regression",
    size = 14, face = "bold"
  )
)

png("Figures/Figure 4.png", width = 9, height = 6, units = "in", res = 300)
ggarrange(fig.4.logit.3, fig.4.poisson.3, ncol = 1)
dev.off()

# ----------------------------------------------
# SUPPLEMENTARY FIGURE 1: SENSITIVITY OF LAMBDA
# ----------------------------------------------

supp.fig.1.logit.ld.s <- get.plot.supp.fig.1(out.supp.fig.1.logit.ld.s, type = "lds")
supp.fig.1.logit.ld.d <- get.plot.supp.fig.1(out.supp.fig.1.logit.ld.d, type = "ldd")
supp.fig.1.logit.hd.s <- get.plot.supp.fig.1(out.supp.fig.1.logit.hd.s, type = "hds", p = 200)

supp.fig.1.logit <- ggarrange(supp.fig.1.logit.ld.s, supp.fig.1.logit.ld.d, supp.fig.1.logit.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.1.logit.3 <- annotate_figure(supp.fig.1.logit,
  top = text_grob("Logistic regression",
    size = 14, face = "bold"
  )
)

supp.fig.1.poisson.ld.s <- get.plot.supp.fig.1(out.supp.fig.1.poisson.ld.s, type = "lds")
supp.fig.1.poisson.ld.d <- get.plot.supp.fig.1(out.supp.fig.1.poisson.ld.d, type = "ldd")
supp.fig.1.poisson.hd.s <- get.plot.supp.fig.1(out.supp.fig.1.poisson.hd.s, type = "hds", p = 80)

supp.fig.1.poisson <- ggarrange(supp.fig.1.poisson.ld.s, supp.fig.1.poisson.ld.d, supp.fig.1.poisson.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.1.poisson.3 <- annotate_figure(supp.fig.1.poisson,
  top = text_grob("Poisson regression",
    size = 14, face = "bold"
  )
)

supp.fig.1.cox.ld.s <- get.plot.supp.fig.1(out.supp.fig.1.cox.ld.s, type = "lds")
supp.fig.1.cox.ld.d <- get.plot.supp.fig.1(out.supp.fig.1.cox.ld.d, type = "ldd")
supp.fig.1.cox.hd.s <- get.plot.supp.fig.1(out.supp.fig.1.cox.hd.s, type = "hds", p = 80)

supp.fig.1.cox <- ggarrange(supp.fig.1.cox.ld.s, supp.fig.1.cox.ld.d, supp.fig.1.cox.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.1.cox.3 <- annotate_figure(supp.fig.1.cox,
  top = text_grob("Cox regression",
    size = 14, face = "bold"
  )
)

png("Figures/Supp Fig 1.png", width = 9, height = 9, units = "in", res = 300)
ggarrange(supp.fig.1.logit.3, supp.fig.1.poisson.3, supp.fig.1.cox.3, ncol = 1)
dev.off()

# -----------------------------------
# SUPPLEMENTARY FIGURE 2: NULL BOUND
# -----------------------------------

gp.bound.supp.fig.2.ind <- get.plot.supp.fig.2.bound(
  data = supp.fig.2.ind, cor = "ind"
)

gp.capture.supp.fig.2.ind <- get.plot.supp.fig.2.capture(
  data = supp.fig.2.ind, cor = "ind"
)

gp.bound.supp.fig.2.medium <- get.plot.supp.fig.2.bound(
  data = supp.fig.2.medium, cor = "medium"
)

gp.capture.supp.fig.2.medium <- get.plot.supp.fig.2.capture(
  data = supp.fig.2.medium, cor = "medium"
)

gp.bound.supp.fig.2.high <- get.plot.supp.fig.2.bound(
  data = supp.fig.2.high, cor = "high"
)

gp.capture.supp.fig.2.high <- get.plot.supp.fig.2.capture(
  data = supp.fig.2.high, cor = "high"
)

bound.supp.fig.2.p <- ggarrange(gp.bound.supp.fig.2.ind,
  gp.bound.supp.fig.2.medium,
  gp.bound.supp.fig.2.high,
  ncol = 3, common.legend = T, legend = "bottom"
)

bound.supp.fig.2 <- annotate_figure(bound.supp.fig.2.p,
  top = text_grob("Null bound size, logistic regression, p = 20, s = 4",
    size = 14, face = "bold"
  )
)

capture.supp.fig.2.p <- ggarrange(gp.capture.supp.fig.2.ind,
  gp.capture.supp.fig.2.medium,
  gp.capture.supp.fig.2.high,
  ncol = 3, common.legend = T, legend = "bottom"
)

capture.supp.fig.2 <- annotate_figure(capture.supp.fig.2.p,
  top = text_grob("Support recovery, logistic regression, p = 20, s = 4",
    size = 14, face = "bold"
  )
)

png("Figures/Supp Fig 2.png", width = 9, height = 6, units = "in", res = 300)
ggarrange(bound.supp.fig.2,
  capture.supp.fig.2,
  ncol = 1
)
dev.off()

# ----------------------------------------------
# SUPPLEMENTARY FIGURE 3: ML AND JEFFREYS PRIOR
# ----------------------------------------------

supp.fig.3.rate.low.s <- get.plot.supp.fig.3(
  data = out.supp.fig.3.low.s, 
  outcome = "rate", x.n = seq(2,40,2), xbreaks = seq(2,40,4)) 

supp.fig.3.pe.low.s <- get.plot.supp.fig.3(
  data = out.supp.fig.3.low.s, pe.flag = "ld.s",
  outcome = "pe", x.n = seq(2,40,2), xbreaks = seq(2,40,4)) 

supp.fig.3.pr.low.s <- get.plot.supp.fig.3(
  data = out.supp.fig.3.low.s, 
  outcome = "pr", x.n = seq(2,40,2), xbreaks = seq(2,40,4)) 

supp.fig.3.rate.low.d <- get.plot.supp.fig.3(
  data = out.supp.fig.3.low.d, 
  outcome = "rate", x.n = seq(2,40,2), xbreaks = seq(2,40,4)) 

supp.fig.3.pe.low.d <- get.plot.supp.fig.3(
  data = out.supp.fig.3.low.d, pe.flag = "ld.d",
  outcome = "pe", x.n = seq(2,40,2), xbreaks = seq(2,40,4)) 

supp.fig.3.pr.low.d <- get.plot.supp.fig.3(
  data = out.supp.fig.3.low.d, 
  outcome = "pr", x.n = seq(2,40,2), xbreaks = seq(2,40,4)) 

supp.fig.3.rate.high.s <- get.plot.supp.fig.3(
  data = out.supp.fig.3.high.s, 
  outcome = "rate", x.n = seq(200,800,40), xbreaks = seq(200,800,100)) 

supp.fig.3.pe.high.s <- get.plot.supp.fig.3(
  data = out.supp.fig.3.high.s, pe.flag = "hd.s",
  outcome = "pe", x.n = seq(200,800,40), xbreaks = seq(200,800,100)) 

supp.fig.3.pr.high.s <- get.plot.supp.fig.3(
  data = out.supp.fig.3.high.s, 
  outcome = "pr", x.n = seq(200,800,40), xbreaks = seq(200,800,100)) 


low.s.supp.fig.3.p <- ggarrange(supp.fig.3.rate.low.s,
                                supp.fig.3.pe.low.s,
                                supp.fig.3.pr.low.s,
                                ncol = 3, common.legend = T, legend = "bottom"
)

low.s.supp.fig.3 <- annotate_figure(low.s.supp.fig.3.p,
                                    top = text_grob("Low-dimensional sparse signals, p = 20, s = 4",
                                                    size = 14, face = "bold"
                                    )
)

low.d.supp.fig.3.p <- ggarrange(supp.fig.3.rate.low.d,
                                supp.fig.3.pe.low.d,
                                supp.fig.3.pr.low.d,
                                ncol = 3, common.legend = T, legend = "bottom"
)

low.d.supp.fig.3 <- annotate_figure(low.d.supp.fig.3.p,
                                    top = text_grob("Low-dimensional dense signals, p = 20, s = 14",
                                                    size = 14, face = "bold"
                                    )
)

high.s.supp.fig.3.p <- ggarrange(supp.fig.3.rate.high.s,
                                supp.fig.3.pe.high.s,
                                supp.fig.3.pr.high.s,
                                ncol = 3, common.legend = T, legend = "bottom"
)

high.s.supp.fig.3 <- annotate_figure(high.s.supp.fig.3.p,
                                    top = text_grob("High-dimensional sparse signals, n = 200, s = 4",
                                                    size = 14, face = "bold"
                                    )
)

png("Figures/Supp Fig 3.png", width = 9, height = 9, units = "in", res = 300)
ggarrange(low.s.supp.fig.3,
          low.d.supp.fig.3,
          high.s.supp.fig.3,
          ncol = 1
)
dev.off()


# --------------------------------------------
# SUPPLEMENTARY FIGURE 4: TYPE I AND II ERROR
# --------------------------------------------

supp.fig.4.logit.ld.s <- get.plot.supp.fig.4(out.main.logit.ld.s[c(3, 14, 25, 36), ],
  out.main.logit.ld.s[c(2, 13, 24, 35), ],
  type = "lds"
)
supp.fig.4.logit.ld.d <- get.plot.supp.fig.4(out.main.logit.ld.d[c(3, 14, 25, 36), ],
  out.main.logit.ld.d[c(2, 13, 24, 35), ],
  type = "ldd"
)
supp.fig.4.logit.hd.s <- get.plot.supp.fig.4(out.main.logit.hd.s[c(3, 14, 25, 36), ],
  out.main.logit.hd.s[c(2, 13, 24, 35), ],
  type = "hds", p = 200
)

supp.fig.4.logit <- ggarrange(supp.fig.4.logit.ld.s, supp.fig.4.logit.ld.d, supp.fig.4.logit.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.4.logit.3 <- annotate_figure(supp.fig.4.logit,
  top = text_grob("Logistic regression",
    size = 14, face = "bold"
  )
)

supp.fig.4.poisson.ld.s <- get.plot.supp.fig.4(out.main.poisson.ld.s[c(3, 14, 25, 36), ],
  out.main.poisson.ld.s[c(2, 13, 24, 35), ],
  type = "lds"
)
supp.fig.4.poisson.ld.d <- get.plot.supp.fig.4(out.main.poisson.ld.d[c(3, 14, 25, 36), ],
  out.main.poisson.ld.d[c(2, 13, 24, 35), ],
  type = "ldd"
)
supp.fig.4.poisson.hd.s <- get.plot.supp.fig.4(out.main.poisson.hd.s[c(3, 14, 25, 36), ],
  out.main.poisson.hd.s[c(2, 13, 24, 35), ],
  type = "hds", p = 80
)

supp.fig.4.poisson <- ggarrange(supp.fig.4.poisson.ld.s, supp.fig.4.poisson.ld.d, supp.fig.4.poisson.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.4.poisson.3 <- annotate_figure(supp.fig.4.poisson,
  top = text_grob("Poisson regression",
    size = 14, face = "bold"
  )
)

supp.fig.4.cox.ld.s <- get.plot.supp.fig.4(out.main.cox.ld.s[c(3, 14, 25, 36), ],
  out.main.cox.ld.s[c(2, 13, 24, 35), ],
  type = "lds"
)
supp.fig.4.cox.ld.d <- get.plot.supp.fig.4(out.main.cox.ld.d[c(3, 14, 25, 36), ],
  out.main.cox.ld.d[c(2, 13, 24, 35), ],
  type = "ldd"
)
supp.fig.4.cox.hd.s <- get.plot.supp.fig.4(out.main.cox.hd.s[c(3, 14, 25, 36), ],
  out.main.cox.hd.s[c(2, 13, 24, 35), ],
  type = "hds", p = 80
)

supp.fig.4.cox <- ggarrange(supp.fig.4.cox.ld.s, supp.fig.4.cox.ld.d, supp.fig.4.cox.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.4.cox.3 <- annotate_figure(supp.fig.4.cox,
  top = text_grob("Cox regression",
    size = 14, face = "bold"
  )
)

png("Figures/Supp Fig 4.png", width = 9, height = 9, units = "in", res = 300)
ggarrange(supp.fig.4.logit.3, supp.fig.4.poisson.3, supp.fig.4.cox.3, ncol = 1)
dev.off()


# -------------------------------------
# SUPPLEMENTARY FIGURE 5: FDR AND FNDR
# -------------------------------------

supp.fig.5.logit.ld.s <- get.plot.supp.fig.5(out.main.logit.ld.s[c(4, 15, 26, 37), ],
  out.main.logit.ld.s[c(5, 16, 27, 38), ],
  type = "lds"
)
supp.fig.5.logit.ld.d <- get.plot.supp.fig.5(out.main.logit.ld.d[c(4, 15, 26, 37), ],
  out.main.logit.ld.d[c(5, 16, 27, 38), ],
  type = "ldd"
)
supp.fig.5.logit.hd.s <- get.plot.supp.fig.5(out.main.logit.hd.s[c(4, 15, 26, 37), ],
  out.main.logit.hd.s[c(5, 16, 27, 38), ],
  type = "hds", p = 200
)

supp.fig.5.logit <- ggarrange(supp.fig.5.logit.ld.s, supp.fig.5.logit.ld.d, supp.fig.5.logit.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.5.logit.3 <- annotate_figure(supp.fig.5.logit,
  top = text_grob("Logistic regression",
    size = 14, face = "bold"
  )
)

supp.fig.5.poisson.ld.s <- get.plot.supp.fig.5(out.main.poisson.ld.s[c(4, 15, 26, 37), ],
  out.main.poisson.ld.s[c(5, 16, 27, 38), ],
  type = "lds"
)
supp.fig.5.poisson.ld.d <- get.plot.supp.fig.5(out.main.poisson.ld.d[c(4, 15, 26, 37), ],
  out.main.poisson.ld.d[c(5, 16, 27, 38), ],
  type = "ldd"
)
supp.fig.5.poisson.hd.s <- get.plot.supp.fig.5(out.main.poisson.hd.s[c(4, 15, 26, 37), ],
  out.main.poisson.hd.s[c(5, 16, 27, 38), ],
  type = "hds", p = 80
)

supp.fig.5.poisson <- ggarrange(supp.fig.5.poisson.ld.s, supp.fig.5.poisson.ld.d, supp.fig.5.poisson.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.5.poisson.3 <- annotate_figure(supp.fig.5.poisson,
  top = text_grob("Poisson regression",
    size = 14, face = "bold"
  )
)

supp.fig.5.cox.ld.s <- get.plot.supp.fig.5(out.main.cox.ld.s[c(4, 15, 26, 37), ],
  out.main.cox.ld.s[c(5, 16, 27, 38), ],
  type = "lds"
)
supp.fig.5.cox.ld.d <- get.plot.supp.fig.5(out.main.cox.ld.d[c(4, 15, 26, 37), ],
  out.main.cox.ld.d[c(5, 16, 27, 38), ],
  type = "ldd"
)
supp.fig.5.cox.hd.s <- get.plot.supp.fig.5(out.main.cox.hd.s[c(4, 15, 26, 37), ],
  out.main.cox.hd.s[c(5, 16, 27, 38), ],
  type = "hds", p = 80
)

supp.fig.5.cox <- ggarrange(supp.fig.5.cox.ld.s, supp.fig.5.cox.ld.d, supp.fig.5.cox.hd.s,
  ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.5.cox.3 <- annotate_figure(supp.fig.5.cox,
  top = text_grob("Cox regression",
    size = 14, face = "bold"
  )
)

png("Figures/Supp Fig 5.png", width = 9, height = 9, units = "in", res = 300)
ggarrange(supp.fig.5.logit.3, supp.fig.5.poisson.3, supp.fig.5.cox.3, ncol = 1)
dev.off()

# ---------------------------------------------
# SUPPLEMENTARY FIGURE 6: LOW-D DENSE LOGISTIC
# ---------------------------------------------

cap.medium.cor <- get.plot.supp.fig.6.cap(data = out.supp.fig.6.1)
pe.medium.cor <- get.plot.supp.fig.6.pe(data = out.supp.fig.6.1, ycap = 1, ybreaks = seq(0,1,0.2))
pr.medium.cor <- get.plot.supp.fig.6.pr(data = out.supp.fig.6.1)

cap.high.cor <- get.plot.supp.fig.6.cap(data = out.supp.fig.6.2)
pe.high.cor <- get.plot.supp.fig.6.pe(data = out.supp.fig.6.2, ycap = 1, ybreaks = seq(0,1,0.2))
pr.high.cor <- get.plot.supp.fig.6.pr(data = out.supp.fig.6.2)

supp.fig.6.medium <- ggarrange(cap.medium.cor,
                               pe.medium.cor,
                               pr.medium.cor,
                               ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.6.medium.p <- annotate_figure(supp.fig.6.medium,
                                       top = text_grob("Logistic regression, p = 20, s = 14, rho = 0.35",
                                                       size = 14, face = "bold"
                                       )
)

supp.fig.6.high <- ggarrange(cap.high.cor,
                             pe.high.cor,
                             pr.high.cor,
                             ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.6.high.p <- annotate_figure(supp.fig.6.high,
                                     top = text_grob("Logistic regression, p = 20, s = 14, rho = 0.7",
                                                     size = 14, face = "bold"
                                     )
)

png("Figures/Supp Fig 6.png", width = 9, height = 6, units = "in", res = 300)
ggarrange(supp.fig.6.medium.p,
          supp.fig.6.high.p,
          ncol = 1
)
dev.off()


# -----------------------------------------
# SUPPLEMENTARY FIGURE 7: COMPUTATION TIME
# -----------------------------------------

supp.fig.7.logit.ld.s <- get.plot.supp.fig.7(
  data = out.supp.fig.7.logit.ld.s, type = "lds", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "n")

supp.fig.7.logit.ld.d <- get.plot.supp.fig.7(
  data = out.supp.fig.7.logit.ld.d, type = "ldd", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "n")

supp.fig.7.logit.hd.s <- get.plot.supp.fig.7(
  data = out.supp.fig.7.logit.hd.s, type = "hds", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "p")

supp.fig.7.poisson.ld.s <- get.plot.supp.fig.7(
  data = out.supp.fig.7.poisson.ld.s, type = "lds", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "n")

supp.fig.7.poisson.ld.d <- get.plot.supp.fig.7(
  data = out.supp.fig.7.poisson.ld.d, type = "ldd", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "n")

supp.fig.7.poisson.hd.s <- get.plot.supp.fig.7(
  data = out.supp.fig.7.poisson.hd.s, type = "hds", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "p")

supp.fig.7.cox.ld.s <- get.plot.supp.fig.7(
  data = out.supp.fig.7.cox.ld.s, type = "lds", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "n")

supp.fig.7.cox.ld.d <- get.plot.supp.fig.7(
  data = out.supp.fig.7.cox.ld.d, type = "ldd", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "n")

supp.fig.7.cox.hd.s <- get.plot.supp.fig.7(
  data = out.supp.fig.7.cox.hd.s, type = "hds", cap=1.2, ybreaks = seq(0,1.2,0.3), np = "p")

supp.fig.7.logit <- ggarrange(supp.fig.7.logit.ld.s, 
                              supp.fig.7.logit.ld.d, 
                              supp.fig.7.logit.hd.s,
                              ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.7.logit.3 <- annotate_figure(supp.fig.7.logit,
                                        top = text_grob("Logistic regression",
                                                        size = 14, face = "bold"
                                        )
)

supp.fig.7.poisson <- ggarrange(supp.fig.7.poisson.ld.s, 
                              supp.fig.7.poisson.ld.d, 
                              supp.fig.7.poisson.hd.s,
                              ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.7.poisson.3 <- annotate_figure(supp.fig.7.poisson,
                                      top = text_grob("Poisson regression",
                                                      size = 14, face = "bold"
                                      )
)

supp.fig.7.cox <- ggarrange(supp.fig.7.cox.ld.s, 
                                supp.fig.7.cox.ld.d, 
                                supp.fig.7.cox.hd.s,
                                ncol = 3, common.legend = T, legend = "bottom"
)

supp.fig.7.cox.3 <- annotate_figure(supp.fig.7.cox,
                                        top = text_grob("Cox regression",
                                                        size = 14, face = "bold"
                                        )
)

png("Figures/Supp Fig 6.png", width = 9, height = 9, units = "in", res = 300)
ggarrange(supp.fig.7.logit.3, supp.fig.7.poisson.3, supp.fig.7.cox.3, ncol = 1)
dev.off()


# ------------------------
# REAL WORLD ANALYSIS
# ------------------------

# ------------------------
# CORRELATION AND CLUSTERS
# ------------------------

png("Figures/Supp Fig 8.png",
  width = 8, height = 8, units = "in", res = 300
)
spine %>%
  correlate() %>%
  network_plot(min_cor = 0.1, colours = c("red", "white", "blue"))
dev.off()

# -------------------------
# MODEL SIZE 
# -------------------------

png("Figures/Supp Fig 9.png",
    width = 8, height = 6, units = "in", res = 300
)
spine.hist
dev.off()


# -------------------------
# PREDICTION PERFORMANCE
# -------------------------

png("Figures/Supp Fig 10.png",
    width = 6, height = 6, units = "in", res = 300
)
spine.box
dev.off()


