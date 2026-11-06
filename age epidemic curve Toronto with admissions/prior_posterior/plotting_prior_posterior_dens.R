library(ggplot2)
dir <- './age epidemic curve Toronto with admissions/prior_posterior/'

# Load posterior draws
posterior <- read.csv(paste0(dir,"admissions_posterior_draws.csv"))

# Prior parameters
phi <- log(2) / 0.4

# Density of logit-normal distribution
dlogitnorm <- function(x, meanlogit, sdlogit = 0.5) {
  y <- numeric(length(x))
  valid <- x > 0 & x < 1
  
  y[valid] <- dnorm(
    qlogis(x[valid]),
    mean = meanlogit,
    sd = sdlogit
  ) / (x[valid] * (1 - x[valid]))
  
  y
}

# Function to plot one parameter
plot_prior_posterior <- function(variable, prior_fun, xlim, xlab) {
  
  ggplot() +
    stat_function(
      fun = prior_fun,
      xlim = xlim,
      colour = "black",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_density(
      data = posterior,
      aes(x = .data[[variable]]),
      colour = "black",
      linewidth = 1,
      bounds = xlim
    ) +
    scale_x_continuous(
      limits = xlim,
      breaks = scales::pretty_breaks(n = 6)
    ) +
    labs(
      x = xlab,
      y = "Density"
    ) +
    theme_bw()+
    theme(axis.title.x = element_text(size=14))
}

# Reporting RW standard deviation
gg_pi <- plot_prior_posterior(
  "sd_p",
  function(x) dexp(x, rate = phi),
  c(0, 1.5),
  expression(sigma[pi*"'"])
)

# Hospitalization AR(1) standard deviation
gg_tau <- plot_prior_posterior(
  "sd_tau",
  function(x) dexp(x, rate = phi),
  c(0, 1.5),
  expression(sigma[tau])
)

# Hospitalization AR(1) correlation
gg_rho <- plot_prior_posterior(
  "rho_tau",
  function(x) dbeta(x, 1, 1),
  c(0, 1),
  expression(rho[tau])
)

# Baseline hospitalization probabilities
gg_tau1 <- plot_prior_posterior(
  "tau1",
  function(x) dlogitnorm(x, -4),
  c(0, 0.15),
  expression(logit^{"-1"} ~ group("(", alpha[0-39], ")"))
)

gg_tau2 <- plot_prior_posterior(
  "tau2",
  function(x) dlogitnorm(x, -4),
  c(0, 0.15),
  expression(logit^{"-1"} ~ group("(", alpha[40-59], ")"))
)

gg_tau3 <- plot_prior_posterior(
  "tau3",
  function(x) dlogitnorm(x, -2),
  c(0, 0.35),
  expression(logit^{"-1"} ~ group("(", alpha[60-79], ")"))
)

gg_tau4 <- plot_prior_posterior(
  "tau4",
  function(x) dlogitnorm(x, -2),
  c(0, 0.35),
  expression(logit^{"-1"} ~ group("(", alpha["80+"], ")"))
)

# View plots
gg_pi
gg_tau
gg_rho
gg_tau1
gg_tau2
gg_tau3
gg_tau4

# Save individual figures
plots <- list(
  sigma_pi_prime = gg_pi,
  sigma_tau = gg_tau,
  rho_tau = gg_rho,
  baseline_0_39 = gg_tau1,
  baseline_40_59 = gg_tau2,
  baseline_60_79 = gg_tau3,
  baseline_80_plus = gg_tau4
)

for (nm in names(plots)) {
  ggsave(
    paste0(dir, nm, "_prior_posterior.pdf"),
    plot = plots[[nm]],
    width = 4,
    height = 2
  )
}
