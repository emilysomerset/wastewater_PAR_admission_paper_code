library(ggplot2)
dir <- './age epidemic curve Toronto/prior_posterior/'

# Load posterior draws
posterior <- read.csv(paste0(dir,"case_only_posterior_draws.csv"))

# Prior-posterior plot
gg1 <- ggplot() +
  stat_function(
    fun = dexp, 
    args = list(rate = log(2)/0.4), 
    xlim = c(0, 1.5),
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  )+
  geom_density(
    data = posterior,
    aes(x = sigma_pi),
    colour = "black",
    linewidth = 1,
    bounds = c(0, Inf)
  ) +
  scale_x_continuous(
    name = expression(sigma[pi]),
    limits = c(0, 1.5),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  labs(
    x = expression(sigma[pi]),
    y = "Density"
  ) +
  theme_bw()+ 
  theme(axis.title.x = element_text(size=14))

gg1

# Save figure
ggsave(
  paste0(dir,"case_only_prior_posterior.pdf"),
  plot = gg1,
  width = 4,
  height = 2
)
