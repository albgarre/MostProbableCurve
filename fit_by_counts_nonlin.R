
library(tidyverse)
library(readxl)
library(rstan)
library(coda)
library(cowplot)

##

my_data <- read_excel("./data/Inactivation _55C_original.xlsx", sheet = "FBR14_forR")

my_wes <- wes_palette("FantasticFox1", 5)

my_data %>%
  mutate(i = row_number()) %>%
  group_by(i) %>%
  mutate(average = mean(c(`plate_I`, `plate_II`), na.rm = TRUE)) %>%
  mutate(count = average/vol*1000*10^`dil factor`) %>%
  mutate(logN = log10(count)) %>%
  ggplot(aes(x = time, y = logN, colour = factor(rep))) +
  geom_point() 

## Fit linear regression

d_0 <- my_data %>%
  # filter(rep == 1) %>%
  select(time, plate_I, plate_II, vol, `dil factor`) %>%
  pivot_longer(starts_with("plate")) %>%
  mutate(count = value/vol*1000*10^`dil factor`) %>%
  filter(!is.na(count)) %>%
  mutate(count = as.integer(count),
         logN = log10(count)) %>%
  filter(is.finite(logN))

set.seed(1242)

regression_model <- stan("nonlineal_regression.stan",
                         data = list(t = d_0$time, 
                                     logN = d_0$logN, 
                                     n = nrow(d_0)
                         ),
                         chains = 1, iter = 4000, warmup = 1000
)

regression_model
# plot(regression_model)
post_regression <- As.mcmc.list(regression_model)
# pairs(regression_model)

# as.data.frame(post_regression[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                logN = .$logN0 - (t/.$D)^.$p,
#                N = 10^logN,
#                N0p1 = 10^qnorm(.1, logN, .$sigma),
#                N0p9 = 10^qnorm(.9, logN, .$sigma))
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   ggplot(aes(x = t, colour = i)) +
#   geom_ribbon(aes(ymin = N0p1, ymax = N0p9), fill = NA) +
#   geom_point(aes(x = time, y = count), 
#              inherit.aes = FALSE, data = d) +
#   theme(legend.position = "none") +
#   scale_y_log10()

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

set.seed(24124)

# p_ci_regression <- as.data.frame(post_regression[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                logN = .$logN0 - (t/.$D)^.$p,
#                N = 10^logN,
#                N0p1 = 10^qnorm(.05, logN, .$sigma),
#                N0p9 = 10^qnorm(.95, logN, .$sigma))
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   group_by(t) %>%
#   summarize(q10 = median(N0p1),
#             q90 = median(N0p9)) %>%
#   ggplot() +
#   geom_ribbon(aes(x = t, ymin = q10, ymax = q90), alpha = .5,
#               colour = "black") +
#   geom_point(aes(x = time, y = count), 
#              inherit.aes = FALSE, data = d_0,
#              shape = 1) +
#   theme_cowplot() +
#   scale_y_log10() +
#   xlab("Treatment time (min)") +
#   ylab("Count in the media (CFU/ml)")

p_ci_regression <- as.data.frame(post_regression[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               logN = .$logN0 - (t/.$D)^.$p,
               N = 10^logN,
               N0p1 = 10^qnorm(.05, logN, .$sigma),
               N0p9 = 10^qnorm(.95, logN, .$sigma))
  ) %>%
  imap_dfr(., ~ mutate(.x, i = .y)) %>%
  group_by(t) %>%
  summarize(q10 = median(N0p1),
            q90 = median(N0p9),
            q05_mean = quantile(N, 0.05),
            q95_mean = quantile(N, .95)
            ) %>%
  ggplot() +
  geom_ribbon(aes(x = t, ymin = q10, ymax = q90), alpha = 0,
              linetype = 2,
              colour = my_wes[1],
              fill = my_wes[1]
              ) +
  geom_ribbon(aes(x = t, ymin = q05_mean, ymax = q95_mean), alpha = .5,
              colour = my_wes[1],
              fill = my_wes[1]
              ) +
  geom_point(aes(x = time, y = count), 
             inherit.aes = FALSE, data = d_0,
             shape = 1, size = 3) +
  theme_cowplot() +
  scale_y_log10(breaks = c(1e0, 1e2, 1e4, 1e6, 1e8),
                labels = fancy_scientific(c(1e0, 1e2, 1e4, 1e6, 1e8))) +
  xlab("Treatment time (min)") +
  ylab("Microbial count (CFU/ml)")


## Data for the Poisson fits

d <- my_data %>%
  # filter(rep == 1) %>%
  select(time, plate_I, plate_II, vol, `dil factor`) %>%
  pivot_longer(starts_with("plate")) %>%
  mutate(count = value/vol*1000*10^`dil factor`) %>%
  filter(!is.na(count)) %>%
  mutate(count = as.integer(count),
         logN = log10(count))

# ## Fit the Poisson model
# 
# # set.seed(1241)
# 
# pois_model <- stan("poisson_inactivation.stan",
#                    data = list(t = d$time, 
#                                count = d$count, 
#                                n = nrow(d)
#                    ),
#                    chains = 1, iter = 2000
# )
# pois_model
# plot(pois_model)
# post_pois <- As.mcmc.list(pois_model)
# pairs(pois_model)
# plot(post_pois)
# 
# as.data.frame(post_pois[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                # N = .$N0*10^(-t*.$beta),
#                # logN = log10(N),
#                logN = .$logN0 - t/.$D
#       )
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   ggplot(aes(x = t, y = logN, colour = i)) +
#   geom_line() +
#   geom_point(aes(x = time, y = log10(count)), 
#              inherit.aes = FALSE, data = d) +
#   theme(legend.position = "none")
# 
# as.data.frame(post_pois[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                logN = .$logN0 - t/.$D,
#                N = 10^logN,
#                N0p1 = qpois(.1, N),
#                N0p9 = qpois(.9, N))
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   ggplot(aes(x = t, colour = i)) +
#   geom_ribbon(aes(ymin = N0p1, ymax = N0p9), fill = NA) +
#   geom_point(aes(x = time, y = count), 
#              inherit.aes = FALSE, data = d) +
#   theme(legend.position = "none") +
#   scale_y_log10()

## Fit the non-lineal Poisson model

set.seed(1241)

nonlinpois_model <- stan("nonlin_poisson_inactivation.stan",
                   data = list(t = d$time, 
                               count = d$count, 
                               n = nrow(d)
                   ),
                   chains = 1, iter = 4000, warmup = 1000
)

nonlinpois_model
# plot(nonlinpois_model)
post_nonlinpois <- As.mcmc.list(nonlinpois_model)
# pairs(nonlinpois_model)
# plot(post_nonlinpois)

# as.data.frame(post_nonlinpois[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                # N = .$N0*10^(-t*.$beta),
#                # logN = log10(N),
#                logN = .$logN0 - (t/.$D)^.$p
#       )
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   ggplot(aes(x = t, y = logN, colour = i)) +
#   geom_line() +
#   geom_point(aes(x = time, y = log10(count)), 
#              inherit.aes = FALSE, data = d) +
#   theme(legend.position = "none")

# as.data.frame(post_nonlinpois[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                logN = .$logN0 - (t/.$D)^.$p,
#                N = 10^logN,
#                N0p1 = qpois(.05, N),
#                N0p9 = qpois(.95, N))
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   ggplot(aes(x = t, colour = i)) +
#   geom_ribbon(aes(ymin = N0p1, ymax = N0p9), fill = NA) +
#   geom_point(aes(x = time, y = count), 
#              inherit.aes = FALSE, data = d) +
#   theme(legend.position = "none") +
#   scale_y_log10()

set.seed(214)

d_plot <- d %>%
  mutate(zero = is.infinite(logN)) 

p_ci_pois <- as.data.frame(post_nonlinpois[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               logN = .$logN0 - (t/.$D)^.$p,
               N = 10^logN,
               N0p1 = qpois(.05, N),
               N0p9 = qpois(.95, N))
  ) %>%
  imap_dfr(., ~ mutate(.x, i = .y)) %>%
  group_by(t) %>%
  summarize(q10 = median(N0p1),
            q90 = median(N0p9),
            q05_mean = quantile(N, 0.05),
            q95_mean = quantile(N, .95)
  ) %>%
  ggplot() +
  geom_ribbon(aes(x = t, ymin = q10, ymax = q90), alpha = 0,
              linetype = 2,
              colour = my_wes[2],
              fill = my_wes[2]
              ) +
  geom_ribbon(aes(x = t, ymin = q05_mean, ymax = q95_mean), alpha = .5,
              colour = my_wes[2],
              fill = my_wes[2]
              ) +
  geom_point(aes(x = time, y = count,
                 shape = zero, size = zero), 
             inherit.aes = FALSE, data = d_plot) +
  theme_cowplot() +
  scale_y_log10(breaks = c(1e0, 1e2, 1e4, 1e6, 1e8),
                labels = fancy_scientific(c(1e0, 1e2, 1e4, 1e6, 1e8))) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(3, 3)) +
  xlab("Treatment time (min)") +
  ylab("Microbial count (CFU/ml)") +
  theme(legend.position = "none")

## Fit the non-lineal gamma-Poisson model

set.seed(1241)

nonlin_negbinom_model <- stan("neg_binom_nonlin.stan",
                         data = list(t = d$time, 
                                     count = d$count, 
                                     n = nrow(d)
                         ),
                         chains = 1, iter = 4000, warmup = 1000
)

nonlin_negbinom_model
# plot(nonlin_negbinom_model)
post_nonlin_negbinom <- As.mcmc.list(nonlin_negbinom_model)
# pairs(nonlin_negbinom_model)
# plot(post_nonlin_negbinom)

# as.data.frame(post_nonlin_negbinom[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                # N = .$N0*10^(-t*.$beta),
#                # logN = log10(N),
#                logN = .$logN0 - (t/.$D)^.$p
#       )
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   ggplot(aes(x = t, y = logN, colour = i)) +
#   geom_line() +
#   geom_point(aes(x = time, y = log10(count)), 
#              inherit.aes = FALSE, data = d) +
#   theme(legend.position = "none")

set.seed(121442)

p_ci_negbinom <- as.data.frame(post_nonlin_negbinom[[1]]) %>%
  mutate(i = row_number()) %>%
  sample_n(1000, replace = FALSE) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 20),
               logN = .$logN0 - (t/.$D)^.$p,
               N = 10^logN,
               theta = .$theta,
               N0p1 = qnbinom(.05, mu = N, size = theta),
               N0p9 = qnbinom(.95, mu = N, size = theta)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, i = .y)) %>%
  group_by(t) %>%
  summarize(q10 = median(N0p1),
            q90 = median(N0p9),
            q05_mean = quantile(N, 0.05),
            q95_mean = quantile(N, .95)
  ) %>%
  ggplot() +
  geom_ribbon(aes(x = t, ymin = q10, ymax = q90), alpha = 0,
              linetype = 2,
              colour = my_wes[3],
              fill = my_wes[3]
              ) +
  geom_ribbon(aes(x = t, ymin = q05_mean, ymax = q95_mean), alpha = .5,
              colour = my_wes[3],
              fill = my_wes[3]
              ) +
  geom_point(aes(x = time, y = count,
                 shape = zero, size = zero),
             inherit.aes = FALSE, data = d_plot) +
  theme_cowplot() +
  scale_y_log10(breaks = c(1e0, 1e2, 1e4, 1e6, 1e8),
                labels = fancy_scientific(c(1e0, 1e2, 1e4, 1e6, 1e8))) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(3, 3)) +
  xlab("Treatment time (min)") +
  ylab("Microbial count (CFU/ml)") +
  theme(legend.position = "none")

# ## Fit the model by plates
# 
# set.seed(14212)
# 
# plate_model <- stan("plate_inactivation.stan",
#                     data = list(t = d$time, 
#                                 count = d$value, 
#                                 dil = as.integer(d$`dil factor`),
#                                 vol = d$vol/1000,  # to ml
#                                 n = nrow(d)
#                     ),
#                     chains = 1, iter = 2000
# )
# 
# plate_model
# # plot(plate_model)
# post_plate <- As.mcmc.list(plate_model)
# # pairs(plate_model)
# # plot(post_plate)
# 
# ## Confidence interval in the media
# 
# as.data.frame(post_plate[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                Nmedia = .$N0*10^(-t/.$D),
#                vol = .1,
#                dil = 5,
#                lambda = Nmedia*vol*.1^dil
#       )
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   group_by(t) %>%
#   summarize(q05 = quantile(Nmedia, .05),
#             q95 = quantile(Nmedia, .95),
#             q50 = median(Nmedia)) %>%
#   ggplot() +
#   geom_ribbon(aes(x = t, ymin = q05, ymax = q95)) +
#   geom_line(aes(x = t, y = q50)) +
#   geom_point(aes(x = time, y = count),
#              data = d,
#              inherit.aes = FALSE) +
#   theme(legend.position = "none") +
#   scale_y_log10()
# 
# ## Variation in the plates
# 
# my_sims <- as.data.frame(post_plate[[1]]) %>%
#   sample_n(100) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ mutate(d,
#                Nmedia = .$N0*10^(-time/.$D),
#                lambda = Nmedia*vol/1000*.1^`dil factor`
#       )
#   ) %>%
#   imap_dfr(., ~ mutate(.x, iter = .y)) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map_dfr(.,
#           ~ tibble(plate = rpois(10, .$lambda),
#                    Nobs = plate/.$vol*1000*10^.$`dil factor`,
#                    iter = .$iter,
#                    time = .$time
#           )
#   ) 
# 
# ggplot(my_sims) +
#   geom_point(aes(x = time, y = Nobs, colour = iter)) +
#   geom_point(aes(x = time, y = count), data = d,
#              inherit.aes = FALSE,
#              shape = 1, size = 5) +
#   theme(legend.position = "none") +
#   scale_y_log10()
# 
# ggplot(my_sims) +
#   geom_boxplot(aes(x = factor(time), y = Nobs)) +
#   theme(legend.position = "none") +
#   scale_y_log10()

## Fit the non-lineal plate model

set.seed(14212)

nonlin_plate_model <- stan("nonline_plate_inactivation.stan",
                    data = list(t = d$time, 
                                count = d$value, 
                                dil = as.integer(d$`dil factor`),
                                vol = d$vol/1000,  # to ml
                                n = nrow(d)
                    ),
                    chains = 1, iter = 4000, warmup = 1000
)

nonlin_plate_model
# plot(nonlin_plate_model)
post_nonlinplate <- As.mcmc.list(nonlin_plate_model)
# pairs(nonlin_plate_model)
# plot(post_nonlinplate)

## Confidence interval in the media

# as.data.frame(post_nonlinplate[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                Nmedia = .$N0*10^(-(t/.$D)^.$p),
#                vol = .1,
#                dil = 5,
#                lambda = Nmedia*vol*.1^dil
#       )
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   group_by(t) %>%
#   summarize(q05 = quantile(Nmedia, .05),
#             q95 = quantile(Nmedia, .95),
#             q50 = median(Nmedia)) %>%
#   ggplot() +
#   geom_ribbon(aes(x = t, ymin = q05, ymax = q95)) + 
#   geom_line(aes(x = t, y = q50)) +
#   geom_point(aes(x = time, y = count),
#              data = d,
#              inherit.aes = FALSE) +
#   theme(legend.position = "none") +
#   scale_y_log10()


set.seed(12412)

p_ci_plate <- as.data.frame(post_nonlinplate[[1]]) %>%
  mutate(i = row_number()) %>%
  # sample_n(1000) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               Nmedia = .$N0*10^(-(t/.$D)^.$p),
               vol = .1,
               dil = 5,
               lambda = Nmedia*vol*.1^dil,
               N0p1 = qpois(.05, lambda = Nmedia),
               N0p9 = qpois(.95, lambda = Nmedia)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, i = .y)) %>%
  group_by(t) %>%
  summarize(q10 = median(N0p1),
            q90 = median(N0p9),
            q05_mean = quantile(Nmedia, 0.05),
            q95_mean = quantile(Nmedia, .95)
  ) %>%
  ggplot() +
  geom_ribbon(aes(x = t, ymin = q10, ymax = q90), alpha = .5,
              colour = my_wes[5],
              fill = my_wes[5]
              ) +
  geom_ribbon(aes(x = t, ymin = q05_mean, ymax = q95_mean), alpha = .5,
              colour = my_wes[5],
              fill = my_wes[5]
              ) +
  geom_point(aes(x = time, y = count,
                 shape = zero, size = zero),
             inherit.aes = FALSE, data = d_plot) +
  theme_cowplot() +
  scale_y_log10(breaks = c(1e0, 1e2, 1e4, 1e6, 1e8),
                labels = fancy_scientific(c(1e0, 1e2, 1e4, 1e6, 1e8))) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(3, 3)) +
  xlab("Treatment time (min)") +
  ylab("Microbial count (CFU/ml)") +
  theme(legend.position = "none")



my_line <- as.data.frame(post_nonlinplate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               Nmedia = .$N0*10^(-(t/.$D)^.$p),
               vol = .1,
               ideal_dil = round(log10(Nmedia), 0)-2,
               dil = ifelse(ideal_dil >= 0, ideal_dil, 0),
               lambda = Nmedia*vol*.1^dil,
               q05 = qpois(.05, lambda),
               q95 = qpois(.95, lambda)
      )
  ) %>% 
  imap_dfr(.,
           ~ mutate(.x,
                    iter = .y,
                    N05 = q05/vol*10^dil,
                    N95 = q95/vol*10^dil
           )
  ) %>%
  group_by(t) %>%
  summarize(q05 = median(N05),
            q95 = median(N95)) %>%
  geom_ribbon(aes(x = t, ymin = q05, ymax = q95), alpha = 0,
              colour = my_wes[5],
              linetype = 2,
              data = .,
              inherit.aes = FALSE) 


p_ci_plate <- p_ci_plate + my_line

# p_ci_plate <- as.data.frame(post_nonlinplate[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 250, length = 100),
#                Nmedia = .$N0*10^(-(t/.$D)^.$p),
#                vol = .1,
#                dil = 5,
#                lambda = Nmedia*vol*.1^dil
#       )
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   group_by(t) %>%
#   summarize(q05 = quantile(Nmedia, .05),
#             q95 = quantile(Nmedia, .95),
#             q50 = median(Nmedia)) %>%
#   ggplot() +
#   geom_ribbon(aes(x = t, ymin = q05, ymax = q95), alpha = .5, colour = "black") + 
#   geom_point(aes(x = time, y = count),
#              data = d, shape = 1,
#              inherit.aes = FALSE) +
#   theme_cowplot() +
#   scale_y_log10() +
#   xlab("Treatment time (min)") +
#   ylab("Count in the media (CFU/ml)")

## Figure 3

plot_grid(p_ci_regression + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p_ci_pois + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p_ci_negbinom + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p_ci_plate + coord_cartesian(ylim = c(1e-1, 1e8)),
          labels = "AUTO")

## Variation in the plates

my_sims <- as.data.frame(post_nonlinplate[[1]]) %>%
  sample_n(100) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ mutate(d,
               Nmedia = .$N0*10^(-(time/.$D)^.$p),
               lambda = Nmedia*vol/1000*.1^`dil factor`
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, iter = .y)) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map_dfr(.,
          ~ tibble(plate = rpois(10, .$lambda),
                   Nobs = plate/.$vol*1000*10^.$`dil factor`,
                   iter = .$iter,
                   time = .$time
          )
  ) 

ggplot(my_sims) +
  geom_point(aes(x = time, y = Nobs, colour = iter)) +
  geom_point(aes(x = time, y = count), data = d,
             inherit.aes = FALSE,
             shape = 1, size = 5) +
  theme(legend.position = "none") +
  scale_y_log10()

ggplot(my_sims) +
  geom_boxplot(aes(x = factor(time), y = Nobs)) +
  theme(legend.position = "none") +
  scale_y_log10()

## Confidence interval in the plate

as.data.frame(post_nonlinplate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               Nmedia = .$N0*10^(-(t/.$D)^.$p),
               vol = .1,
               ideal_dil = round(log10(Nmedia), 0)-2,
               dil = ifelse(ideal_dil >= 0, ideal_dil, 0),
               lambda = Nmedia*vol*.1^dil,
               q05 = qpois(.05, lambda),
               q95 = qpois(.95, lambda)
      )
  ) %>% 
  imap_dfr(.,
           ~ mutate(.x,
                    iter = .y,
                    N05 = q05/vol*10^dil,
                    N95 = q95/vol*10^dil
           )
  ) %>%
  group_by(t, dil) %>%
  summarize(q05 = median(N05),
            q95 = median(N95)) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = factor(dil)), alpha = .5) +
  geom_point(aes(x = time, y = count), data = d,
             inherit.aes = FALSE) +
  scale_y_log10()



## Confidence interval in the plate

set.seed(21412)

p1 <- as.data.frame(post_nonlinplate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               Nmedia = .$N0*10^(-(t/.$D)^.$p),
               vol = 1,
               ideal_dil = round(log10(Nmedia), 0)-2,
               dil = ifelse(ideal_dil >= 0, ideal_dil, 0),
               lambda = Nmedia*vol*.1^dil,
               q05 = qpois(.05, lambda),
               q95 = qpois(.95, lambda)
      )
  ) %>% 
  imap_dfr(.,
           ~ mutate(.x,
                    iter = .y,
                    N05 = q05/vol*10^dil,
                    N95 = q95/vol*10^dil
           )
  ) %>%
  group_by(t, dil) %>%
  summarize(q05 = median(N05),
            q95 = median(N95)) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = factor(dil)), alpha = .5) +
  # geom_point(aes(x = time, y = count), data = d,
  #            inherit.aes = FALSE) +
  scale_y_log10() +
  theme_cowplot() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  ylab("Estimated microbial concentration (log CFU/g)") +
  xlab("Treatment time (min)")

set.seed(21412)

p2 <- as.data.frame(post_nonlinplate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               Nmedia = .$N0*10^(-(t/.$D)^.$p),
               vol = .1,
               ideal_dil = round(log10(Nmedia), 0)-2,
               dil = ifelse(ideal_dil >= 0, ideal_dil, 0),
               lambda = Nmedia*vol*.1^dil,
               q05 = qpois(.05, lambda),
               q95 = qpois(.95, lambda)
      )
  ) %>% 
  imap_dfr(.,
           ~ mutate(.x,
                    iter = .y,
                    N05 = q05/vol*10^dil,
                    N95 = q95/vol*10^dil
           )
  ) %>%
  group_by(t, dil) %>%
  summarize(q05 = median(N05),
            q95 = median(N95)) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = factor(dil)), alpha = .5) +
  # geom_point(aes(x = time, y = count), data = d,
  #            inherit.aes = FALSE) +
  scale_y_log10() +
  theme_cowplot() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  ylab("Estimated microbial concentration (log CFU/g)") +
  xlab("Treatment time (min)")

set.seed(21412)

p3 <- as.data.frame(post_nonlinplate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               Nmedia = .$N0*10^(-(t/.$D)^.$p),
               vol = .05,
               ideal_dil = round(log10(Nmedia), 0)-2,
               dil = ifelse(ideal_dil >= 0, ideal_dil, 0),
               lambda = Nmedia*vol*.1^dil,
               q05 = qpois(.05, lambda),
               q95 = qpois(.95, lambda)
      )
  ) %>% 
  imap_dfr(.,
           ~ mutate(.x,
                    iter = .y,
                    N05 = q05/vol*10^dil,
                    N95 = q95/vol*10^dil
           )
  ) %>%
  group_by(t, dil) %>%
  summarize(q05 = median(N05),
            q95 = median(N95)) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = factor(dil)), alpha = .5) +
  # geom_point(aes(x = time, y = count), data = d,
  #            inherit.aes = FALSE) +
  scale_y_log10() +
  theme_cowplot() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  ylab("Estimated microbial concentration (log CFU/g)") +
  xlab("Treatment time (min)")

set.seed(21412)

p4 <- as.data.frame(post_nonlinplate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               Nmedia = .$N0*10^(-(t/.$D)^.$p),
               vol = 1,
               ideal_dil = round(log10(Nmedia), 0)-1,
               dil = ifelse(ideal_dil >= 0, ideal_dil, 0),
               lambda = Nmedia*vol*.1^dil,
               q05 = qpois(.05, lambda),
               q95 = qpois(.95, lambda)
      )
  ) %>% 
  imap_dfr(.,
           ~ mutate(.x,
                    iter = .y,
                    N05 = q05/vol*10^dil,
                    N95 = q95/vol*10^dil
           )
  ) %>%
  group_by(t, dil) %>%
  summarize(q05 = median(N05),
            q95 = median(N95)) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = factor(dil)), alpha = .5) +
  # geom_point(aes(x = time, y = count), data = d,
  #            inherit.aes = FALSE) +
  scale_y_log10() +
  theme_cowplot() +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  ylab("Estimated microbial concentration (log CFU/g)") +
  xlab("Treatment time (min)")

## Supp. Figure 1

plot_grid(p1, p2, p3, p4,
          labels = "AUTO")

## Data for the variability fits

d_var <- my_data %>%
  # filter(rep == 1) %>%
  select(time, plate_I, plate_II, vol, `dil factor`, day, rep) %>%
  pivot_longer(starts_with("plate")) %>%
  mutate(count = value/vol*1000*10^`dil factor`) %>%
  filter(!is.na(value)) %>%
  mutate(count = as.integer(count),
         logN = log10(count),
         day = as.integer(day),
         rep = as.integer(rep))

## Fit the model by plates with variability

# set.seed(14212)
# 
# plate_var_model <- stan("plate_inactivation_variability.stan",
#                         data = list(t = d_var$time, 
#                                     count = d_var$value, 
#                                     dil = as.integer(d_var$`dil factor`),
#                                     vol = d_var$vol/1000,  # to ml
#                                     n = nrow(d_var),
#                                     nbio = length(unique(d_var$day)),
#                                     bio = as.integer(d_var$day)
#                         ),
#                         chains = 1, iter = 8000, warmup = 1000
# )
# 
# plate_var_model
# # plot(plate_var_model)
# post_plate_var <- As.mcmc.list(plate_var_model)
# # pairs(plate_var_model)
# # plot(post_plate_var)
# 
# ## Confidence interval in the media
# 
# as.data.frame(post_plate_var[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 180, length = 100),
#                Nmedia1 = .$`N0[1]`*10^(-t/.$`D[1]`),
#                Nmedia2 = .$`N0[2]`*10^(-t/.$`D[2]`),
#                Nmedia3 = .$`N0[3]`*10^(-t/.$`D[3]`),
#                vol = .1,
#                dil = 5,
#       )
#   ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   group_by(t) %>%
#   summarize(q05_1 = quantile(Nmedia1, .05),
#             q95_1 = quantile(Nmedia1, .95),
#             q50_1 = median(Nmedia1),
#             q05_2 = quantile(Nmedia2, .05),
#             q95_2 = quantile(Nmedia2, .95),
#             q50_2 = median(Nmedia2),
#             q05_3 = quantile(Nmedia3, .05),
#             q95_3 = quantile(Nmedia3, .95),
#             q50_3 = median(Nmedia3)
#   ) %>%
#   ggplot() +
#   geom_ribbon(aes(x = t, ymin = q05_1, ymax = q95_1), alpha = .5, fill = "red") + 
#   geom_line(aes(x = t, y = q50_1), colour = "red") +
#   geom_ribbon(aes(x = t, ymin = q05_2, ymax = q95_2), alpha = .5, fill = "green") +
#   geom_line(aes(x = t, y = q50_2), colour = "green") +
#   geom_ribbon(aes(x = t, ymin = q05_3, ymax = q95_3), alpha = .5, fill = "blue") +
#   geom_line(aes(x = t, y = q50_3), colour = "blue") +
#   geom_point(aes(x = time, y = count, colour = factor(day)),
#              data = d_var, shape = 1,
#              inherit.aes = FALSE) +
#   # theme(legend.position = "none") +
#   scale_y_log10()

## Fit the non-linear model by plates with variability

set.seed(14212)

nonlinplate_var_model <- stan("nonlin_plate_inactivation_variability.stan",
                        data = list(t = d_var$time, 
                                    count = d_var$value, 
                                    dil = as.integer(d_var$`dil factor`),
                                    vol = d_var$vol/1000,  # to ml
                                    n = nrow(d_var),
                                    nbio = length(unique(d_var$day)),
                                    bio = as.integer(d_var$day)
                        ),
                        chains = 1, iter = 8000, warmup = 1000
)

nonlinplate_var_model
# plot(nonlinplate_var_model)
post_nonlinplate_var <- As.mcmc.list(nonlinplate_var_model)
# pairs(nonlinplate_var_model)
# plot(post_nonlinplate_var)

## Confidence interval in the media per repetition - Figure 6

set.seed(24112)

# my_c <- wes_palette("Darjeeling1", 3)
my_c <- wes_palette("FantasticFox1", 3)

as.data.frame(post_nonlinplate_var[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 280, length = 100),
               Nmedia1 = .$`N0[1]`*10^(-(t/.$`D[1]`)^.$`p[1]`),
               Nmedia2 = .$`N0[2]`*10^(-(t/.$`D[2]`)^.$`p[2]`),
               Nmedia3 = .$`N0[3]`*10^(-(t/.$`D[3]`)^.$`p[3]`),
               vol = .1,
               dil = 5,
               pred05_1 = qpois(.05, Nmedia1),
               pred95_1 = qpois(.95, Nmedia1),
               pred05_2 = qpois(.05, Nmedia2),
               pred95_2 = qpois(.95, Nmedia2),
               pred05_3 = qpois(.05, Nmedia3),
               pred95_3 = qpois(.95, Nmedia3)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, i = .y)) %>%
  group_by(t) %>%
  summarize(q05_1 = quantile(Nmedia1, .05),
            q95_1 = quantile(Nmedia1, .95),
            # q50_1 = median(Nmedia1),
            q05_2 = quantile(Nmedia2, .05),
            q95_2 = quantile(Nmedia2, .95),
            # q50_2 = median(Nmedia2),
            q05_3 = quantile(Nmedia3, .05),
            q95_3 = quantile(Nmedia3, .95),
            # q50_3 = median(Nmedia3),
            pred05_1 = median(pred05_1),
            pred95_1 = median(pred95_1),
            pred05_2 = median(pred05_2),
            pred95_2 = median(pred95_2),
            pred05_3 = median(pred05_3),
            pred95_3 = median(pred95_3)
  ) %>%
  ggplot() +
  geom_ribbon(aes(x = t, ymin = q05_1, ymax = q95_1), alpha = .5, fill = my_c[1]) + 
  geom_ribbon(aes(x = t, ymin = pred05_1, ymax = pred95_1), alpha = .5,
              fill = my_c[1], colour = my_c[1]) + 
  # geom_line(aes(x = t, y = q50_1), colour = "red") +
  geom_ribbon(aes(x = t, ymin = q05_2, ymax = q95_2), alpha = .5, fill = my_c[2]) +
  geom_ribbon(aes(x = t, ymin = pred05_2, ymax = pred95_2), alpha = .5, 
              fill = my_c[2], colour = my_c[2]) + 
  # geom_line(aes(x = t, y = q50_2), colour = "green") +
  geom_ribbon(aes(x = t, ymin = q05_3, ymax = q95_3), alpha = .5, fill = my_c[3]) +
  geom_ribbon(aes(x = t, ymin = pred05_3, ymax = pred95_3), alpha = .5, 
              fill = my_c[3], colour = my_c[3]) + 
  # geom_line(aes(x = t, y = q50_3), colour = "blue") +
  geom_point(aes(x = time, y = count, colour = factor(day)),
             data = d_var, shape = 1,
             inherit.aes = FALSE) +
  # theme(legend.position = "none") +
  theme_cowplot() +
  scale_y_log10() +
  xlab("Treatment time (min)") +
  ylab("Count in the media (CFU/ml)") +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  scale_color_manual(values = my_c) +
  coord_cartesian(ylim = c(1e-1, 1e8))

## Model fitted to var 2

d_var %>%
  filter(day == 2) %>%
  filter(is.finite(logN)) %>%
  nls(logN ~ logN0 - (time/D)^p,
      data = .,
      start = list(logN0 = 8,
                   D = 10,
                   p = 1))

## Variability in the D-value

as.data.frame(post_nonlinplate_var[[1]]) %>%
  select(m_logD, s_logD) %>%
  sample_n(10) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(logD = seq(-3, 5, length = 100),
               p = dnorm(logD, .$m_logD, .$s_logD)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  ggplot() +
  geom_line(aes(x = logD, y = p, colour = sim)) +
  theme(legend.position="none")


set.seed(21442)

p_var_logD <- as.data.frame(post_nonlinplate_var[[1]]) %>%
  select(m_logD, s_logD) %>%
  # sample_n(10) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(logD = seq(-3, 5, length = 100),
               p = dnorm(logD, .$m_logD, .$s_logD)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  group_by(logD) %>%
  summarize(q50 = median(p),
            q05 = quantile(p, .05),
            q95 = quantile(p, .95)) %>%
  ggplot() +
  geom_line(aes(x = logD, y = q50), colour = my_wes[5]) +
  geom_ribbon(aes(x = logD, ymin = q05, ymax = q95), alpha = .5, fill = my_wes[5]) +
  theme_bw(base_size = 14) +
  xlab("logarithm of the D-value of strain FBR14 (log min)") +
  ylab("Probability density (·)")


## Variability in the p-value

as.data.frame(post_nonlinplate_var[[1]]) %>%
  select(m_logp, s_logp) %>%
  sample_n(10) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(logp = seq(-2, 2, length = 100),
               p = dnorm(logp, .$m_logp, .$s_logp)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  ggplot() +
  geom_line(aes(x = logp, y = p, colour = sim)) +
  theme(legend.position="none")


set.seed(21442)

p_var_logp <- as.data.frame(post_nonlinplate_var[[1]]) %>%
  select(m_logp, s_logp) %>%
  # sample_n(10) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(logp = seq(-3, 2, length = 100),
               p = dnorm(logp, .$m_logp, .$s_logp)
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  group_by(logp) %>%
  summarize(q50 = median(p),
            q05 = quantile(p, .05),
            q95 = quantile(p, .95)) %>%
  ggplot() +
  geom_line(aes(x = logp, y = q50), colour = my_wes[5]) +
  geom_ribbon(aes(x = logp, ymin = q05, ymax = q95), alpha = .5, fill = my_wes[5]) +
  theme_bw(base_size = 14) +
  xlab("logarithm of the beta-value of strain FBR14 (log min)") +
  ylab("Probability density (·)") 


## Prediction interval with logD variability

dist_logD = as.data.frame(post_nonlinplate_var[[1]]) %>%
  select(m_logD, s_logD, m_logp, s_logp) %>%
  summarize(m_logD = median(m_logD),
            s_logD = median(s_logD),
            m_logp = median(m_logp),
            s_logp = median(s_logp))

p <- tibble(logD = rnorm(n = 1000, mean = dist_logD$m_logD, sd = dist_logD$s_logD),
       logp = rnorm(n = 1000, mean = dist_logD$m_logp, sd = dist_logD$s_logp),
       ) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               logN = 8 - (t/10^.$logD)^10^.$logp)
  ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  group_by(t) %>%
  summarize(m_logN = median(logN),
            q05 = quantile(logN, .05),
            q95 = quantile(logN, .95)) %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = m_logN)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = .5)

my_sims <- as.data.frame(post_nonlinplate_var[[1]]) %>%
  select(m_logD, s_logD, m_logp, s_logp) %>%
  sample_n(100) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(
        logD = rnorm(10, .$m_logD, .$s_logD),
        logp = rnorm(10, .$m_logp, .$s_logp)
      ) %>%
        mutate(j = row_number()) %>%
        split(.$j) %>%
        map(.,
            ~ tibble(t = seq(0, 180, length = 100),
                     logN = 8 - (t/10^.$logD)^10^.$logp)
        ) %>%
        imap_dfr(., ~ mutate(.x, j = .y))
  ) %>%
  imap_dfr(., ~ mutate(.x, i = .y))

# my_sims %>%
#   unite(k, c(i, j)) %>%
#   ggplot() +
#   geom_line(aes(x = t, y = logN, colour = factor(k))) +
#   theme(legend.position = "none")

p1 <- my_sims %>%
  group_by(t) %>%
  summarize(m_logN = median(logN),
            q05 = quantile(logN, .05),
            q95 = quantile(logN, .95)) %>%
  geom_ribbon(aes(x = t, ymin = q05, ymax = q95), alpha = .5, data = .)

p_pred_interval <- p + p1 +
  coord_cartesian(ylim = c(0, 8)) +
  theme_cowplot() +
  xlab("Treatment time (min)") +
  ylab("Log-count in the media (log CFU/ml)")

## Figure 7

plot_grid(p_var_logD, p_var_logp, 
          # p_pred_interval,
          labels = "AUTO",
          ncol = 1)

## Sup. Figure 2

p_pred_interval

## Variation in the plates

my_sims <- as.data.frame(post_nonlinplate_var[[1]]) %>%
  sample_n(100) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ mutate(d,
               Nmedia = .$`N0[2]`*10^(-(time/.$`D[2]`)^.$`p[2]`),
               lambda = Nmedia*vol/1000*.1^`dil factor`
      )
  ) %>%
  imap_dfr(., ~ mutate(.x, iter = .y)) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map_dfr(.,
          ~ tibble(plate = rpois(10, .$lambda),
                   Nobs = plate/.$vol*1000*10^.$`dil factor`,
                   iter = .$iter,
                   time = .$time
          )
  ) 

ggplot(my_sims) +
  geom_point(aes(x = time, y = Nobs, colour = iter)) +
  geom_point(aes(x = time, y = count), 
             data = filter(d_var, day == 2),
             inherit.aes = FALSE,
             shape = 1, size = 5) +
  theme(legend.position = "none") +
  scale_y_log10()

## Confidence interval in the plate

as.data.frame(post_nonlinplate_var[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 250, length = 100),
               Nmedia = .$`N0[3]`*10^(-(t/.$`D[3]`)^.$`p[3]`),
               vol = .1,
               ideal_dil = round(log10(Nmedia), 0)-2,
               dil = ifelse(ideal_dil >= 0, ideal_dil, 0),
               lambda = Nmedia*vol*.1^dil,
               q05 = qpois(.05, lambda),
               q95 = qpois(.95, lambda)
      )
  ) %>% 
  imap_dfr(.,
           ~ mutate(.x,
                    iter = .y,
                    N05 = q05/vol*10^dil,
                    N95 = q95/vol*10^dil
           )
  ) %>%
  group_by(t, dil) %>%
  summarize(q05 = median(N05),
            q95 = median(N95)) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = factor(dil)), alpha = .5) +
  geom_point(aes(x = time, y = count), 
             data = filter(d_var, day == 3),
             inherit.aes = FALSE) +
  scale_y_log10()




