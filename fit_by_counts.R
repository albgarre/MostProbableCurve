
library(tidyverse)
library(readxl)
library(rstan)
library(coda)
library(cowplot)
library(wesanderson)

##

my_data <- read_excel("./data/Inactivation _55C_original.xlsx", sheet = "FBR16_forR")

my_wes <- wes_palette("FantasticFox1", 5)

# my_data %>%
#   mutate(i = row_number()) %>%
#   group_by(i) %>%
#   mutate(average = mean(c(`plate_I`, `plate_II`), na.rm = TRUE)) %>%
#   mutate(count = average/vol*1000*10^`dil factor`) %>% 
#   mutate(count - `cfu/ml`) %>% View()
  

my_data %>%
  mutate(i = row_number()) %>%
  group_by(i) %>%
  mutate(average = mean(c(`plate_I`, `plate_II`), na.rm = TRUE)) %>%
  mutate(count = average/vol*1000*10^`dil factor`) %>%
  mutate(logN = log10(count)) %>%
  ggplot(aes(x = time, y = logN, colour = factor(rep))) +
  geom_point() 
  # geom_line()

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

regression_model <- stan("lineal_regression.stan",
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
#       ~ tibble(t = seq(0, 180, length = 100),
#                logN = .$logN0 - t/.$D,
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

p_ci_regression <- as.data.frame(post_regression[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               logN = .$logN0 - t/.$D,
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
  
## Fit the Poisson model

set.seed(8212)

pois_model <- stan("poisson_inactivation.stan",
                   data = list(t = d$time, 
                               count = d$count, 
                               n = nrow(d)
                               ),
                   chains = 1, iter = 4000, warmup = 1000
                   )
pois_model
# plot(pois_model)
post_pois <- As.mcmc.list(pois_model)
# pairs(pois_model)
# plot(post_pois)

# as.data.frame(post_pois[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 180, length = 100),
#                # N = .$N0*10^(-t*.$beta),
#                # logN = log10(N),
#                logN = .$logN0 - t/.$D
#                )
#       ) %>%
#   imap_dfr(., ~ mutate(.x, i = .y)) %>%
#   ggplot(aes(x = t, y = logN, colour = i)) +
#   geom_line() +
#   geom_point(aes(x = time, y = log10(count)), 
#              inherit.aes = FALSE, data = d) +
#   theme(legend.position = "none")

# as.data.frame(post_pois[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 180, length = 100),
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

set.seed(21412)

d_plot <- d %>%
  mutate(zero = is.infinite(logN)) 

p_ci_pois <- as.data.frame(post_pois[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               logN = .$logN0 - t/.$D,
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
              colour = my_wes[2],
              fill = my_wes[2],
              linetype = 2
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

## Fit the gamma-Poisson model

set.seed(42124)

negbin_model <- stan("neg_binomial.stan",
                   data = list(t = d$time, 
                               count = d$count, 
                               n = nrow(d)
                   ),
                   chains = 1, iter = 4000, warmup = 1000
)

# tibble(x = seq(0, 10, length = 100),
#        y = dgamma(x, 1, 1)) %>%
#   ggplot() +geom_line(aes(x, y))

negbin_model
# plot(negbin_model)
post_negbin <- As.mcmc.list(negbin_model)
# pairs(negbin_model)

# as.data.frame(post_negbin[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 180, length = 100),
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

set.seed(12412)

p_ci_negbinom <- as.data.frame(post_negbin[[1]]) %>%
  mutate(i = row_number()) %>%
  sample_n(1000, replace = FALSE) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 20),
               logN = .$logN0 - t/.$D,
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

## Fit the model by plates

set.seed(14212)

plate_model <- stan("plate_inactivation.stan",
                   data = list(t = d$time, 
                               count = d$value, 
                               dil = as.integer(d$`dil factor`),
                               vol = d$vol/1000,  # to ml
                               n = nrow(d)
                   ),
                   chains = 1, iter = 4000, warmup = 1000
)

plate_model
# plot(plate_model)
post_plate <- As.mcmc.list(plate_model)
# pairs(plate_model)
# plot(post_plate)

## Confidence interval in the media

# as.data.frame(post_plate[[1]]) %>%
#   mutate(i = row_number()) %>%
#   split(.$i) %>%
#   map(.,
#       ~ tibble(t = seq(0, 180, length = 100),
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

set.seed(12412)

p_ci_plate <- as.data.frame(post_plate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               Nmedia = .$N0*10^(-t/.$D),
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

my_line <- as.data.frame(post_plate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               Nmedia = .$N0*10^(-t/.$D),
               vol = 0.1,
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
              linetype = 2,
              colour = my_wes[5],
              data = .,
              inherit.aes = FALSE)

p_ci_plate <- p_ci_plate + my_line

## Figure 1

plot_grid(p_ci_regression + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p_ci_pois + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p_ci_negbinom + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p_ci_plate + coord_cartesian(ylim = c(1e-1, 1e8)),
          labels = "AUTO")

## Variation in the plates

my_sims <- as.data.frame(post_plate[[1]]) %>%
  sample_n(100) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ mutate(d,
             Nmedia = .$N0*10^(-time/.$D),
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
               time = .$time,
               dil = .$`dil factor`,
               vol = .$vol
               )
      ) 

## Supp. Figure ??

p1 <- ggplot(my_sims) +
  geom_point(aes(x = time, y = Nobs, colour = iter)) +
  geom_point(aes(x = time, y = count), data = d,
             inherit.aes = FALSE,
             shape = 1, size = 5) +
  theme(legend.position = "none") +
  scale_y_log10() +
  ylab("Calculated microbial concentration (CFU/ml)") +
  xlab("Treatment time")

p2 <- my_sims %>%
  # filter(time %in% c(0, 2 100)) %>%
  ggplot() +
  # geom_histogram(aes(Nobs, colour = factor(dil))) +
  geom_density(aes(Nobs, colour = factor(dil)), kernel = "rectangular") +
  geom_vline(aes(xintercept = count, colour = factor(`dil factor`)), 
             data = d, 
             linetype = 2) +
  facet_wrap("time", scales = "free") +
  xlab("Calculated microbial concentration (CFU/ml)") + 
  ylab("")

# my_sims %>%
#   # filter(time %in% c(0, 2 100)) %>%
#   ggplot() +
#   geom_histogram(aes(Nobs)) +
#   # geom_density(aes(Nobs, colour = dil))
#   # geom_density(aes(Nobs), kernel = "triangular") +
#   geom_vline(aes(xintercept = count), data = d, # filter(d, time %in% c(0, 100)),
#              inherit.aes = FALSE,
#              colour = "red", linetype = 2) +
#   facet_wrap("time", scales = "free") +
#   xlab("Calculated microbial concentration (CFU/ml)") + 
#   ylab("")

plot_grid(p1, p2,
          ncol = 1)

ggplot(my_sims) +
  geom_boxplot(aes(x = factor(time), y = Nobs)) +
  theme(legend.position = "none") +
  scale_y_log10()

## Confidence interval in the plate

set.seed(21412)

p1 <- as.data.frame(post_plate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               Nmedia = .$N0*10^(-t/.$D),
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
  ylab("Estimated microbial concentration (CFU/ml)") +
  xlab("Treatment time (min)")

set.seed(21412)

p2 <- as.data.frame(post_plate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               Nmedia = .$N0*10^(-t/.$D),
               vol = 0.1,
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
  ylab("Estimated microbial concentration (CFU/ml)") +
  xlab("Treatment time (min)")

set.seed(21412)

p3 <- as.data.frame(post_plate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               Nmedia = .$N0*10^(-t/.$D),
               vol = 0.05,
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
  ylab("Estimated microbial concentration (CFU/ml)") +
  xlab("Treatment time (min)")

set.seed(21412)

p4 <- as.data.frame(post_plate[[1]]) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               Nmedia = .$N0*10^(-t/.$D),
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
  ylab("Estimated microbial concentration (CFU/ml)") +
  xlab("Treatment time (min)")

## Figure 2

plot_grid(p1 + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p2 + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p3 + coord_cartesian(ylim = c(1e-1, 1e8)), 
          p4 + coord_cartesian(ylim = c(1e-1, 1e8)),
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

set.seed(14212)

plate_var_model <- stan("plate_inactivation_variability.stan",
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

plate_var_model
# plot(plate_var_model)
post_plate_var <- As.mcmc.list(plate_var_model)
# pairs(plate_var_model)
# plot(post_plate_var)

## Confidence interval in the media per repetition - Figure 4

set.seed(12412)

# my_c <- wes_palette("Darjeeling1", 3)
my_c <- wes_palette("FantasticFox1", 3)[c(1, 3, 2)]

as.data.frame(post_plate_var[[1]]) %>%
  sample_n(1000) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               Nmedia1 = .$`N0[1]`*10^(-t/.$`D[1]`),
               Nmedia2 = .$`N0[2]`*10^(-t/.$`D[2]`),
               Nmedia3 = .$`N0[3]`*10^(-t/.$`D[3]`),
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
  ylab("Estimated microbial concentration in the media (CFU/ml)") +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  scale_color_manual(values = my_c) +
  coord_cartesian(ylim = c(1e-1, 1e8))

## Variability in the D-value

as.data.frame(post_plate_var[[1]]) %>%
  select(m_logD, s_logD) %>%
  sample_n(10) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(logD = seq(0, 2.5, length = 100),
               p = dnorm(logD, .$m_logD, .$s_logD)
               )
      ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  ggplot() +
  geom_line(aes(x = logD, y = p, colour = sim)) +
  theme(legend.position="none")

set.seed(21442)

p_var_logD <- as.data.frame(post_plate_var[[1]]) %>%
  select(m_logD, s_logD) %>%
  # sample_n(10) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(logD = seq(0.5, 2, length = 100),
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
  geom_ribbon(aes(x = logD, ymin = q05, ymax = q95), alpha = .5,
              fill = my_wes[5]
              # colour = my_wes[5]
              ) +
  theme_bw(base_size = 14) +
  xlab("logarithm of the D-value of strain FBR16 (log min)") +
  ylab("Probability density (·)")

## Prediction interval with logD variability

dist_logD = as.data.frame(post_plate_var[[1]]) %>%
  select(m_logD, s_logD) %>%
  summarize(m_logD = median(m_logD),
            s_logD = median(s_logD))
  
set.seed(241)

p <- rnorm(n = 1000, mean = dist_logD$m_logD, sd = dist_logD$s_logD) %>%
  map(.,
      ~ tibble(t = seq(0, 180, length = 100),
               logN = 8 - t/10^.)
      ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  group_by(t) %>%
  summarize(m_logN = median(logN),
            q05 = quantile(logN, .05),
            q95 = quantile(logN, .95)) %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = m_logN), colour = my_wes[5]) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = .5,
              fill = my_wes[5])



my_sims <- as.data.frame(post_plate_var[[1]]) %>%
  select(m_logD, s_logD) %>%
  sample_n(100) %>%
  mutate(i = row_number()) %>%
  split(.$i) %>%
  map(.,
      ~ tibble(
               logD = rnorm(10, .$m_logD, .$s_logD)
               ) %>%
        mutate(j = row_number()) %>%
        split(.$j) %>%
        map(.,
            ~ tibble(t = seq(0, 180, length = 100),
                     logN = 8 - t/10^.$logD)
            ) %>%
        imap_dfr(., ~ mutate(.x, j = .y))
  ) %>%
  imap_dfr(., ~ mutate(.x, i = .y))

my_sims %>%
  unite(k, c(i, j)) %>%
  ggplot() +
  geom_line(aes(x = t, y = logN, colour = factor(k))) +
  theme(legend.position = "none")

p1 <- my_sims %>%
  group_by(t) %>%
  summarize(m_logN = median(logN),
            q05 = quantile(logN, .05),
            q95 = quantile(logN, .95)) %>%
  geom_ribbon(aes(x = t, ymin = q05, ymax = q95), 
              fill = my_wes[5], alpha = .5, data = .)

p_CI_var <- p + p1 +
  # geom_point(aes(x = time, y = logN, colour = factor(day)),
  #            data = d_var,
  #            inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0, 8)) +
  theme_cowplot() +
  geom_point(aes(x = time, y = log10(count),
                 shape = zero, size = zero), 
             inherit.aes = FALSE, data = d_plot) +
  scale_shape_manual(values = c(1, 16)) +
  scale_size_manual(values = c(3, 3)) +
  xlab("Treatment time (min)") +
  ylab("Log-count in the media (log CFU/ml)") +
  theme(legend.position = "none")

## Figure 5

plot_grid(p_var_logD, p_CI_var,
          nrow = 2,
          labels = "AUTO")

  
  


  