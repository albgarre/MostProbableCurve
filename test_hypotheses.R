
library(tidyverse)
library(cowplot)

## Random walk within an experiment/chance

simulate_1_exp <- function(delta_t, D, N0, max_time) {
  
  t <- seq(0, max_time, by = delta_t)
  red_factor <- -delta_t/D
  
  N <- matrix(NA, nrow = length(t), ncol = 1)
  N[1] <- N0
  
  for (i in 1:(length(t) - 1)) {  # This is painful
    
    Ni <- N[i]
    lambda <- Ni*10^red_factor
    Nip1 <- rpois(1, lambda)
    N[i+1] <- Nip1
    
  }
  
  as.double(N)
  
}

delta_t <- 2
D <- 2
N0 <- 1e8
max_time <- 16

n_sims <- 50

set.seed(12412)

my_sims <- c(1:n_sims) %>%
  map_dfc(.,
      ~ simulate_1_exp(delta_t, D, N0, max_time)
      ) %>%
  set_names(., paste0("sim_", 1:n_sims)) %>%
  mutate(t = seq(0, max_time, by = delta_t)) %>%
  select(t, everything())

my_sims %>%
  pivot_longer(-t, names_to = "sim", values_to = "N") %>%
  mutate(logN = log10(N)) %>%
  ggplot(aes(x = t, y = logN, colour = sim)) +
  geom_point(size = 1) +
  geom_line(size = .1) +
  theme(legend.position = "none")

#- Effect of delta t

D <- 2
N0 <- 1e8
max_time <- 16

n_sims <- 10000

set.seed(12412)

my_sims_d2 <- c(1:n_sims) %>%
  map_dfc(.,
          ~ simulate_1_exp(delta_t = 2, D, N0, max_time)
  ) %>%
  set_names(., paste0("sim_", 1:n_sims)) %>%
  mutate(t = seq(0, max_time, by = 2)) %>%
  select(t, everything()) %>%
  pivot_longer(-t, names_to = "sim", values_to = "N") %>%
  mutate(logN = log10(N))

my_sims_d1 <- c(1:n_sims) %>%
  map_dfc(.,
          ~ simulate_1_exp(delta_t = 1, D, N0, max_time)
  ) %>%
  set_names(., paste0("sim_", 1:n_sims)) %>%
  mutate(t = seq(0, max_time, by = 1)) %>%
  select(t, everything()) %>%
  pivot_longer(-t, names_to = "sim", values_to = "N") %>%
  mutate(logN = log10(N))

my_sims_d4 <- c(1:n_sims) %>%
  map_dfc(.,
          ~ simulate_1_exp(delta_t = 4, D, N0, max_time)
  ) %>%
  set_names(., paste0("sim_", 1:n_sims)) %>%
  mutate(t = seq(0, max_time, by = 4)) %>%
  select(t, everything()) %>%
  pivot_longer(-t, names_to = "sim", values_to = "N") %>%
  mutate(logN = log10(N))

my_sims_d8 <- c(1:n_sims) %>%
  map_dfc(.,
          ~ simulate_1_exp(delta_t = 8, D, N0, max_time)
  ) %>%
  set_names(., paste0("sim_", 1:n_sims)) %>%
  mutate(t = seq(0, max_time, by = 8)) %>%
  select(t, everything()) %>%
  pivot_longer(-t, names_to = "sim", values_to = "N") %>%
  mutate(logN = log10(N))

my_sims_d0p2 <- c(1:n_sims) %>%
  map_dfc(.,
          ~ simulate_1_exp(delta_t = .2, D, N0, max_time)
  ) %>%
  set_names(., paste0("sim_", 1:n_sims)) %>%
  mutate(t = seq(0, max_time, by = .2)) %>%
  select(t, everything()) %>%
  pivot_longer(-t, names_to = "sim", values_to = "N") %>%
  mutate(logN = log10(N))

list(`delta = 1` = my_sims_d1,
     `delta = 2` = my_sims_d2,
     `delta = 4` = my_sims_d4,
     `delta = 8` = my_sims_d8,
     `delta = .2` = my_sims_d0p2
     ) %>%
  map(., ~ filter(., t == max_time)) %>%
  imap_dfr(., ~ mutate(.x, delta = .y)) %>%
  ggplot() +
  geom_histogram(aes(N, fill = delta), position = "dodge")
  
## Error in the plating

simulate_plating <- function(delta_t, D, N0, max_time, 
                             dil_f, n_dils, plated_v,
                             n_plates = 1) {
  
  N <- simulate_1_exp(delta_t, D, N0, max_time)
  
  lambdas <- N*dil_f^n_dils*plated_v
  
  c(1:n_plates) %>%
    map(.,
        ~tibble(colonies = rpois(length(N), lambdas),
                lambda = lambdas,
                N_obs = colonies/plated_v/(dil_f^n_dils))
        )
  
}

delta_t <- 2
D <- 2
N0 <- 1e8
max_time <- 16

dil_f <- .1
n_dils <- c(6, 5, 4, 3, 2, 1, 0, 0, 0)
plated_v <- .1
n_plates <- 10

# n_sims <- 50
set.seed(12412)

simulate_plating(delta_t, D, N0, max_time,
                 dil_f, n_dils, plated_v,
                 n_plates
                 ) %>%
  map(., ~ mutate(., 
                  t = seq(0, max_time, by = delta_t),
                  logNobs = log10(N_obs)),
      ) %>%
  imap_dfr(., ~ mutate(.x, sim = .y)) %>%
  ggplot(aes(x = t, y = logNobs, colour = factor(sim))) +
  geom_point(size = 1) +
  geom_line(size = .2) +
  theme(legend.position = "none")

## Combine both

delta_t <- 2
D <- 2
N0 <- 1e8
max_time <- 16

dil_f <- .1
n_dils <- c(6, 5, 4, 3, 2, 1, 0, 0, 0)
plated_v <- .1
n_plates <- 3

n_sims <- 5
set.seed(12412)

c(1:n_sims) %>%
  map(., 
      ~ simulate_plating(delta_t, D, N0, max_time,
                         dil_f, n_dils, plated_v,
                         n_plates
                         )
      ) %>%
  map(.,
      ~ map(., ~ mutate(., 
                        t = seq(0, max_time, by = delta_t),
                        logNobs = log10(N_obs)),
            )
      ) %>%
  map(.,
      ~ imap_dfr(., ~ mutate(.x, plate = as.character(.y)))
      ) %>%
  imap_dfr(.,
           ~ mutate(.x, rep = as.character(.y))
           ) %>%
  ggplot(aes(x = t, y = logNobs, colour = rep)) +
  geom_point(aes(shape = plate), size = 1) +
  geom_line(aes(linetype = plate), size = .2) +
  theme(legend.position = "none")

## Resistance varying between experiments

delta_t <- 4
mean_logD <- log10(2)
sd_logD <- .05
N0 <- 1e8
max_time <- 16

dil_f <- .1
n_dils <- c(6, 4, 2, 0, 0)
plated_v <- .1
n_plates <- 3

n_sims <- 5
# set.seed(12412)

rnorm(n_sims, mean_logD, sd_logD) %>%
  map(., 
      ~ simulate_plating(delta_t, D=10^., N0, max_time,
                         dil_f, n_dils, plated_v,
                         n_plates
      )
  ) %>%
  map(.,
      ~ map(., ~ mutate(., 
                        t = seq(0, max_time, by = delta_t),
                        logNobs = log10(N_obs)),
      )
  ) %>%
  map(.,
      ~ imap_dfr(., ~ mutate(.x, plate = as.character(.y)))
  ) %>%
  imap_dfr(.,
           ~ mutate(.x, rep = as.character(.y))
  ) %>%
  ggplot(aes(x = t, y = logNobs, colour = rep)) +
  geom_point(aes(shape = plate), size = 1) +
  geom_line(aes(linetype = plate), size = .2) +
  theme(legend.position = "none")

## N0 varying between experiments

delta_t <- 4
mean_logD <- log10(2)
sd_logD <- .05
mean_logNmax <- 10
sd_logNmax <- .1
dil_N0 <- 1e-2
max_time <- 16

dil_f <- .1
n_dils <- c(6, 4, 2, 0, 0)
plated_v <- .1
n_plates <- 3

n_sims <- 5
set.seed(12412)

tibble(logNmax = rnorm(n_sims, mean_logNmax, sd_logNmax),
       logD = rnorm(n_sims, mean_logD, sd_logD),
       N0 = rpois(n_sims, 10^logNmax*dil_N0)
       ) %>%
  mutate(rep = row_number()) %>%
  split(.$rep) %>%
  map(., 
      ~ simulate_plating(delta_t, D=10^.$logD, .$N0, max_time,
                         dil_f, n_dils, plated_v,
                         n_plates
      )
  ) %>%
  map(.,
      ~ map(., ~ mutate(., 
                        t = seq(0, max_time, by = delta_t),
                        logNobs = log10(N_obs)),
      )
  ) %>%
  map(.,
      ~ imap_dfr(., ~ mutate(.x, plate = as.character(.y)))
  ) %>%
  imap_dfr(.,
           ~ mutate(.x, rep = as.character(.y))
  ) %>%
  ggplot(aes(x = t, y = logNobs, colour = rep)) +
  geom_point(aes(shape = plate), size = 1) +
  geom_line(aes(linetype = plate), size = .2) +
  theme(legend.position = "none")


 



