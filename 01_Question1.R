# Question 1----------------------------------------
library(tidyverse)

gety <- function(seed, 
                 n = 15) {
  set.seed(seed)
  X <- runif(n, min = 0, max = 1)
  beta <- 2 * exp(-2 * X) 
  gamma <- 2.2
  
  y1 <- map(gamma * beta, ~rpois(1, .)) %>% unlist()
  y2 <- map(beta, ~rpois(1, .)) %>% unlist()
  
  simy <- cbind(y1 = y1, y2 = y2)
  return(simy)
}

## simulate 5000 y1 y2 set.
sim_num <- 5000
y_ss15 <- map(1:sim_num, ~gety(seed = ., n = 15))
y_ss100 <- map(1:sim_num, ~gety(seed = ., n = 100))
save(y_ss15, y_ss100, file = "bios_qexam_data_1.Rdata")

y1 <- map(y_ss15, ~mean(.[, 1])) %>% unlist()
y2 <- map(y_ss15, ~mean(.[, 2])) %>% unlist()
mean(y1); mean(y2); mean(y1) / mean(y2); sqrt((mean(y1) + mean(y2))^3 * mean(y2) / (mean(y1)^3 * 15))
mean(y1) / mean(y2) + 1.96 * sqrt((mean(y1) + mean(y2))^3 * mean(y2) / (mean(y1)^3 * 15))
mean(y1) / mean(y2) - 1.96 * sqrt((mean(y1) + mean(y2))^3 * mean(y2) / (mean(y1)^3 * 15))
power.t.test(n = 15, delta = 2.2 -1, sd = 0.421, sig.level = 0.05)

y1 <- map(y_ss100, ~mean(.[, 1])) %>% unlist()
y2 <- map(y_ss100, ~mean(.[, 2])) %>% unlist()
mean(y1); mean(y2); mean(y1) / mean(y2); sqrt((mean(y1) + mean(y2))^3 * mean(y2) / (mean(y1)^3 * 100))
mean(y1) / mean(y2) + 1.96 * sqrt((mean(y1) + mean(y2))^3 * mean(y2) / (mean(y1)^3 * 100))
mean(y1) / mean(y2) - 1.96 * sqrt((mean(y1) + mean(y2))^3 * mean(y2) / (mean(y1)^3 * 100))
power.t.test(n = 100, delta = 1.2, sd = 0.163, sig.level = 0.05)


s15 <- rnorm(n = 5000, mean = 2.2, sd = 0.163)
s100 <- rnorm(n = 5000, mean = 2.2, sd = 0.163)
ggplot() +
  geom_freqpoly(aes(s15), color = "orange") +
  geom_freqpoly(aes(s100), color = "red") +
  geom_vline(xintercept = 1,
             linetype="dashed", 
             color = "black", 
             size = 1) +
  ylab("gamma") + 
  theme_bw()

load("bios_qexam_data_1.Rdata")

## function for question 1d
get_gamma <- function(y,
                      alpha = 0.05,
                      gamma_min = 0,
                      gamma_max = 5.5,
                      gamma_step = 0.01,
                      N = 2000) {
  y <- as.data.frame(y)
  ## sum of y1 y2
  y1 <- y[, 1]
  y2 <- y[, 2]
  y12 <- y1 + y2
  gamma_int <- seq(gamma_min, gamma_max, by = gamma_step)
  gamma_hat <- sum(y1) / sum(y2)
  ## emptry matrix to store the values
  gamma_mle <- matrix(NA, nrow = N, ncol = length(gamma_int))
  
  for (j in seq_along(1:N)) {
    for (i in seq_along(1:length(gamma_int))) {
      phat <- gamma_int[i] / (1 + gamma_int[i])
      ## simulate y1 | y1+y2 2000 times
      # ya3 <- map(y12, ~rbinom(n = 1, size = ., 
      #                         prob = phat)) %>% unlist()
      
      ya3 <- rbinom(length(y1),
                    size = unlist(y1 + y2),
                    prob = phat)
      gamma_mle[j, i] <- sum(ya3) / sum(y12 - ya3) 
    }
  }
  
  ## get F(0.25) and F(0.975)
  q025 <- apply(gamma_mle, 2, 
                quantile, probs = alpha / 2,
                na.rm = TRUE)
  
  q975 <- apply(gamma_mle, 2, 
                quantile, probs = 1 - alpha / 2, 
                na.rm = TRUE)
  ## check whether gamma_int is located inside F025 F975
  gamma_ind <- ifelse((gamma_hat > q025) & (gamma_hat < q975), 1, 0)
  ## find the smallest and largest gamma_int
  lower <- min(gamma_int[which(gamma_ind == 1)])
  upper <- max(gamma_int[which(gamma_ind == 1)])
  return(c(lower, upper))
}

get_gamma(y_ss100[[1]])

## parallel programming
library(furrr)
library(parallel)
detectCores()

## use 6 cores
plan(multisession, workers = 6)
finall_5k_N15 <- future_map_dfc(y_ss15, ~get_gamma(y = .))
finall_5k_N100 <- future_map_dfc(y_ss100, ~get_gamma(y = .))

save(finall_5k_N15, finall_5k_N100, finall_per,
     file = "bios_qexam_question1.Rdata")

load("bios_qexam_question1.Rdata")
library(matrixStats)
rowMeans(as.matrix(finall_5k_N15))
rowMeans(as.matrix(finall_5k_N100))

library(tidyverse)
y <- read_csv("data/MVPA_intervention.csv")
y1 <- y$Y1
y2 <- y$Y2
y1_bar <- mean(y1)
y2_bar <- mean(y2)
## mle 
(mean <- y1_bar / y2_bar)
## sd 
## use observed information
(gamma <- y1_bar / y2_bar)
Inof <- ((2 * gamma + 1) / (gamma^2 * (1 + gamma)^2) * sum(y1) - 1 / (gamma + 1)^2 * sum(y2)) / 15
(sd <- (1 / Inof) %>% sqrt()) / sqrt(15); sd
## use phat
sd0 <- sqrt((y1_bar + y2_bar) * y1_bar / y2_bar^3) / sqrt(15)
## 95% CI
mean - 1.96 * sd0
mean + 1.96 * sd0 

get_gamma(y = y, gamma_max = 15)

## use permutation to get the confidence interval 
y_permu <- map(1:1000, ~sample_n(y, size = 15, replace = TRUE))

## parallel programming
library(furrr)
library(parallel)

## use 6 cores
plan(multisession, workers = 6)
finall_per <- future_map_dfc(y_permu, ~get_gamma(y = ., gamma_max = 10))
rowMedians(as.matrix(finall_per))
## [1] 1.67 9.95
## [1] 1.825 8.640