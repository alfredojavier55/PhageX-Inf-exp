library(SimInf)
## Modelling Phage-X first model ----
## Infection experiment ----
## There is the time of the experiment 60 days ----
tspan <- seq(from = 1,
             to = 60,
             by = 1)

## There are the compartments on the model ----
compartments <- c("Sp", "Ip", "Sa", "Ia")

## There are the propensities between each compartment ----
transitions <- c("Sp -> beta * Sp * (Ip + Ia)/(Sp + Ip + Sa + Ia) -> Ip",      # piglet get infected
                 "Sa -> beta * Sa * (Ip + Ia)/(Sp + Ip + Sa + Ia) -> Ia",      # Adult get infected
                 "Ia -> gamma * Ia -> Sa",           # recovery in adults
                 "Ip -> gamma * Ip -> Sp")           # recovery in piglets

## Parameters
gdata <- c(beta  = 0.2,      # rate of infection
           gamma = 1 / 17)     # animals recover in about 17 days

## Number of nodes pens on the barn ----
n <- 1

## This is the starting condition of the model ----
u0 <- data.frame(Sp = rep(0, n),
                 Ip = rep(0, n),
                 Sa = rep(1, n),
                 Ia = rep(0, n))

## E matrix (select matrix, to do external transfer) ----
E <- matrix(c(1,0,0,0,
              1,0,1,0,
              0,0,1,0,
              1,1,1,1), nrow = 4,
            dimnames = list(c("Sp","Ip","Sa","Ia"), c("births", "susceptibles", "susceptible sows", "all")))

## N matrix (shift) ----
N <- matrix(c(1,0,1,0), nrow = 4,
            dimnames = list(c("Sp","Ip","Sa","Ia"), c("seed infection")))

## Birth piglets on day 10
add1 <- data.frame(event = "enter",
                  time = 10,
                  node = seq_len(n),
                  dest = 0,
                  n = 12,
                  proportion = 0,
                  select = 1, # from the E (Only susceptible piglets)
                  shift = 0 # don't shift them at all
                  )

## Infect the sows on day 7
add2 <- data.frame(event = "intTrans",
                  time = 7,
                  node = seq_len(n),
                  dest = 0,
                  n = 1,
                  proportion = 0,
                  select = 3, # Susceptible sows
                  shift = 1 # shift the sows down according to the 1st column of N
)

events <- rbind(add1,add2)

events
## Build the model ----
model <- mparse(transitions = transitions,
                compartments = compartments,
                gdata = gdata,
                tspan = tspan,
                E = E,
                N = N,
                events = events,
                u0 = u0)

## Run the model ----
result <- run(model)

## Find appropriate parameters for the model given some observed prevalence

## Lets say in the experiment we observe that the prevlance in a pen
## is 50% on 10 days after the piglets are born:
distance <- function(result, ...) {
    abs(prevalence(result, Ia + Ip ~ Ia + Ip + Sa + Ia)$prevalence[20] - 0.50)
}

fit <- abc(model = model,
           priors = c(beta ~ uniform(0, 1)),
           ninit = 1000,
           npart = 200,
           distance = distance)

## This is the estimate of the appropriate beta given the system and the target prevalance (Kind of shit :) )
pdf("fit.pdf")
plot(fit)
dev.off()

## This is a sample from the posterior beta

pdf("posterior.pdf")
plot(run(fit))
dev.off()
