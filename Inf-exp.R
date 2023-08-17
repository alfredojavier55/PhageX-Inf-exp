# Modelling Phage-X first model ----
# Infection experiment ----
# There is the time of the experiment 60 days ----
tspan <- seq(from = 1,
             to = 60,
             by = 1)

## There are the compartments on the model ----
compartments <- c("Sp", "Ip", "Sa", "Ia")

## There are the propensities between each compartment ----
transitions <- c("@ -> mu * Sp -> Sp",              # birth rate of Sp
                 "Sp -> Sp * mu -> @",              # natural-dead of Sp
                 "Ip -> Ip * mu -> @",              # natural-dead of Ip
                 "Sa -> Sa * mu -> @",              # natural-dead of Sa
                 "Ia -> Ia * mu -> @",              # natural-dead of Sa
                 "Sp -> beta * Sp * Ip -> Ip",      # piglet get infected Ip
                 "Sp -> alpha * Sp -> Sa",          # Susceptible piglet becomes susceptible adult
                 "Ip -> alpha * Ip -> Ia",          # Infected piglet becomes Infected adult
                 "Sa -> beta * Sa * Ia -> Ia",      # Adult get infected (horizontal) Ia
                 "Sa -> gamma * Sa -> @",           # Susceptible adult go to slaughter
                 "Ia -> gamma * Ia -> @")           # Infected adult go to slaughter

# There are the parameters ----
gdata <- c(beta  = 1 / 2,
           beta2 = 0.01,   
           alpha = 1 / 70,      # time the piglets become adults
           mu    = 1 / 365*3,   # Simulation is not going to last maybe 6 pregnancies 2.2 pregnancies per year
           gamma = 1 / 150)     # animals go to slaughter

# Number of nodes pens on the barn ----
n <- 5

## This is the starting condition of the model ----
u0 <- data.frame(Sp = rep(10,n),
                 Ip = rep(0, n),
                 Sa = rep(1, n),
                 Ia = rep(0, n))

# E matrix (select matrix, to do external transfer) ----
E <- matrix(c(1,0,0,0,
              1,0,1,0,
              1,1,1,1), nrow = 4, 
            dimnames = list(c("Sp","Ip","Sa","Ia"), c("births", "susceptibles", "all")))

# N matrix (shift) ----
N <- matrix(c(1,0,1,0), nrow = 4, 
            dimnames = list(c("Sp","Ip","Sa","Ia"), c("get-seed infected")))

# Two forms of saving the events: ----
# 1 Movement between compartments (alone, each single event at a time)
add1 <- data.frame(event = "intTrans", 
                  time = 5,
                  node = 1,
                  dest = 0,
                  n = 0,
                  proportion = 0.1, 
                  select = 2, # from the E
                  shift = 1 # from N number of columns with means the states of the compartments
                  )

# 1 Movement between compartments and nodes(alone, each single event at a time)
add2 <- data.frame(event = "extTrans",
                  time = 10,
                  node = 1,
                  dest = 2,
                  n = 0,
                  proportion = 0.5,
                  select = 3,
                  shift = 0
)

events <- rbind(add1,add2)

#2 Movements between compartments two types insted of doing the rbind
add <- data.frame(event = c("intTrans", "extTrans"), 
                  time = c(5, 10),
                  node = c(1,1),
                  dest = c(0,2),
                  n = 0,
                  proportion = c(0.1, 0.5), 
                  select = c(2, 3), # from the E
                  shift = c(1,0) # from N number of columns with means the states of the compartments
)

events <- add

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

# Plot results
plot(result)
head(trajectory(result), n=70)
tail(trajectory(result))
plot(result, range = 0.95)

p <- prevalence(result, Ip+Ia ~ . )
plot(p)
