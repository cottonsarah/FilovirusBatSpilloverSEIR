# Set up ------------------------------------
library(tidyverse)
library(ggplot2)
library(deSolve)


# Define parameters -------------------------------------

# bat rates
beta_bat <- c(0.0045, .45) # paper beta vs what worked for us 
death_bat_juvenile <- 0.00231 # juvenile death
death_bat_adult <- 0.00051 # adult death
aging_bat <- 0.0027 # aging juvenile to adult 
incubation_bat <- c(1/7, 1/21) # incubation rate 
seroconversion_bat <- 0.14 # rate of seroconversion = recovery 

# birth parameters 
birth_timing <- 0 # control birth timing
birth_pulses <- c(0.01, 0.025) # annual birth pulse
birth_scalar <- 0.0041 # scalar control birth rate
birth_synchrony <- 14.3 # annual birth sychrony
carrying_capacity <- 100000 # carrying capacity
birth_offset <-  (2*pi)/365 # not sure where this best defined, used table in citation [23]


# humans params / rates : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4169395/
# based on models of Ebola outbreaks in Guinea, Sierra Leone, and Liberia
beta_human <- 0.333 #human to human infection rate
incubation_human <- 1/5.3 # rate of becoming infectious after exposure
serconversion_human <- 1/5.61 #rate of recovery from infectious
disease_death_human <- 0.643 #case fatality rate
natural_birth_human <- 0.03507 #birth/death rate

interaction_term <- c(0.0085, 0.5006) # how often are fruit bats hunted by humans


# Create coupled modeled for bat populations ------------------------------------------------------

coupled_sir_function <- function(t, y, param){
  
  # all bats 
  beta_b = param["beta_b"] # transmission
  sigma_b = param["sigma_b"] # incubation rate
  tau_b = param["tau_b"] # recovery / seroconversion
  episilon_b = param["episilon_b"] # aging
  
  pulse = param["pulse"] # number of pulses
  
  birth_rate = birth_scalar*sqrt(birth_synchrony/pi)*exp(-birth_synchrony*cos(t*pulse - birth_offset^2))
  
  # juvenile bats 
  death_j = param["death_j"] # death_juvenile
  
  # adult bats
  death_a = param["death_a"] # death_adult
  
  m = param["m"] # migration in example couple...  this will be our biting rate?
  
  #humans
  beta_h = param["beta_h"] # transmission
  sigma_h <- param["sigma_h"] # rate of becoming infectious after exposure
  tau_h <- param["tau_h"] #rate of recovery from infectious
  death_h <- param["death_h"]
  
  # redefining states for BATS
  S1 = y["S1"]
  E1 = y["E1"] # adding exposed compartment 
  I1 = y["I1"]
  R1 = y["R1"]
  
  S2 = y["S2"]
  E2 = y["E2"] # adding exposed compartment 
  I2 = y["I2"]
  R2 = y["R2"]
  
  # defining states for humans
  
  SH = y["SH"]
  EH = y["EH"] # adding exposed compartment 
  IH = y["IH"]
  RH = y["RH"]
  
  
  N_Huma = SH + EH + IH + RH
  N_Bat = S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2
  N_a = S2 + E2+ I2 + R2
  
  # these are the equations for juvenile bats
  dS1 =  -(beta_b * S1 *(I1+I2) ) - (episilon_b*S1) + (birth_rate*N_a*(1-N_Bat/carrying_capacity)) - (death_j*S1*(N_Bat/carrying_capacity)) # -transmission - aging + birth - death
  dE1 = (beta_b * S1 *(I1+I2)) - (episilon_b*E1) - sigma_b*E1 - (death_j*E1*(N_Bat/carrying_capacity)) # +transmission - aging - incubation - death
  dI1 = sigma_b*E1 - (tau_b * I1) - (episilon_b*I1) - (death_j*I1*(N_Bat/carrying_capacity)) # +incubation - recovery - aging - death
  dR1 = (tau_b * I1) - (episilon_b*R1) - (death_j*R1*(N_Bat/carrying_capacity)) # +recovery - aging - death
  
  # these are the equations for adult bats
  dS2 = -(beta_b*S2*(I1+I2)) + (episilon_b*S1) - (death_a*S2*(N_Bat/carrying_capacity)) # -transmission + aging - death
  dE2 = (beta_b*S2*(I1+I2)) + (episilon_b*E1) - (sigma_b*E2) - (death_a*E2*(N_Bat/carrying_capacity)) # +transmission + aging - incubation - death
  dI2 = (sigma_b*E2) - (tau_b * I2) + (episilon_b*I1) - (death_a*I2*(N_Bat/carrying_capacity)) # +incubation - recovery + aging - death
  dR2 = (tau_b * I2) + (episilon_b*R1) - (death_a*R2*(N_Bat/carrying_capacity)) # +recovery + aging - death
  
  # these are the equations for adult humans 
  dSH = -(beta_h*SH*(IH+(m*I2))) + natural_birth_human*N_Huma - natural_birth_human*SH
  dEH = (beta_h*SH*(IH+(m*I2))) - (sigma_h*EH) - natural_birth_human*EH
  dIH = (sigma_h*EH) - (tau_h*IH) - natural_birth_human*IH
  dRH = (1-death_h)*(tau_h*IH) - natural_birth_human*RH 
  
  # returning list of current system of eq for time t
  return(list(c(dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2, dSH, dEH, dIH, dRH))) 
}


run_coupled_sir_model <- function (beta_b, sigma_b, tau_b, episilon_b, death_j, death_a, m, beta_h, sigma_h, tau_h, death_h, pulse, initial.state, max.time, freq.dependent) { # will need to add human params here
  
  # account for freq dependent transmission
  beta_divisor_1 <- ifelse(freq.dependent == TRUE, initial.state["S1"] + initial.state["E1"] + initial.state["I1"] + initial.state["R1"] + initial.state["S2"] + initial.state["E2"] + initial.state["I2"] + initial.state["R2"], 1)
  
  beta_divisor_h <- ifelse(freq.dependent == TRUE, initial.state["SH"] + initial.state["EH"] + initial.state["IH"] + initial.state["RH"], 1)
  
  
  # add params to list 
  param <- c(beta_b = beta_b/beta_divisor_1, sigma_b = sigma_b, tau_b = tau_b, episilon_b = episilon_b, 
             death_j = death_j, death_a = death_a, m = m, beta_h=beta_h/beta_divisor_h, sigma_h=sigma_h, 
             tau_h=tau_h, death_h=death_h, pulse = pulse)#  #beta_human = beta_human/beta_divisor_2, # K = K,   
  
  
  
  times <- seq(0, max.time, 1)
  
  # solving system of equations
  coupled_sir_output <- deSolve::lsoda(initial.state, times, coupled_sir_function, param)
  
  return(as.data.frame(coupled_sir_output))
}

# Set model paramaters -----

# bat terms 
tau_b <- seroconversion_bat
episilon_b <- aging_bat

death_j <- death_bat_juvenile
death_a <- death_bat_adult


# human terms
m<- interaction_term
beta_h <- beta_human
sigma_h <- incubation_human
tau_h <- serconversion_human
death_h <- disease_death_human
pulse <- birth_pulses[[1]]


freq.dependent <- TRUE
max.time <- 5000
initial.state <- c(S1=20000, E1 = 0, I1=500, R1=0, S2=20000, E2 = 0, I2=500, R2=0, SH = 10000, EH = 0, IH = 0, RH = 0)

# Run model, no stochastic -----------------
combined_simulations <- data.frame()

for (k in interaction_term){
  for(x in beta_bat){
    for(b in birth_pulses){
      for(i in incubation_bat){
        beta_b <- x
        sigma_b <- i
        pulse <- b
        m <- k
        
        simu <- run_coupled_sir_model(beta_b, sigma_b, tau_b, episilon_b, death_j, death_a, m,beta_h, sigma_h,tau_h,death_h,pulse, initial.state, max.time, freq.dependent) #birth_rate, death_j, death_a, K, episilon_b,
        
        simu <- simu %>% 
          mutate(incubation_period = sigma_b,
                 birth_pulse = as.factor(b),
                 bat_transmission = x, 
                 interaction = k,
                 N = S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2)
        
        combined_simulations <- rbind(combined_simulations, simu)
        
      }
    }
  }
}

# clean up code and exploratory visualization 
combined_simulations <- combined_simulations %>%
  mutate(birth_pulse = ifelse(birth_pulse == (0.01), "1 birth pulse / year", "2 birth pulses / year"),
         bat_transmission = ifelse(bat_transmission == (0.0045), "original beta", "high beta"),
         incubation_period = ifelse(incubation_period == (1/7), "7 days", "21 days"),
         interaction = ifelse(interaction==(0.0085), "low contact", 'high contact'))

combined_simulations %>%  
  ggplot(aes(x = time)) +
  geom_line(aes(y = S1 + S2), color = "gold") +
  geom_line(aes(y = E1 + E2), color = "purple") +
  geom_line(aes(y = I1 + I2), color = "red") +
  geom_line(aes(y = R1 + R2), color = "green") +
  geom_line(aes(y = N), color = "black") +
  labs(title = "SEIR Over Time in Bat Populations",
       subtitle = "Varying birth pulse and transmission rate",
       x = "Time",
       y = "Population",
       color = "Legend")+ 
  facet_wrap(~birth_pulse*bat_transmission*interaction, nrow = 3, scales = "free_y")


# Run simulations with stochastic transmission in humans ------------

# Number of simulations, this takes a long time to run FYI 
num_simulations <- 100

simulations_for_tracking_all <- data.frame()

for (k in interaction_term){
  for(x in beta_bat){
    for(b in birth_pulses){
      for(i in incubation_bat){
        for (n in 1:num_simulations) {
    
            beta_b <- x # assign per each simulation 
            sigma_b <- i
            pulse <- b
            
            
            # add stochastic elements
            shape1 <- 3 
            shape2 <- shape1 * (1/beta_human - 1)
            
            stochastic_beta <- rbeta(1, shape1, shape2) # human beta is stochastic follows mean of beta_human
            
            beta_h <- stochastic_beta
            
            simu <- run_coupled_sir_model(beta_b, sigma_b, tau_b, episilon_b, death_j, death_a, m, beta_h, sigma_h,tau_h,death_h,pulse, initial.state, max.time, freq.dependent) #birth_rate, death_j, death_a, K, episilon_b,
            
            simu <- simu %>% 
              mutate(incubation_period = sigma_b,
                     birth_pulse = as.factor(b),
                     bat_transmission = x, 
                     interaction = k,
                     Bat_N = S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2,
                     Human_N = SH + IH + EH + RH,
                     human_transmission = stochastic_beta,
                     simulation = n)
            
            # Store the results of the simulation
            simulations_for_tracking_all <- rbind(simulations_for_tracking_all, simu)
        }
      }
    }
  }
}
      

#write_parquet(simulations_for_tracking_all, "final-project/data/simulations_100_master.parquet")

# tidy up 
simulations_for_tracking_all <- simulations_for_tracking_all %>%
  mutate(birth_pulse = ifelse(birth_pulse == (0.01), "1 birth pulse / year", "2 birth pulses / year"),
         bat_transmission = ifelse(bat_transmission == (0.0045), "original beta", "high beta"),
         incubation_period = ifelse(incubation_period == (1/7), "7 days", "21 days"),
         interaction = ifelse(interaction==(0.0085), "low contact", 'high contact')
         )

############# PART 2 creating data visualizations -----------------


# set up ----------------------------
library(arrow)
library(ggplot2)
library(tidyverse)

# read in data ---------------------------

# insert csv w/ local directory or uncomment sim line below 
sims <- simulations_for_tracking_all
# or load data below 
#sims <- read_csv("data/simulation_results.csv")

sims_analysis <- sims %>% 
  mutate(S_Bat = S1 + S2, 
         E_Bat = E1 + E2,
         I_Bat = I1 + I2,
         R_Bat = R1 + R2,
         time_period = ifelse(time < 1000, "First Wave", "Secondary Waves")) %>% 
  select(-`...1`)

# human infections over time
sims_analysis %>% 
  ggplot(aes(x = time, color = human_transmission)) +
  geom_line(aes(y = IH)) +
  labs(title = "SEIR Over Time in Human Populations",
       subtitle = "Varying birth rate and incubation period",
       x = "Time",
       y = "Population",
       color = "Legend")+ 
  facet_wrap(~birth_pulse*bat_transmission*incubation_period, scales = "free_y") 



# bat population varying birth pulse 

fig1 <- sims_analysis %>%  
  filter(bat_transmission == "original beta", 
         incubation_period == "7 days") %>% 
  select(time, Bat_N, birth_pulse) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = Bat_N, color = Bat_N)) +
  labs(title = "Bat Population Over Time by Birth Pulses [1,2]",
       x = "Time",
       y = "Population Size",
       color = "Population Growth")+ 
  facet_wrap(~birth_pulse, nrow = 2, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradient2(low = "darkgreen", mid = "darkorange", high = "#C41E3A", midpoint = 70000,
                        labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw() 


# infected and susceptible bats over time

fig2 <- sims_analysis %>% 
  filter(incubation_period == "7 days",
         birth_pulse == "1 birth pulse / year") %>% 
  select(time, S1, I1, I2, bat_transmission) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = S1, linetype = "Susceptible Juveniles"), color = "darkblue") +
  geom_line(aes(y = I1, linetype = "Infected Juveniles"), color = "#C41E3A") +
  geom_line(aes(y = I2, linetype = "Infected Adults"), color = "darkorange") +
  labs(title = "Bat Infections and Susceptibility Over Time by Bat Transmission",
       subtitle = "Assuming 7 day incubation period and 1 birth pulse / year",
       x = "Time",
       y = "Population Size",
       color = "Population Growth")+ 
  facet_wrap(~bat_transmission, nrow = 2, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(labels = scales::comma) +
  scale_linetype_manual(name = "Line Type", 
                        values = c("solid", "solid", "dashed"), 
                        labels = c("Infected Adults", "Infected Juveniles", "Susceptible Juveniles")) +
  theme_bw() 

# infected young and adult bats over time 
fig3 <- sims_analysis %>% 
  filter(incubation_period == "21 days",
         birth_pulse == "1 birth pulse / year") %>% 
  select(time, S1, I1, I2, bat_transmission) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = I1, linetype = "Infected Juveniles"), color = "#C41E3A") +
  geom_line(aes(y = I2, linetype = "Infected Adults"), color = "darkorange") +
  labs(title = "Bat Infections Over Time by Bat Transmission",
       subtitle = "Assuming 21 day incubation period and 1 birth pulse / year",
       x = "Time",
       y = "Population Size",
       color = "Population Growth")+ 
  facet_wrap(~bat_transmission, nrow = 2, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_continuous(labels = scales::comma) +
  scale_linetype_manual(name = "Line Type", 
                        values = c("solid", "solid"), 
                        labels = c("Infected Adults", "Infected Juveniles")) +
  theme_bw() 


#### visualizing spillover 

all_combinations <- expand.grid(
  time_period = unique(sims_analysis$time_period),
  incubation_period = unique(sims_analysis$incubation_period),
  birth_pulse = unique(sims_analysis$birth_pulse),
  bat_transmission = unique(sims_analysis$bat_transmission)
)

# Calculate spillover percentages

# create spillover summary dataset by params
spillover <- sims_analysis %>%
  filter(IH > 10) %>%
  group_by(time_period, incubation_period, birth_pulse, bat_transmission) %>%
  summarise(spillover_percent = ifelse(n_distinct(simulation) > 1, (n_distinct(simulation)/100)*100, 0)) %>%
  right_join(all_combinations, by = c("time_period", "incubation_period", "birth_pulse", "bat_transmission")) %>%
  mutate(spillover_percent = ifelse(is.na(spillover_percent), 0, spillover_percent))

# model human infections depending on bat transmissions 
fig4 <- sims_analysis %>% 
  select(time, S1, IH, EH, bat_transmission, time_period) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = IH, linetype = "Infected Humans", color = IH)) +
  labs(title = "Humans Infections Over Time by Bat Transmission",
       x = "Time",
       y = "Infections")+ 
  facet_wrap(~time_period*bat_transmission, nrow = 2, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradient2(low = "darkgreen", mid = "darkorange", high = "#C41E3A", midpoint = 800) +
  scale_linetype_manual(name = "Line Type", 
                        values = c("solid"), 
                        labels = c("Infected Humans")) +
  theme_bw() 


# summarize variability of spillover into human pops 
fig5 <- spillover %>%  
  mutate(incubation_period = ifelse(incubation_period == "21 days", "21 day incubation", "7 day incubation")) %>% 
  ggplot(aes(x = as.factor(bat_transmission), y = spillover_percent, fill = spillover_percent)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "forestgreen", high = "#C41E3A", labels = c("0%", "25%", "50%", "75%", "100%")) +
  facet_wrap(~time_period*birth_pulse*incubation_period, nrow = 2) + 
  labs(x = "Bat Transmission Values", y = "% of Simulations", fill = "% of Simulations with Spillover",
       subtitle = "Varying Birth Pulse & Incubation Period") +
  ggtitle("Variability in Spillover into Human Populations (Infected > 10) in First & Secondary Waves") +
  theme_bw()


# variability in bat transmission by params 
fig8 <- sims_analysis %>% 
  mutate(incubation_period = ifelse(incubation_period == "21 days", "21 day incubation", "7 day incubation")) %>% 
  filter(simulation == 1) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = I_Bat, linetype = "Infected Bats", color = I_Bat)) +
  labs(title = "Infections Over Time in Bat Populations",
       subtitle = "Varying transmission rate, birth pulse, and incubation period",
       x = "Time",
       y = "Infections",
       color = "Legend")+ 
  facet_wrap(~bat_transmission*birth_pulse*incubation_period, nrow = 2, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_color_gradient2(low = "darkgreen", mid = "darkorange", high = "#C41E3A", midpoint = 7500) +
  scale_linetype_manual(name = "Line Type", 
                        values = c("solid"), 
                        labels = c("Infected Bats")) +
  theme_bw() 




