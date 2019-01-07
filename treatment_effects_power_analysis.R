rm(list=ls(all=TRUE))

library(ggplot2)
library(cowplot)
library(pwr)
library(data.table)
library(fitdistrplus)
library(truncdist)

############## define some functions for use in simulations ############

combine_dfs <- function(sim_results_list_item){
  temp_df <- data.frame(treatment=numeric(),
                        deformed=numeric(),
                        undeformed=numeric()
  )
  for(i in 1:length(sim_results_list_item)){
    temp_df <- rbind(temp_df, sim_results_list_item[[i]])
  }
  return(temp_df)
}

model_and_extract <- function(sim_results){
  deformity_results <- as.matrix(sim_results[,c("deformed", "undeformed")])
  treatment <- sim_results[,"treatment"]
  mod <- glm(deformity_results~treatment, family=binomial(link=logit))
  treatment_z <- coef(summary(mod))[6]
  treatment_p <- coef(summary(mod))[8]
  treatment_result <- c(treatment_z, treatment_p)
  names(treatment_result) <- c("z_value", "p_value")
  return(treatment_result)
}

summarise_p_values <- function(z_and_p_val_vector_list, alpha_val=0.05){
  p_vals <- unlist(lapply(z_and_p_val_vector_list, "[[", 2))
  under_alpha <- length(subset(p_vals, p_vals<alpha_val))
  power <- under_alpha/length(z_and_p_val_vector_list)
  output <- c(under_alpha, length(z_and_p_val_vector_list), alpha_val, power)
  names(output) <- c("number_sig_results", "number_simulations", "alpha", "power")
  return(output)
}

# functions for fitting a truncated poisson distribution
# adapted from 
# https://stackoverflow.com/questions/16947799/fitting-a-lognormal-distribution-to-truncated-data-in-r
dtruncated_poisson <- function(x, lambda, number_individuals) {
  dtrunc(x, "pois", a = 1, b = number_individuals, lambda)
}
ptruncated_poisson <- function(x, lambda, number_individuals) {
  ptrunc(x, "pois", a = 1, b = number_individuals, lambda)
}
rtruncated_poisson <- function(n, lambda, number_individuals) {
  rtrunc(n, "pois", a = 1, b = number_individuals, lambda)
}


simulate_experiment <- function(number_simulations,
                                number_beakers,
                                number_individuals,
                                recovery="full",
                                lambda_estimate,
                                treatment_x,
                                treatment_y){
  
  allsim_results_list <- list()
  
  for(number_sim in 1:number_simulations){
    #onesim_results_list <- list()
    all_treatment_results <- list()
    
    for(i in 1:length(treatment_x)){
      temp_x <- treatment_x[i]
      temp_y <- treatment_y[i]
      all_beaker_results <- list()
      
      for(b in 1:number_beakers){
        
        #set number of individuals for this beaker and treatment
         if(recovery=="full"){
           number_ind <- number_individuals
         } else{
        number_ind <- rtruncated_poisson(n=1, lambda=lambda_estimate,
                                         number_individuals=number_individuals)
        }
        
        per_beaker_results <- rbinom(n=number_ind, prob=temp_y, size=1)
        all_beaker_results[[b]] <- per_beaker_results
      }
      
      per_treatment_deformed <- unlist(lapply(all_beaker_results, sum))
      per_treatment_undeformed <- unlist(lapply(all_beaker_results, function(x){length(x)-sum(x)}))
      all_treatment_results[[i]] <- data.frame(treatment=temp_x,
                                               deformed=per_treatment_deformed,
                                               undeformed=per_treatment_undeformed)
    }
      allsim_results_list[[number_sim]] <- all_treatment_results
  }

  # apply model to each item in allsim_results_list
  allsim_results <- lapply(allsim_results_list, combine_dfs)
  model_results <- lapply(allsim_results, model_and_extract)
  power_results <- summarise_p_values(model_results)
  return(power_results)
}

#########################

#set the working directory
setwd("~/Documents/Students_and_Collaborators/bryant_power_analysis/deformities_power_analysis")

####read in data for the overall stress analysis####

data <- (read.csv("combined2.1.csv", header=TRUE))

#read in variables
condition <- data$Condition
mat <- matrix(c(data$Deformed,data$OK), ncol=2)
trans <- data$prop_trans

# calculate number recovered
recovered <- mat[,1]+mat[,2]

# estimating parameters for the distribution of no. individuals recovered
# based on the real data

fitted_poisson <- fitdistrplus::fitdist(recovered, distr="truncated_poisson", start = list(lambda = 8.5)) 
lambda_estimate <- fitted_poisson$estimate
lambda_sd <- fitted_poisson$sd


############ calculate mean deformity proportion in control

control_mean_deformed <- mean(data$Proportion_deformed[which(data$Condition2=="Control")])

############ for Copper experiment

Cdata <- data[ which(data$Copper=='YES'), ]
Cdf <- matrix(c(Cdata$Deformed,Cdata$OK), ncol=2)
Ctrans <- Cdata$prop_trans
Ctreatment <- Cdata$copper_concentration
Clin <- glm(Cdf~Ctreatment, family=binomial(link=logit), data=Cdata) #POWER TEST NEEDED
par(mfrow=c(2,2))
plot(Clin)
plot(Cdata$Proportion_deformed~Cdata$copper_concentration)
treatment_value <- Cdata$copper_concentration
treatment_x <- as.numeric(levels(as.factor(treatment_value)))

# set treatment means so that treatment 4 is twice deformity freq of
# the control (and linear relationship)
slope_calc_veryweak <- control_mean_deformed/treatment_x[5]
treatment_y_veryweak <- slope_calc_veryweak*treatment_x+control_mean_deformed

# set treatment means so that treatment 3 is twice deformity freq of 
# the control (and linear relationship)
slope_calc_weak <- control_mean_deformed/treatment_x[4]
treatment_y_weak <- slope_calc_weak*treatment_x+control_mean_deformed

# set treatment means so that treatment 2 is twice deformity freq of
# control (and linear relationship)
slope_calc_med <- control_mean_deformed/treatment_x[3]
treatment_y_med <- slope_calc_med*treatment_x+control_mean_deformed

# set treatment means so that treatment 1 is twice deformity freq of
# control (and linear relationship)
slope_calc_strong <- control_mean_deformed/treatment_x[2]
treatment_y_strong <- slope_calc_strong*treatment_x+control_mean_deformed

###############
# simulations #
###############

slope_range <- seq(0.1, 1, by=0.05)
power_output <- c()

for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, no death
  treatment_y <- each_slope*treatment_x+control_mean_deformed
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=3,
                                     number_individuals=10,
                                     recovery="full",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
nodeath_3beakers <- cbind(slope_range, power_output)

power_output <- c()
for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, 
  # death according to truncated Poisson distribution modelled from real data
  treatment_y <- each_slope*treatment_x+control_mean_deformed
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=3,
                                     number_individuals=10,
                                     recovery="partial",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
somedeath_3beakers <- cbind(slope_range, power_output)

power_output <- c()

for(each_slope in slope_range){
  # simulation with 4 beakers, 10 ind per beaker, no death
  treatment_y <- each_slope*treatment_x+control_mean_deformed
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=4,
                                     number_individuals=10,
                                     recovery="full",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
nodeath_4beakers <- cbind(slope_range, power_output)

power_output <- c()
for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, 
  # death according to truncated Poisson distribution modelled from real data
  treatment_y <- each_slope*treatment_x+control_mean_deformed
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=4,
                                     number_individuals=10,
                                     recovery="partial",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
somedeath_4beakers <- cbind(slope_range, power_output)

copper_results <- list(nodeath_3beakers, somedeath_3beakers,
                       nodeath_4beakers, somedeath_4beakers)
names(copper_results) <- c("nodeath_3beakers", "somedeath_3beakers",
                         "nodeath_4beakers", "somedeath_4beakers")


par(mfrow=c(1,1))
plot(power_output~slope_range, data=nodeath_3beakers, type="l",
     xlab="Slope (proportion deformed / unit copper)",
     ylab="Power (from 10000 simulations)")
lines(power_output~slope_range, data=somedeath_3beakers, type="l", lty=2)
lines(power_output~slope_range, data=nodeath_4beakers, type="l",
      col="lightblue")
lines(power_output~slope_range, data=somedeath_4beakers, type="l",
      col="lightblue", lty=2)
# lines(power_output~slope_range, data=nodeath_3beakers_20ind, type="l",
#       col="red")
# lines(power_output~slope_range, data=uptohalfdeath_3beakers_20ind, type="l",
#       col="red", lty=2)
# abline(v = slope_calc_veryweak, lty=3)
# text(x=slope_calc_veryweak, y=0.9, labels="Treatment 4", cex=0.8)
# abline(v = slope_calc_weak, lty=3)
# text(x=slope_calc_weak, y=0.7, labels="Treatment 3", cex=0.8)
# abline(v = slope_calc_med, lty=3)
# text(x=slope_calc_med, y=0.2, labels="Treatment 2", cex=0.8)
# abline(v = slope_calc_strong, lty=3)
# text(x=slope_calc_strong, y=0.2, labels="Treatment 1", cex=0.8)
abline(h=0.8, col="red", lty=3)

# Plot of real data including line showing hypothetical correlation
# for which you had 80% power
plot(Cdata$Proportion_deformed~Cdata$copper_concentration,
     pch=19, col=rgb(255, 0, 0, 100, maxColorValue=255))
lines(x=c(0, 1), y=c(0+control_mean_deformed, 0.645+control_mean_deformed))

############
############
# For imidacloprid experiment
############
############

Idata <- data[ which(data$Imidacloprid=='YES'), ]
Idf <- matrix(c(Idata$Deformed,Idata$OK), ncol=2)
Itrans <- Idata$prop_trans
Itreatment <- Idata$imidacloprid_concentration
Ilin <- glm(Idf~Itreatment, family=binomial(link=logit), data=Idata) #POWER TEST NEEDED
par(mfrow=c(2,2))
plot(Ilin)
plot(Idata$Proportion_deformed~Idata$imidacloprid_concentration)
treatment_value <- Idata$imidacloprid_concentration
treatment_x <- as.numeric(levels(as.factor(treatment_value)))

# set treatment means so that treatment 4 is twice deformity freq of
# the control (and linear relationship)
slope_calc_veryweak <- control_mean_deformed/treatment_x[5]
treatment_y_veryweak <- slope_calc_veryweak*treatment_x+control_mean_deformed

# set treatment means so that treatment 3 is twice deformity freq of 
# the control (and linear relationship)
slope_calc_weak <- control_mean_deformed/treatment_x[4]
treatment_y_weak <- slope_calc_weak*treatment_x+control_mean_deformed

# set treatment means so that treatment 2 is twice deformity freq of
# control (and linear relationship)
slope_calc_med <- control_mean_deformed/treatment_x[3]
treatment_y_med <- slope_calc_med*treatment_x+control_mean_deformed

# set treatment means so that treatment 1 is twice deformity freq of
# control (and linear relationship)
slope_calc_strong <- control_mean_deformed/treatment_x[2]
treatment_y_strong <- slope_calc_strong*treatment_x+control_mean_deformed

slope_range <- seq(0.05, 0.5, by=0.05)

for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, no death
  treatment_y <- each_slope*treatment_x+control_mean_deformed
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=3,
                                     number_individuals=10,
                                     recovery="full",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
nodeath_3beakers <- cbind(slope_range, power_output)

power_output <- c()
for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, 
  # death according to truncated Poisson distribution modelled from real data
  treatment_y <- each_slope*treatment_x+control_mean_deformed
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=3,
                                     number_individuals=10,
                                     recovery="partial",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
somedeath_3beakers <- cbind(slope_range, power_output)

power_output <- c()

for(each_slope in slope_range){
  # simulation with 4 beakers, 10 ind per beaker, no death
  treatment_y <- each_slope*treatment_x+control_mean_deformed
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=4,
                                     number_individuals=10,
                                     recovery="full",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
nodeath_4beakers <- cbind(slope_range, power_output)

power_output <- c()
for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, 
  # death according to truncated Poisson distribution modelled from real data
  treatment_y <- each_slope*treatment_x+control_mean_deformed
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=4,
                                     number_individuals=10,
                                     recovery="partial",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
somedeath_4beakers <- cbind(slope_range, power_output)

imidacloprid_results <- list(nodeath_3beakers, somedeath_3beakers,
                       nodeath_4beakers, somedeath_4beakers)
names(imidacloprid_results) <- c("nodeath_3beakers", "somedeath_3beakers",
                         "nodeath_4beakers", "somedeath_4beakers")


par(mfrow=c(1,1))
plot(power_output~slope_range, data=nodeath_3beakers, type="l",
     xlab="Slope (proportion deformed / unit imidacloprid)",
     ylab="Power (from 10000 simulations)")
lines(power_output~slope_range, data=somedeath_3beakers, type="l", lty=2)
lines(power_output~slope_range, data=nodeath_4beakers, type="l",
      col="lightblue")
lines(power_output~slope_range, data=somedeath_4beakers, type="l",
      col="lightblue", lty=2)
abline(h=0.8, col="red", lty=3)
#abline(v=0.18)

# Plot of real data including line showing hypothetical correlation
# for which you had 80% power
plot(Idata$Proportion_deformed~Idata$imidacloprid_concentration,
     pch=19, col=rgb(255, 0, 0, 100, maxColorValue=255))
##CHECK SLOPE INCLUDED HERE!!
lines(x=c(0, 2), y=c(0+control_mean_deformed, 0.36+control_mean_deformed))


############
############
# For food limitation experiment
############
############

Fdata <- data[ which(data$Food=='YES'), ]
Fdf <- matrix(c(Fdata$Deformed,Fdata$OK), ncol=2)
Ftrans <- Fdata$prop_trans
Ftreatment <- Fdata$food_amount
Flin <- glm(Fdf~Ftreatment, family=binomial(link=logit), data=Fdata) #POWER TEST NEEDED
par(mfrow=c(2,2))
plot(Flin)
plot(Fdata$Proportion_deformed~Fdata$food_amount)
treatment_value <- Fdata$food_amount
treatment_x <- as.numeric(levels(as.factor(treatment_value)))

# NB for this experiment, Treatment 6 is really the control (full food)
# and we expect a negative relationship between food and deformity rate

# set treatment means so that treatment 1 is twice deformity freq of
# the control (and linear relationship)
slope_calc_veryweak <- control_mean_deformed*(-1)/(treatment_x[6]-treatment_x[1])
intercept_calc_veryweak <- control_mean_deformed-slope_calc_veryweak*treatment_x[6]
treatment_y_veryweak <- slope_calc_veryweak*treatment_x+intercept_calc_veryweak

# set treatment means so that treatment 2 is twice deformity freq of 
# the control (and linear relationship)
slope_calc_weak <- control_mean_deformed*(-1)/(treatment_x[6]-treatment_x[2])
intercept_calc_weak <- control_mean_deformed-slope_calc_weak*treatment_x[6]
treatment_y_weak <- slope_calc_weak*treatment_x+intercept_calc_weak

# set treatment means so that treatment 4 is twice deformity freq of
# control (and linear relationship)
slope_calc_med <- control_mean_deformed*(-1)/(treatment_x[6]-treatment_x[4])
intercept_calc_med <- control_mean_deformed-slope_calc_med*treatment_x[6]
treatment_y_med <- slope_calc_med*treatment_x+intercept_calc_med

# set treatment means so that treatment 5 is twice deformity freq of
# control (and linear relationship)
slope_calc_strong <- control_mean_deformed*(-1)/(treatment_x[6]-treatment_x[5])
intercept_calc_strong <- control_mean_deformed-slope_calc_strong*treatment_x[6]
treatment_y_strong <- slope_calc_strong*treatment_x+intercept_calc_strong

slope_range <- seq(-0.00002, -0.0002, by=-0.00002)

power_output <- c()
for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, no death
  intercept_calc <- control_mean_deformed-each_slope*treatment_x[6]
  treatment_y <- each_slope*treatment_x+intercept_calc
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=3,
                                     number_individuals=10,
                                     recovery="full",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
nodeath_3beakers <- cbind(slope_range, power_output)

power_output <- c()
for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, 
  # death according to truncated Poisson distribution modelled from real data
  intercept_calc <- control_mean_deformed-each_slope*treatment_x[6]
  treatment_y <- each_slope*treatment_x+intercept_calc
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=3,
                                     number_individuals=10,
                                     recovery="partial",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
somedeath_3beakers <- cbind(slope_range, power_output)

power_output <- c()

for(each_slope in slope_range){
  # simulation with 4 beakers, 10 ind per beaker, no death
  intercept_calc <- control_mean_deformed-each_slope*treatment_x[6]
  treatment_y <- each_slope*treatment_x+intercept_calc
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=4,
                                     number_individuals=10,
                                     recovery="full",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
nodeath_4beakers <- cbind(slope_range, power_output)

power_output <- c()
for(each_slope in slope_range){
  # simulation with 3 beakers, 10 ind per beaker, 
  # death according to truncated Poisson distribution modelled from real data
  intercept_calc <- control_mean_deformed-each_slope*treatment_x[6]
  treatment_y <- each_slope*treatment_x+intercept_calc
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=4,
                                     number_individuals=10,
                                     recovery="partial",
                                     lambda_estimate=lambda_estimate,
                                     treatment_x=treatment_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
somedeath_4beakers <- cbind(slope_range, power_output)

food_results <- list(nodeath_3beakers, somedeath_3beakers,
                             nodeath_4beakers, somedeath_4beakers)
names(food_results) <- c("nodeath_3beakers", "somedeath_3beakers",
                         "nodeath_4beakers", "somedeath_4beakers")

############
############
# For timecourse experiment
############
############

setwd ("~/Documents/Students_and_Collaborators/bryant_power_analysis/deformities_power_analysis")
TCdata <- (read.csv("timecourse_data.csv", header=TRUE))

#read in data
day <- TCdata$Day
OK <- TCdata$OK_animals
def <- TCdata$Deformed_Animals
TCfreq <- TCdata$Proportion_deformed

TCfreq <- TCdata$Proportion_deformed
TCtrans <- TCdata$prop_trans
TCJfreq <- jitter(TCfreq, factor = 6, amount = NULL)
TCmat <- matrix(c(def,OK), ncol=2)

#analysis for day effects
model <- glm(TCmat ~ day, family=binomial(link=logit), data=TCdata) #POWER TEST NEEDED
summary(model)

plot(TCfreq~day)

day_x <- as.numeric(levels(as.factor(day)))
day_ <- 0:4

# calculate number recovered
TCrecovered <- TCmat[,1]+TCmat[,2]

# estimating parameters for the distribution of no. individuals recovered
# based on the real data

TC_fitted_poisson <- fitdistrplus::fitdist(TCrecovered, distr="truncated_poisson", start = list(lambda = 8.5)) 
TC_lambda_estimate <- TC_fitted_poisson$estimate
TC_lambda_sd <- TC_fitted_poisson$sd

Day5_mean_deformity <- mean(TCfreq[which(day==5)])



#set treatment means so that Day 9 is twice deformity freq of
#the control (and linear relationship)
slope_calc_veryweak <- Day5_mean_deformity/(9-5)

# set treatment means so that Day 7 is twice deformity freq of
# control (and linear relationship)
slope_calc_strong <- Day5_mean_deformity/(7-5)

# 
slope_range <- seq(0.03, 0.08, by=0.005)
# 
power_output <- c()
for(each_slope in slope_range){
  # simulation with 4 beakers, 10 ind per beaker, no death
  treatment_y <- each_slope*day_x+Day5_mean_deformity
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=4,
                                     number_individuals=10,
                                     recovery="full",
                                     lambda_estimate=TC_lambda_estimate,
                                     treatment_x=day_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
nodeath_4beakers <- cbind(slope_range, power_output)

power_output <- c()
for(each_slope in slope_range){
  # simulation with 4 beakers, 10 ind per beaker,
  # death according to truncated Poisson distribution modelled from real data
  treatment_y <- each_slope*day_x+Day5_mean_deformity
  print(paste("simulation for slope", each_slope))
  temp_result <- simulate_experiment(number_simulations = 10000,
                                     number_beakers=4,
                                     number_individuals=10,
                                     recovery="partial",
                                     lambda_estimate=TC_lambda_estimate,
                                     treatment_x=day_x,
                                     treatment_y=treatment_y)
  power_output <- c(power_output, temp_result[4])
}
somedeath_4beakers <- cbind(slope_range, power_output)


tc_results <- list(nodeath_4beakers, somedeath_4beakers)
names(tc_results) <- c("nodeath_4beakers", "somedeath_4beakers")


##############################
# plots of all 4 experiments
##############################

pdf(file="Power_analysis_plots_2018-10-29.pdf",
    width=7.5, height=11)
par(mfrow=c(4,2))

plot(power_output~slope_range, data=copper_results$nodeath_3beakers, type="l",
     xlab="Slope (proportion deformed / mg/L copper)",
     ylab="Power")
text(x=0.1, y=0.95, labels="A")
lines(power_output~slope_range, data=copper_results$somedeath_3beakers, type="l", lty=2)
# lines(power_output~slope_range, data=copper_results$nodeath_4beakers, type="l",
#       col="lightblue")
# lines(power_output~slope_range, data=copper_results$somedeath_4beakers, type="l",
#       col="lightblue", lty=2)
abline(h=0.8, col="red", lty=3)


# Plot of real data including line showing hypothetical correlation
# for which you had 80% power
treatment_value <- jitter(Cdata$copper_concentration)
treatment_x <- as.numeric(levels(as.factor(treatment_value)))

plot(Cdata$Proportion_deformed~Cdata$copper_concentration,
     pch=19, col=rgb(255, 0, 0, 100, maxColorValue=255),
     ylab="Proportion deformed",
     xlab="Copper concentration (mg/L)")
text(x=0, y=0.47, labels="B")
intercept_calc <- control_mean_deformed
abline(a = intercept_calc, b = 0.645)


plot(power_output~slope_range, data=imidacloprid_results$nodeath_3beakers, type="l",
     xlab="Slope (proportion deformed / ug/L imidacloprid)",
     ylab="Power")
text(x=0.05, y=0.95, labels="C")
lines(power_output~slope_range, data=imidacloprid_results$somedeath_3beakers, type="l", lty=2)
# lines(power_output~slope_range, data=imidacloprid_results$nodeath_4beakers, type="l",
#       col="lightblue")
# lines(power_output~slope_range, data=imidacloprid_results$somedeath_4beakers, type="l",
#       col="lightblue", lty=2)
abline(h=0.8, col="red", lty=3)

# Plot of real data including line showing hypothetical correlation
# for which you had 80% power
treatment_value <- Idata$imidacloprid_concentration
treatment_x <- as.numeric(levels(as.factor(treatment_value)))

plot(Idata$Proportion_deformed~Idata$imidacloprid_concentration,
     pch=19, col=rgb(255, 0, 0, 100, maxColorValue=255),
     ylab="Proportion deformed",
     xlab="Imidacloprid concentration (ug/L)")
text(x=0, y=0.37, labels="D")
intercept_calc <- control_mean_deformed
abline(a = intercept_calc, b = 0.18)

plot(power_output~slope_range, data=food_results$nodeath_3beakers, type="l",
     xlab="Slope (proportion deformed / uL food)",
     ylab="Power")
lines(power_output~slope_range, data=food_results$somedeath_3beakers, type="l", lty=2)
text(x=-0.0002, y=0.7, labels="E")
abline(h=0.8, col="red", lty=3)
#abline(v=-0.00017)

# Plot of real data including line showing hypothetical correlation
# for which you had 80% power
treatment_value <- Fdata$food_amount
treatment_x <- as.numeric(levels(as.factor(treatment_value)))

plot(Fdata$Proportion_deformed~Fdata$food_amount,
     pch=19, col=rgb(255, 0, 0, 100, maxColorValue=255),
     ylab="Proportion deformed",
     xlab="Food amount (uL)")
text(x=200, y=0.25, labels="F")
plotting_intercept <- control_mean_deformed-(-0.00017)*treatment_x[6]
abline(a = plotting_intercept, b = -0.00017)
#lines(x=c(0, 2), y=c(0+control_mean_deformed, -0.00006+control_mean_deformed))

plot(power_output~slope_range, data=tc_results$nodeath_4beakers, type="l",
     xlab="Slope (proportion deformed / day)",
     ylab="Power")
text(x=0.03, y=0.95, labels="G")
lines(power_output~slope_range, data=tc_results$somedeath_4beakers, type="l", lty=2)
abline(h=0.8, col="red", lty=3)

# Plot of real data including line showing hypothetical correlation
# for which you had 80% power
treatment_value <- TCdata$Day
day_x <- as.numeric(levels(as.factor(treatment_value)))

plot(TCdata$Proportion_deformed~TCdata$Day,
     pch=19, col=rgb(255, 0, 0, 100, maxColorValue=255),
     ylab="Proportion deformed",
     xlab="Day",
     xlim=c(5, 9), ylim=c(0, 1))
text(x=5, y=0.95, labels="H")
intercept_calc <- Day5_mean_deformity
lines(x=c(5, 9), y=c(Day5_mean_deformity, (0.0655*9)+Day5_mean_deformity))

dev.off()
