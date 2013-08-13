#### Finish the analysis by linking the CV_precip to the elasticity to variance 

rm(list = ls())

setwd("~/Documents/Courses/Matrix Population Models/final_project/analyses/")
source('Ellis_et_al_deterministic_analysis.R')
source('Ellis_et_al_stochastic_analysis.R')
source('find_cv_precip_for_Ellis_et_al_pops.R')

species_info = read.table('Ellis_et_al_Species_Information.txt', sep = '\t', header = T)


names(species_info)

species_info
deterministic_analysis = read.csv('deterministic_analysis_results.csv', header = T)
stochastic_analysis = read.csv('stochastic_analysis_results.csv', header = T)
precip_cv_data = read.csv('prec_cv_data.csv')


all_data = cbind(stochastic_analysis, deterministic_analysis, prec_cv = precip_cv_data$cv_vals)

all_data = all_data[, -11]

seed_bank = species_info[, c(1, 18)]
seed_bank

all_data = merge(all_data,seed_bank, by.all = "SPP")
all_data

all_data$E_var_proportion = abs(all_data$sum_E_var)/(all_data$sum_E_mu + abs(all_data$sum_E_var))
all_data$SE_var_proportion = abs(all_data$sum_SE_var)/(all_data$sum_E_mu + abs(all_data$sum_E_var))
all_data$FE_var_proportion = abs(all_data$sum_FE_var)/(all_data$sum_E_mu + abs(all_data$sum_E_var))


summary(all_data$E_var_proportion)
hist(all_data$E_var_proportion)
dev.copy2pdf(file = "hist_Elasticity_variance.pdf")

summary(all_data$sum_SE_var)
hist(all_data$sum_SE_var)
dev.copy2pdf(file = "hist_Elasticity_survival_variance.pdf")

summary(all_data$sum_FE_var)
hist(all_data$sum_FE_var)
dev.copy2pdf(file = "hist_Elasticity_fecundity_variance.pdf")

hist(all_data$sum_E_mu)
hist(all_data$prec_cv)

summary(all_data$Lambda)
hist(all_data$Lambda)
summary(all_data$prec_cv)

summary(all_data$Lambda_s)
summary(all_data$Lambda)
######### GEN FIGURE 1

xlab1 = "CV of Precipitation (%)"
ylab1 = "Proportion Elasticity to Variance"
ylab2 = "Proportion Elasticity to Variance of Survival"
ylab3 = "proportion Elasticity to Variance of Fecundity"

par (mfrow = c(3,1))
#### panel a
plot(E_var_proportion ~ prec_cv, data = all_data, xlab = xlab1, ylab = ylab1)
cortesta = lm(all_data$E_var_proportion ~ all_data$prec_cv)
summary(cortesta)

#### panel b
plot(SE_var_proportion ~ prec_cv, data = all_data, xlab = xlab1, ylab = ylab2)
cortestb = lm(all_data$SE_var_proportion ~ all_data$prec_cv)
summary(cortestb)

#### panel c
plot(FE_var_proportion ~ prec_cv, data = all_data, xlab = xlab1, ylab = ylab3)
cortestc = lm(all_data$FE_var_proportion ~ all_data$prec_cv)
summary(cortestc)
dev.copy2pdf(file = "figure_1.pdf")


#### Gen Figure 2
par (mfrow = c(1,1))
plot(log(LE2) ~ prec_cv, data = all_data, type = 'n', xlab = xlab1, ylab = "log life expectancy")
points(log(LE2) ~ prec_cv, data = subset(all_data, Seed.bank %in% c('yes', 'unknown')), pch = 1, col =1)
points(log(LE2) ~ prec_cv, data = subset(all_data, Seed.bank == 'no'), pch = 2, col = 2)
legend("topright", c("With Seedbank", "No Seedbank"), col = c(1,2), pch = c(1,2))
cor1 = cor(log(all_data$LE2), all_data$prec_cv, use = 'complete.obs')
cortest2 = cor.test(log(all_data$LE2), all_data$prec_cv, use = 'complete.obs')
text1 = paste("r = ", round(cor1, 2), paste("; p =", round(as.numeric(cortest2[3]), 2)))
text(23,6.5, text1)
dev.copy2pdf(file = "fig_2.pdf")


boxplot(E_var_proportion ~ Seed.bank, all_data, ylab = "Relative Elasticity to Variance", xlab = "Seed Bank Status")
dev.copy2pdf(file = "fig_3.pdf")

plot(E_var_proportion ~ log(LE2), all_data, type = 'n', ylab = ylab1, xlab = "Log life expectancy")
points(E_var_proportion ~ log(LE2), data = subset(all_data, Seed.bank %in% c('yes', 'unknown')), pch = 1, col =1)
points(E_var_proportion ~ log(LE2), data = subset(all_data, Seed.bank == 'no'), pch = 2, col = 2)
legend("topright", c("With Seedbank", "No Seedbank"), col = c(1,2), pch = c(1,2))
cor2 = cor(all_data$E_var_proportion,  log(all_data$LE2), use = 'complete.obs')
cortest4 = cor.test(all_data$E_var_proportion,  log(all_data$LE2), use = 'complete.obs')
text2 = paste("r = ", round(cor2, 2), paste("; p =", round(as.numeric(cortest4[3]), 2)))
text(7,0.2, text2)
dev.copy2pdf(file = "fig_4.pdf")




plot(sum_E_var ~ log(R0), all_data)

boxplot(prec_cv ~ Seed.bank, all_data)

plot(log(LE2) ~ prec_cv, all_data)

plot(sum_E_var ~ log(GT), all_data)
plot(Lambda_s ~ prec_cv, all_data)
plot(log(GT) ~ prec_cv, all_data)
plot(log(R0) ~ prec_cv, all_data)
plot(Lambda_s ~ Lambda, all_data)

Lambda_red = all_data$Lambda - all_data$Lambda_s
plot(Lambda_red ~ prec_cv, all_data)
plot(Lambda_red ~ log(LE2), all_data)

plot(Lambda_red ~ log(R0), all_data)
plot(Lambda_red ~ log(GT), all_data)

plot(sum_SE_var ~ sum_FE_var, all_data)
plot(Lambda_red ~ sum_E_var, all_data)

levels(all_data$SPP)


