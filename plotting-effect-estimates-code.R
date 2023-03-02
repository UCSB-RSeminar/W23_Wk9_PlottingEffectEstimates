## Plotting estimates from model effects in R!
## created PAG 2 March 2023
#####
#setup stuff
#####
#load libraries
library(tidyverse)
library(lme4)

#set working directory--I use navigation bar to set to source file location

#load the dataset from PAG et al 2021
d <- as_tibble(read.csv("data-for-plotting.csv"))

#check out the dataset
d

#####
#building our statistical model
#####
#we want to know the proportion of the group that approaches the stimulus according to stimulus type, group size, and location of trial
#note we initially had 2-way interactions but I've removed these b/c they didn't matter.
reduced2mmodel <- glmer(data = prop2mdata, formula = cbind(as.integer(mean.in2m), as.integer(mean.out2m)) ~ 
                          stim + field.grp.sz + loc +
                          (1|grp), family = binomial(link = "logit"))
#singular fit is driven by zero RE estimation. Patrick can describe this later for those interested. 
#test for significant fixed effects
drop1(reduced2mmodel, test = "Chisq") #significant effects of stimulus type and group size!
summary(reduced2mmodel)

#post-hoc tests of specific other effects
#using glht and following directions in documentation for multiple comparisons (mcp)
#test what stimulus differences are driving the overall stimulus effect
summary(glht(reduced2mmodel, linfct = mcp(stim = c("t_scents - c_scents = 0",
                                                   "c_intruders - c_scents = 0",
                                                   "t_intruders - t_scents = 0",
                                                   "t_intruders - c_intruders = 0"))))
#looks like scents and intruders treatment effects are really important. 
#control responses don't differ and treatment responses don't differ

#summary stats of the stimulus effect
prop2mdata%>%
  group_by(stim)%>%
  summarise(mean(stim.mean.prop2m), sd(stim.mean.prop2m))

#####
#plotting raw data for categorical predictors
#you try it first???
#####
######
#Patrick's answer here
#####
#reorder stimulus first
prop2mdata$stim <- factor(prop2mdata$stim, c("c_scents", "t_scents", "c_intruders", "t_intruders"))
ggplot(data = prop2mdata, aes(y = stim.mean.prop2m, x = stim))+
  geom_boxplot()+
  scale_y_continuous(breaks = seq(0,1.2,0.25), limits = c(-0.05, 1.15))+
  xlab("Stimulus")+
  ylab("Proportion approaching stimulus")+
  scale_x_discrete(labels = c("c_scents" = "Control Scents", "t_scents" = "Rival Scents",
                              "c_intruders" = "Control Intruders", "t_intruders" = "Rival Intruders"))+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))
#can add more shit if we want or pretty it up
#but we're only looking at the raw data here! 
#we don't know how the estimates from the model output match up against the raw data
#this is important b/c the model accounts for group ID, other fixed effects, etc....

#####
#plotting model effects for categorical predictors
#####
#We plotted raw data as transparent points and the model means & SEs from the reduced model as darker, larger points w/ associated vertical bars.
#create newdata dataset to be filled in with predicted values from the reduced model
#this essentially creates an empty dataset that we will fill in
newdata <- with(prop2mdata,
                expand.grid(
                  stim = c("c_scents", "t_scents", "c_intruders", "t_intruders"), #add all 4 categories of scent
                  field.grp.sz = mean(field.grp.sz), #add the mean field group size
                  mean.in2m = mean(mean.in2m), #the mean prop w/in 2m
                  mean.out2m = mean(mean.out2m), #mean prop outside 2m
                  loc = c("Core", "Non-core"), #location options
                  stim.mean.prop2m=0 #this will be filled in from the model predictions
                ))
#compare to prop2mdata
prop2mdata

#take the values of response variable from the reduced model above
mm <- model.matrix(terms(reduced2mmodel), newdata)
newdata$stim.mean.prop2m <- predict(reduced2mmodel, newdata, re.form=NA, type="response")
newdata$stim.mean.prop2m.logit <- predict(reduced2mmodel, newdata, re.form=NA) #gives predicted values for the response variable from the model; used for confidence intervals
pvar1 <- diag(mm %*% tcrossprod(vcov(reduced2mmodel),mm)) #variance from the model, used to calculate confidence intervals

#confidence intervals--these have to be calculated by hand
newdata <- data.frame(newdata,
                      plo.ci = newdata$stim.mean.prop2m.logit-2*sqrt(pvar1), #confidence intervals (2*square root of variance)
                      phi.ci = newdata$stim.mean.prop2m.logit+2*sqrt(pvar1),
                      plo.se = newdata$stim.mean.prop2m.logit-sqrt(pvar1), #standard errors (square root of variance)
                      phi.se = newdata$stim.mean.prop2m.logit+sqrt(pvar1))

#convert to response scale
#this function converts logit to probability https://sebastiansauer.github.io/convert_logit2prob/ 
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#add upper and lower confidence intervals, and upper and lower standard errors, to newdata
newdata$uci <- logit2prob(newdata$phi.ci)
newdata$lci <- logit2prob(newdata$plo.ci)
newdata$use <- logit2prob(newdata$phi.se)
newdata$lse <- logit2prob(newdata$plo.se)

#check out newdata
newdata

#average the values across location treatments
#we decided to do this because location did not significant impact responses. 
plotnewdata <- newdata%>%
  group_by(stim)%>%
  summarise(mean(stim.mean.prop2m), mean(lse), mean(use))
names(plotnewdata) <- c("stim", "stim.mean.prop2m", "lse", "use")

#now make the plot!
ggplot(d = prop2mdata, aes(x = stim, y = stim.mean.prop2m))+ #start with the raw data
  geom_jitter(width = 0.05, alpha = 0.25, pch=ifelse(grepl("c_", prop2mdata$stim), 17,19), cex=2.5)+ #use geom_jitter to add raw data points
  geom_point(d = plotnewdata, aes(x = stim, y = stim.mean.prop2m), shape = rep(c(17,19),2), cex = 4)+ #now use geom_point to add the model estimate means
  geom_linerange(d = plotnewdata, aes(x=stim,ymin=lse,ymax=use),size=1.5,position=position_dodge(width=0.25)) + #and geom_linerange to add the standard errors
  #things below this line just make it look pretty
  scale_y_continuous(breaks = seq(0,1.2,0.25), limits = c(-0.05, 1.15))+
  xlab("Stimulus")+
  ylab("Proportion approaching stimulus")+
  scale_x_discrete(labels = c("c_scents" = "Control Scents", "t_scents" = "Rival Scents",
                              "c_intruders" = "Control Intruders", "t_intruders" = "Rival Intruders"))+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))

#####
#YOU CAN MESS W/ THIS PLOT, TOO!
#e.g., how would we plot confidence intervals instead of standard errors?
######
#####
#Patrick's answer here
#####
plotnewdata2 <- newdata%>%
  group_by(stim)%>%
  summarise(mean(stim.mean.prop2m), mean(lci), mean(uci))
names(plotnewdata2) <- c("stim", "stim.mean.prop2m", "lci", "uci")

ggplot(d = prop2mdata, aes(x = stim, y = stim.mean.prop2m))+ #start with the raw data
  geom_jitter(width = 0.05, alpha = 0.25, pch=ifelse(grepl("c_", prop2mdata$stim), 17,19), cex=2.5)+ #use geom_jitter to add raw data points
  geom_point(d = plotnewdata2, aes(x = stim, y = stim.mean.prop2m), shape = rep(c(17,19),2), cex = 4)+ #now use geom_point to add the model estimate means
  geom_linerange(d = plotnewdata2, aes(x=stim,ymin=lci,ymax=uci),size=1.5,position=position_dodge(width=0.25)) + #and geom_linerange to add the standard errors
  #things below this line just make it look pretty
  scale_y_continuous(breaks = seq(0,1.2,0.25), limits = c(-0.05, 1.15))+
  xlab("Stimulus")+
  ylab("Proportion approaching stimulus")+
  scale_x_discrete(labels = c("c_scents" = "Control Scents", "t_scents" = "Rival Scents",
                              "c_intruders" = "Control Intruders", "t_intruders" = "Rival Intruders"))+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))
#just gives you some wider confidence intervals, as expected. 

#####
#plotting raw data for continuous predictors w/ ggplot
#####
ggplot(data = prop2mdata, aes(x = field.grp.sz, y = stim.mean.prop2m))+
  geom_point()+
  geom_smooth(method = "lm", col = "black")+
  scale_y_continuous(breaks = seq(0,1.2,0.25), limits = c(-0.05, 1.15))+
  xlab("Field group size")+
  ylab("Proportion approaching stimulus")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))
#so this looks fine, but again it's not showing the actual estimate from the model!
#e.g., in the model, the error doesn't change across the group size range, so the ribbon won't widen at either end!

#####
#plotting model effects for continuous predictors
#####
#essentially the same process as before, but we've changed just a bit
newdata <- with(prop2mdata,
                expand.grid(
                  stim = c("c_scents", "t_scents", "c_intruders", "t_intruders"),
                  field.grp.sz = seq(min(field.grp.sz), max(field.grp.sz), 1), #field group size is now a sequence of continuous data!
                  mean.in2m = mean(mean.in2m),
                  mean.out2m = mean(mean.out2m),
                  loc = c("Core", "Non-core"),
                  stim.mean.prop2m=0
                ))
newdata #compared to the example above, newdata is a loooooong dataset

#values of response variable
#we now calculate these for EVERY value of field group size!
mm <- model.matrix(terms(reduced2mmodel), newdata)
newdata$stim.mean.prop2m <- predict(reduced2mmodel, newdata, re.form=NA, type="response")
newdata$stim.mean.prop2m.logit <- predict(reduced2mmodel, newdata, re.form=NA)
pvar1 <- diag(mm %*% tcrossprod(vcov(reduced2mmodel),mm))

#confidence intervals
newdata <- data.frame(newdata,
                      plo.ci = newdata$stim.mean.prop2m.logit-2*sqrt(pvar1), #confidence intervals (2*square root of variance)
                      phi.ci = newdata$stim.mean.prop2m.logit+2*sqrt(pvar1),
                      plo.se = newdata$stim.mean.prop2m.logit-sqrt(pvar1), #standard errors (square root of variance)
                      phi.se = newdata$stim.mean.prop2m.logit+sqrt(pvar1))

#convert to response scale
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

newdata$uci <- logit2prob(newdata$phi.ci)
newdata$lci <- logit2prob(newdata$plo.ci)
newdata$use <- logit2prob(newdata$phi.se)
newdata$lse <- logit2prob(newdata$plo.se)

plotnewdata <- newdata%>%
  group_by(field.grp.sz)%>% #grouping by field group size because this is what we are plotting
  summarise(mean(stim.mean.prop2m), mean(lse), mean(use))
names(plotnewdata) <- c("field.grp.sz", "stim.mean.prop2m", "lse", "use")

#now plot in ggplot
ggplot(data = prop2mdata, aes(x = field.grp.sz, y = stim.mean.prop2m))+
  geom_point()+
  geom_line(data = plotnewdata, aes(x = field.grp.sz, y = stim.mean.prop2m), lwd = 1)+
  geom_ribbon(data = plotnewdata, aes(ymin = lse, ymax = use), alpha = 0.25)+
  ylim(c(-0.05, 1.05))+
  xlab("Group size")+
  ylab("Proportion approaching stimulus")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text = element_text(size = 12))


#####
# another technique for other types of models
# original source: https://quantdev.ssri.psu.edu/tutorials/five-ish-steps-create-pretty-interaction-plots-multi-level-model-r
#####
#####
#load data and package
#####
library(effects)
#here we have a dataset on how many intergroup fights ("igis") male mongooses have been in over their lifetime
igis.age <- as_tibble(read.csv("igis-age-data.csv"))
#are older males involved in more fights?

#check data distribution before building model
hist(igis.age$n.igis.tot)
summary(igis.age$n.igis.tot) #really right-skewed!
#not zero-inflated but essentially "one-inflated". Lots of values = 1. 
#we originally did a zero-inflated model in the paper but I've simplified to remove that here. 

#create an observation-level random effect to deal w/ overdispersion
igis.age$olre <- as.factor(row.names(igis.age))

#build the model--a negative binomial model!
#using scaled predictor and obs-level random effect
mod.olre <- glmmTMB(data = igis.age, formula = n.igis.tot~scale(max.age.tot) + (1|olre), family = nbinom2)
hist(residuals(mod.olre))#histogram of residuals looks good
drop1(mod.olre, test = "Chisq") #use drop1 to test significance
summary(mod.olre)

#####
#plot with ggplot?
#####
#fuck if i know how to do that w/ negative binomial model
ggplot(data = igis.age, aes(x = max.age.tot/365, y = n.igis.tot))+ #i've adjusted max.age.tot to be in years, not days
  geom_point()+
  geom_smooth(method = "glm", method.args = list(family = nbinom2), col = "black")+
  ylim(c(0, 200))+
  scale_x_continuous(breaks = c(0,3,6,9,12))+
  xlab("Age (years)")+
  ylab("Number of lifetime intergroup contests")
#I guess it can be done!
#gets a bit weird and has same problems as prior plots

#####
#plot model effect for continuous predictors, v2
#####
#use the effect package
#you basically tell effects the following
#the term/variable you want to plot on the x axis using the language in the model (here, max.age.tot, but w/ scale())
#the "levels" of that variable (here, a sequence from minimum to max, by 100 days)
#the model you want to plot from
#MAKE SURE you recognize it as a data frame!
ef.data <- as.data.frame(effect(term = "scale(max.age.tot)", xlevels = list(max.age.tot = seq(min(igis.age$max.age.tot), max(igis.age$max.age.tot), by = 100)), mod = mod.olre))
#you also have to make a new "fit" variable in your RAW data that is the same as your y-axis variable.
#this is just easier for ggplot to figure out. 
igis.age$fit <- igis.age$n.igis.tot

#now plot it!
ggplot(data = ef.data, aes(x = max.age.tot/365, y = fit))+
  geom_line()+
  geom_ribbon(aes(ymin = fit-se, ymax = fit+se), alpha = 0.25)+
  geom_point(data = igis.age, aes(x = max.age.tot/365, y = fit))+ #note the change in dataset here!
  scale_x_continuous(breaks = c(0,3,6,9,12))+
  xlab("Age (years)")+
  ylab("Number of lifetime intergroup contests")

#ggsave("NIGIs~age.png", dpi=300, height=4, width=6, units="in")

#####
#same package works for simple linear models, too
#####
#even though this is the wrong model to fit, let's just try it
mod.lm <- lm(data = igis.age, formula = n.igis.tot~scale(max.age.tot))
summary(mod.lm)

ef.data.linear <- as.data.frame(effect(term = "scale(max.age.tot)", xlevels = list(max.age.tot = seq(min(igis.age$max.age.tot), max(igis.age$max.age.tot), by = 100)), mod = mod.lm))
#you also have to make a new "fit" variable in your RAW data that is the same as your y-axis variable.
#this is just easier for ggplot to figure out. 
igis.age$fit <- igis.age$n.igis.tot

#now plot it!
ggplot(data = ef.data.linear, aes(x = max.age.tot/365, y = fit))+
  geom_line()+
  geom_ribbon(aes(ymin = fit-se, ymax = fit+se), alpha = 0.25)+
  geom_point(data = igis.age, aes(x = max.age.tot/365, y = fit))+
  scale_x_continuous(breaks = c(0,3,6,9,12))+
  xlab("Age (years)")+
  ylab("Number of lifetime intergroup contests")
