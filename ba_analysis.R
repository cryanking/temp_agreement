## docker run --rm -it -v "~/gitdirs/misc_applied_in_dev/sa_temp/:/research" cryanking/temperature_container R

library(lme4)
library(dplyr)
library(splines)
library(pbkrtest)
library(readr)
library(tidyr)
library(magrittr)
setwd("/research")

set.seed(1002)
normal_sample <- exp(rnorm(2000))

find_loa_factor <- function(sd=1, width=0.95) {
if(! exists("normal_sample")) { normal_sample<- exp(rnorm(2000)) }
uniroot(f= function(x){ mean(pchisq(x*x *normal_sample^sd , df=1)) - 0.95}, lower= 0.01, upper=20 )$root
}

loa_from_fit <- function(local_fit) {

## sd est
  total_sd <- (local_fit %>% summary %>% extract2("varcor") %>% extract2("Case_Number") %>% as.numeric + local_fit %>% summary %>% extract2("sigma") %>% raise_to_power(2) ) %>% sqrt

  
local_boot <- bootMer(local_fit, FUN=function(local_fit){ log(local_fit%>% summary %>% extract2("varcor") %>% extract2("Case_Number") %>% as.numeric + local_fit %>% summary %>% extract2("sigma") %>% raise_to_power(2) ) } , nsim=100L)%>% extract2("t") %>% as.numeric 

boot_sd <- local_boot %>% sd
boot_sd_limits <- local_boot %>% quantile(c(0.025, .975)) %>% exp %>% sqrt

local_loa_mult <- boot_sd %>% find_loa_factor(sd=.)


## extract var(mean) and mean(mean) of fe
prediction_variance <- matrix(colMeans(model.matrix(local_fit, type="fixed")), nrow=1) %*% tcrossprod(pbkrtest::vcovAdj(local_fit) , matrix( colMeans( model.matrix(local_fit, type="fixed")), nrow=1 ) ) %>% as.numeric


## bias_avg, bias_sd
local_output <- c(  mean(fitted(local_fit)), sqrt(prediction_variance ) , total_sd*local_loa_mult ,  boot_sd_limits*local_loa_mult  )


names(local_output ) <- c("bias_avg", "bias_sd", "loa_est", "loa_lower", "loa_upper")

local_output  
}


## main dataset
main_data <- read_csv("ob_temp_working_data.csv")  
main_data %<>% filter(Incl_Excl == "incl")

## pivot to a "long" form; this is a function so that per-patient bootstrapping is easy, but that isn't actually needed in this analysis
functional_longify <- function(data, indicies=NULL) {
  if(! is.null(indicies)) {
    data <- data[indicies,]
    data$Case_Number <- seq.int(nrow(data))
  }
  
  ## require 3/5 the first points to be present 
  data %<>% filter(  is.na(t_drager_40)+is.na(t_drager_30)+is.na(t_drager_20)+is.na(t_drager_10)+is.na(t_drager_00)  < 3 )
  
  data %<>% select( c(one_of("Case_Number"), contains("t_drager"),  contains("t_ir"), contains("t_oral") ) )
  data %<>% select( -t_drager_arr,-t_drager_lve )
 
  ## filter out changes of > 1 degree in the first 10 minutes 
  data %<>% mutate(t_drager_00 = if_else( abs(t_drager_00-t_drager_10) > 1, t_drager_10, t_drager_00  ) )

  data %<>% pivot_longer(cols=contains("t_"), names_to=c("source","timepoint"), names_prefix="t_", names_sep="_", values_to="temper")
  data$timepoint %<>% as.numeric

  data$Case_Number %<>% factor
  return(data)
}

## compute pairwise differences in a slightly wider format
## note that using deltas means that some potential smoothing information is lost on the more frequent signal

main_data %>% functional_longify  %>% pivot_wider(id_cols = one_of("Case_Number", "timepoint" ), names_from="source", values_from="temper"  ) %>% mutate(delta1=drager-ir, delta2=drager-oral, delta3=oral-ir) -> delta_data

##########
##drager vs ir

## with only 4 timepoints, just make time cat
## random intercept -> per person variation in calibration
  drager_ir_fit <- lmer( delta1~factor(timepoint) + (1|Case_Number) , data=delta_data)

  drager_ir_out <- loa_from_fit(drager_ir_fit)
  drager_ir_out %>% round(2)

  delta_data$delta1 %>% quantile(probs=c(.025, .05,.1586, .5, .8413 , .95, .975), na.rm=T) %>% round(2)


##########
##drager vs oral
  drager_oral_fit <- lmer( delta2~factor(timepoint) + (1|Case_Number) , data=delta_data)
  
  drager_oral_out <- loa_from_fit(drager_oral_fit)
  drager_oral_out %>% round(2)


  delta_data$delta2 %>% quantile(probs=c(.025, .05,.1586, .5, .8413 , .95, .975), na.rm=T) %>% round(2)


##########
##oral vs ir
  ir_oral_fit <- lmer( delta3~factor(timepoint) + (1|Case_Number) , data=delta_data)
  ir_oral_out <- loa_from_fit(ir_oral_fit)
  ir_oral_out %>% round(2)


delta_data$delta3 %>% quantile(probs=c(.025, .05,.1586, .5, .8413 , .95, .975), na.rm=T) %>% round(2)


overal_output <- bind_rows(drager_ir_out, drager_oral_out, ir_oral_out)

 overal_output %>% write_csv("agreement_output.csv")

## pictures
mycol <- rgb(0, 0, 0 , max = 255, alpha = 100, names = "black50")
mycol2 <- rgb(255, 0, 0 , max = 255, alpha = 100, names = "red50")
mycol3 <- rgb(0, 0, 255 , max = 255, alpha = 50, names = "blue25")


jpeg("drager_ir_ba.jpg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
delta_data %>% filter(ir > 34) %>% filter(drager > 34) %>% mutate(x=(drager+ir)/2, y=delta1 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Drager versus IR", type="p", xlim=c(34,39), ylim=c(-4,4))}
abline(h=drager_ir_out[["bias_avg"]] , col="black")
abline(h=drager_ir_out[["bias_avg"]] + drager_ir_out[["bias_sd"]]*c(-2,2), col=mycol)

abline(h=drager_ir_out[["bias_avg"]] + drager_ir_out[["loa_est"]]*c(-1,1) , col="red")
abline(h=drager_ir_out[["bias_avg"]] + drager_ir_out[["loa_lower"]]*c(-1,1) , col=mycol2)
abline(h=drager_ir_out[["bias_avg"]] + drager_ir_out[["loa_upper"]]*c(-1,1) , col=mycol2)
polygon(x=matrix(c(25, -.5, 50,-.5, 50, .5, 25,.5  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

delta_data %>% filter(ir > 34) %>% filter(drager > 34) %>% mutate(x=ir, y=drager) %>% {plot(x=.$x, y=.$y, xlab="IR", ylab="Drager", main="Drager versus IR", type="p")}

## repeated measures analysis
rm_glm <- delta_data %>% filter(drager > 34) %>% filter(ir > 34)%>% lmer(drager~ir+(1|Case_Number), data=.)
abline(a= rm_glm %>% summary %>% extract2("coefficients") %>% extract(1,1) , b=rm_glm %>% summary %>% extract2("coefficients") %>% extract(2,1) , col="red")
polygon(x=matrix(c(25,24, 50,49, 50, 51,25,26  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

# delta_data %>% filter(drager > 34) %>% filter(ir > 34)%>% lm(drager~ir, data=.) %>% abline(col="red")
abline(0,1)
dev.off()


 jpeg("drager_oral_ba.jpg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
delta_data  %>% filter(oral > 34) %>% filter(drager > 34) %>% mutate(x=(drager+oral)/2, y=delta2 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Drager versus oral", type="p", xlim=c(34,39), ylim=c(-4,4))}
abline(h=drager_oral_out[["bias_avg"]] , col="black")
abline(h=drager_oral_out[["bias_avg"]] + drager_oral_out[["bias_sd"]]*c(-2,2), col=mycol)

abline(h=drager_oral_out[["bias_avg"]] + drager_oral_out[["loa_est"]]*c(-1,1) , col="red")
abline(h=drager_oral_out[["bias_avg"]] + drager_oral_out[["loa_lower"]]*c(-1,1) , col=mycol2)
abline(h=drager_oral_out[["bias_avg"]] + drager_oral_out[["loa_upper"]]*c(-1,1) , col=mycol2)
polygon(x=matrix(c(25, -.5, 50,-.5, 50, .5, 25,.5  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

delta_data %>% filter(oral > 34) %>% filter(drager > 34) %>% mutate(x=oral, y=drager) %>% {plot(x=.$x, y=.$y, xlab="Oral", ylab="Drager", main="Drager versus Oral", type="p")}

rm_glm <- delta_data %>% filter(drager > 34) %>% filter(oral > 34)%>% lmer(drager~oral+(1|Case_Number), data=.)

abline(a= rm_glm %>% summary %>% extract2("coefficients") %>% extract(1,1) , b=rm_glm %>% summary %>% extract2("coefficients") %>% extract(2,1)  , col="red")
polygon(x=matrix(c(25,24, 50,49, 50, 51,25,26  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

# delta_data %>% filter(oral > 34) %>% filter(drager > 34) %>% lm(drager~oral, data=.) %>%  abline(col="red")
abline(0,1)

dev.off()


jpeg("oral_ir_ba.jpg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
delta_data %>% filter(oral > 34) %>% filter(ir > 34) %>% mutate(x=(oral+ir)/2, y=delta3 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Oral versus IR", type="p", xlim=c(34,39), ylim=c(-4,4))}
abline(h=ir_oral_out[["bias_avg"]] , col="black")
abline(h=ir_oral_out[["bias_avg"]] + ir_oral_out[["bias_sd"]]*c(-2,2), col=mycol)

abline(h=ir_oral_out[["bias_avg"]] + ir_oral_out[["loa_est"]]*c(-1,1) , col="red")
abline(h=ir_oral_out[["bias_avg"]] + ir_oral_out[["loa_lower"]]*c(-1,1) , col=mycol2)
abline(h=ir_oral_out[["bias_avg"]] + ir_oral_out[["loa_upper"]]*c(-1,1) , col=mycol2)
polygon(x=matrix(c(25, -.5, 50,-.5, 50, .5, 25,.5  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

delta_data %>% filter(oral > 34) %>% filter(ir > 34) %>% mutate(x=ir, y=oral) %>% {plot(x=.$x, y=.$y, xlab="IR", ylab="Oral", main="Oral versus IR", type="p")}

rm_glm <- delta_data %>% filter(ir > 34) %>% filter(oral > 34)%>% lmer(oral~ir+(1|Case_Number), data=.)

abline(a= rm_glm %>% summary %>% extract2("coefficients") %>% extract(1,1) , b=rm_glm %>% summary %>% extract2("coefficients") %>% extract(2,1) , col="red")
polygon(x=matrix(c(25,24, 50,49, 50, 51,25,26  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

# delta_data %>% filter(oral > 34) %>% filter(ir > 34) %>% lm(oral~ir, data=.) %>%  abline(col="red")
abline(0,1)

dev.off()




delta_data %>% select(drager, oral, ir) %>% mutate_all( function(x){if_else(x<34, NA_real_, x) } ) %>% as.matrix %>% cor(use="pairwise") %>% round(2) %>% write.csv("pairwise_correlations.csv")

# delta_data %>% select(drager, oral, ir) %>% mutate_all( function(x){if_else(x<34, NA_real_, x) } ) %>% as.matrix %>% cov(use="pairwise") %>% round(2)


###############
## same plots without the 34 degree filter

jpeg("oral_ir_ba_no_filter.jpg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
delta_data  %>% mutate(x=(oral+ir)/2, y=delta3 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Oral versus IR", type="p" )}
abline(h=ir_oral_out[["bias_avg"]] , col="black")
abline(h=ir_oral_out[["bias_avg"]] + ir_oral_out[["bias_sd"]]*c(-2,2), col=mycol)

abline(h=ir_oral_out[["bias_avg"]] + ir_oral_out[["loa_est"]]*c(-1,1) , col="red")
abline(h=ir_oral_out[["bias_avg"]] + ir_oral_out[["loa_lower"]]*c(-1,1) , col=mycol2)
abline(h=ir_oral_out[["bias_avg"]] + ir_oral_out[["loa_upper"]]*c(-1,1) , col=mycol2)
polygon(x=matrix(c(25, -.5, 50,-.5, 50, .5, 25,.5  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

delta_data  %>% mutate(x=ir, y=oral) %>% {plot(x=.$x, y=.$y, xlab="IR", ylab="Oral", main="Oral versus IR", type="p")}

rm_glm <- delta_data %>% lmer(oral~ir+(1|Case_Number), data=.)

abline(a= rm_glm %>% summary %>% extract2("coefficients") %>% extract(1,1) , b=rm_glm %>% summary %>% extract2("coefficients") %>% extract(2,1) , col="red")
polygon(x=matrix(c(25,24, 50,49, 50, 51,25,26  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

# delta_data %>% filter(oral > 34) %>% filter(ir > 34) %>% lm(oral~ir, data=.) %>%  abline(col="red")
abline(0,1)

dev.off()






 jpeg("drager_oral_ba_no_filter.jpg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
delta_data  %>% mutate(x=(drager+oral)/2, y=delta2 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Drager versus oral", type="p")}
abline(h=drager_oral_out[["bias_avg"]] , col="black")
abline(h=drager_oral_out[["bias_avg"]] + drager_oral_out[["bias_sd"]]*c(-2,2), col=mycol)

abline(h=drager_oral_out[["bias_avg"]] + drager_oral_out[["loa_est"]]*c(-1,1) , col="red")
abline(h=drager_oral_out[["bias_avg"]] + drager_oral_out[["loa_lower"]]*c(-1,1) , col=mycol2)
abline(h=drager_oral_out[["bias_avg"]] + drager_oral_out[["loa_upper"]]*c(-1,1) , col=mycol2)
polygon(x=matrix(c(25, -.5, 50,-.5, 50, .5, 25,.5  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

delta_data  %>% mutate(x=oral, y=drager) %>% {plot(x=.$x, y=.$y, xlab="Oral", ylab="Drager", main="Drager versus Oral", type="p")}

rm_glm <- delta_data %>% lmer(drager~oral+(1|Case_Number), data=.)

abline(a= rm_glm %>% summary %>% extract2("coefficients") %>% extract(1,1) , b=rm_glm %>% summary %>% extract2("coefficients") %>% extract(2,1) , col="red")
polygon(x=matrix(c(25,24, 50,49, 50, 51,25,26  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

# delta_data %>% filter(oral > 34) %>% filter(drager > 34) %>% lm(drager~oral, data=.) %>%  abline(col="red")
abline(0,1)

dev.off()



jpeg("drager_ir_ba_no_filter.jpg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
delta_data  %>% mutate(x=(drager+ir)/2, y=delta1 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Drager versus IR", type="p")}
abline(h=drager_ir_out[["bias_avg"]] , col="black")
abline(h=drager_ir_out[["bias_avg"]] + drager_ir_out[["bias_sd"]]*c(-2,2), col=mycol)

abline(h=drager_ir_out[["bias_avg"]] + drager_ir_out[["loa_est"]]*c(-1,1) , col="red")
abline(h=drager_ir_out[["bias_avg"]] + drager_ir_out[["loa_lower"]]*c(-1,1) , col=mycol2)
abline(h=drager_ir_out[["bias_avg"]] + drager_ir_out[["loa_upper"]]*c(-1,1) , col=mycol2)
polygon(x=matrix(c(25, -.5, 50,-.5, 50, .5, 25,.5  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)
delta_data  %>% mutate(x=ir, y=drager) %>% {plot(x=.$x, y=.$y, xlab="IR", ylab="Drager", main="Drager versus IR", type="p")}

## repeated measures analysis
rm_glm <- delta_data %>% lmer(drager~ir+(1|Case_Number), data=.)
abline(a= rm_glm %>% summary %>% extract2("coefficients") %>% extract(1,1) , b=rm_glm %>% summary %>% extract2("coefficients") %>% extract(2,1) , col="red")
polygon(x=matrix(c(25,24, 50,49, 50, 51,25,26  ) , ncol=2, byrow=TRUE) , col=mycol3, lty=0)

# delta_data %>% filter(drager > 34) %>% filter(ir > 34)%>% lm(drager~ir, data=.) %>% abline(col="red")
abline(0,1)
dev.off()



###################
## reshape to calculate the change over the 30 minutes and delta of that
delta_data %>% group_by(Case_Number) %>% summarize( across(starts_with("delta"), function(x) {last(na.omit(x) ) - first(na.omit(x))} )) -> ddeltas


main_data %>% mutate(t_drager_00 = if_else( abs(t_drager_00-t_drager_10) > 1, t_drager_10, t_drager_00  ) ) %>% 
mutate( drgaer_delta = t_drager_30 - t_drager_00 , ir_delta = t_ir_30 - t_ir_00, oral_delta = t_oral_30 - t_oral_00) %>% select(Case_Number, drgaer_delta, ir_delta, oral_delta) %>% mutate(Case_Number = factor(Case_Number)) %>% inner_join(ddeltas, by="Case_Number") -> ddeltas


## this is a standard BA analysis
semean <- function(x){sd(x, na.rm=TRUE)/sqrt(sum(is.finite(x)))}

jpeg("changes_drager_ir_ba_no_filter.jpeg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
ddeltas %>% mutate(x=(drgaer_delta+ir_delta)/2, y=delta1 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Drager versus IR", type="p")}

abline(h=ddeltas[["delta1"]] %>% mean , col="black")
abline(h=ddeltas[["delta1"]] %>% mean + ddeltas[["delta1"]] %>% semean %>% multiply_by(c(-2,2)), col=mycol)

abline(h=ddeltas[["delta1"]] %>% mean + ddeltas[["delta1"]] %>%sd %>% multiply_by(c(-2,2)) , col="red")

ddeltas  %>% mutate(x=ir_delta, y=drgaer_delta) %>% {plot(x=.$x, y=.$y, xlab="IR", ylab="Drager", main="Drager versus IR", type="p")}
abline(0,1)
ddeltas %>% lm(drgaer_delta~ir_delta, data=.) %>%  abline(col="red")

dev.off()

jpeg("changes_frager_oral_ba_no_filter.jpeg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
ddeltas %>% mutate(x=(drgaer_delta+oral_delta)/2, y=delta2 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Drager versus oral", type="p")}

abline(h=ddeltas[["delta2"]] %>% mean , col="black")
abline(h=ddeltas[["delta2"]] %>% mean + ddeltas[["delta2"]] %>% semean %>% multiply_by(c(-2,2)), col=mycol)

abline(h=ddeltas[["delta2"]] %>% mean + ddeltas[["delta2"]] %>%sd %>% multiply_by(c(-2,2)) , col="red")

ddeltas  %>% mutate(x=oral_delta, y=drgaer_delta) %>% {plot(x=.$x, y=.$y, xlab="oral", ylab="Drager", main="Drager versus oral", type="p")}
abline(0,1)
ddeltas %>% lm(drgaer_delta~oral_delta, data=.) %>%  abline(col="red")

dev.off()

jpeg("changes_oral_ir_ba_no_filter.jpeg", res=300, width=12, height=6, units="in")
par(mfrow=c(1,2))
ddeltas %>% mutate(x=(oral_delta+ir_delta)/2, y=delta3 ) %>% {plot(x=.$x, y=.$y, xlab="average", ylab="difference", main="Oral versus IR", type="p")}

abline(h=ddeltas[["delta3"]] %>% mean , col="black")
abline(h=ddeltas[["delta3"]] %>% mean + ddeltas[["delta3"]] %>% semean %>% multiply_by(c(-2,2)), col=mycol)

abline(h=ddeltas[["delta3"]] %>% mean + ddeltas[["delta3"]] %>%sd %>% multiply_by(c(-2,2)) , col="red")

ddeltas  %>% mutate(x=oral_delta, y=ir_delta) %>% {plot(x=.$x, y=.$y, xlab="Oral", ylab="IR", main="Oral versus IR", type="p")}
abline(0,1)
ddeltas %>% lm(ir_delta~oral_delta, data=.) %>%  abline(col="red")

dev.off()

ddeltas %>% select(drgaer_delta, oral_delta, ir_delta)  %>% as.matrix %>% cor(use="pairwise") %>% round(2) %>% write.csv("pairwise_correlations_changes.csv")


library("tableone")

cont_vars <- c(
"Age"                  ,
"Weight"               ,
"Height"    ,
"Gest_Age_compl_weeks" ,
"Gest_age_actual"      ,
"Gravidity"            ,
"Parity" ,
"Ambient_T_theatre"    ,
"T_Recovery_ambient"   ,
"Ephedrine"            ,
"Phenyl_bolus"  ,
"Blood loss"           ,
"APGAR_5min" 
)

main_data %<>% mutate(across(one_of(cont_vars), as.numeric ))
main_data[["Co-morb"]] %<>% tolower

cat_vars <- c(
"Singleton"            ,
"Co-morb" ,
"C/S_Urg_Cat"          ,
"Active_Labour"        ,
"Standard_SA_dose"     ,
"Warming_FAW"          ,
"Warming_Warm_Fluid"  ,
"Oxy_Bolus"            ,
"Oxy_infusion_intraop" ,
"Ocy_infusion_postop"  ,
"Phenyl_infu"          ,
"Spinal_Level"         ,
"Convert_GA"           ,
"Prophyl_Phenyl_inf"   ,
"Shivering"            ,
"Vomiting"             ,
"Neonate_admission"    ,
"Neonatal_CPR"  
)


tab3 <- CreateTableOne(vars = c(cont_vars, cat_vars)  , data = main_data   , factorVars = cat_vars,  includeNA = FALSE, argsExact=list(simulate.p.values=TRUE))

tab3 %>% print(contDigits=0, printToggle=FALSE, nonnormal=TRUE, smd=TRUE) %>% (kableExtra::kbl) %>% kableExtra::save_kable(file="table1_ba.html")
tab3 %>% print(contDigits=0, printToggle=FALSE, nonnormal=TRUE, smd=TRUE) %>% write.csv("table1_ba.csv")
