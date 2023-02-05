

# install.packages("lmeresampler") ## this package would make the bootstrap faster, but doesn't seem to work with other packages right now

setwd("/research")
library(dplyr)
library(tidyr)
library(readr)
library(lme4)
library(splines)
library(boot)
library(parallel)
library(purrr)
library(magrittr)

run_large <- TRUE

set.seed(101)

## this happens to get all the column types that I care about correct
main_data <- read_csv("ob_temp_working_data.csv")  
main_data %<>% filter(Incl_Excl == "incl")
main_data %<>% filter( Convert_GA == 0)
main_data %<>% filter( Age >= 18)

starting_temperature <- 36.8
hypothermia_threshold <- 36.0
hypothermia_threshold2 <- 35.0 
max_time_in_prediction <- 90
use_pacu <- FALSE


max_loss <-  function(x) { min(x-x[1] , na.rm=TRUE ) }

## pivot to a "long" form; this is a function so that per-patient bootstrapping is easy
functional_longify <- function(data, indicies=NULL) {
  if(! is.null(indicies)) {
    data <- data[indicies,]
    data$Case_Number <- seq.int(nrow(data))
  }
  
  ## require the first 3/5 points to be present - could loose this later
  data %<>% filter(  is.na(t_drager_40)+is.na(t_drager_30)+is.na(t_drager_20)+is.na(t_drager_10)+is.na(t_drager_00)  < 3 )

  if(use_pacu){
    pacu_data <- data %>% select( one_of("Case_Number", "t_drager_arr", "t_drager_lve", "t_SA", "t_rec_arrival", "t_rec_leave") ) 
  } 
  
  data %<>% select( c(one_of("Case_Number"), contains("t_drager_") ) )
  data %<>% select( -t_drager_arr,-t_drager_lve )

  
  ## filter out changes of > 1 degree in the first 10 minutes 
  data %<>% mutate(t_drager_00 = if_else( abs(t_drager_00-t_drager_10) > 1, t_drager_10, t_drager_00  ) )

  data %<>% pivot_longer(cols=contains("t_drager_"), names_to="timepoint", names_prefix="t_drager_",values_to="temper")
  data$timepoint %<>% as.numeric

  ## join in the pacu observations
  if(use_pacu){
    pacu_data %<>% mutate(t_rec_leave = difftime(t_rec_leave , t_SA, units="mins") %>% as.numeric,  t_rec_arrival=difftime(t_rec_arrival , t_SA, units="mins") %>% as.numeric)

    data <- bind_rows(data,  
      pacu_data %>% select(Case_Number, temper=t_drager_arr, timepoint=t_rec_arrival) , 
      pacu_data %>% select(Case_Number, temper=t_drager_lve, timepoint=t_rec_leave)  
    )
  }
  
  data$Case_Number %<>% factor
  return(data)
}


## spline functions with appropriate windows
if(use_pacu) {
  bs_fun <- . %>% bs(Boundary.knots = c(0, 250), knots = c(30, 60, 90, 120, 180), degree=1)
  } else {
  bs_fun <- . %>% bs(Boundary.knots = c(0, 100), knots = c(30, 60), degree=1)
}


##################
## random slope version
##################

my.boot.fun <- function(data, indicies=NULL) {
  data <- functional_longify(data, indicies)

  ## get the number of basis elements
  modified_data <-  bs_fun(data$timepoint)
  hold_names <- paste0("ts" , seq.int(ncol(modified_data )) ) 
  modified_data %<>% as.data.frame %>% set_colnames(hold_names)
  data <- bind_cols(data, modified_data  )
  
  ts_form <- formula(paste("temper ~ (1|Case_Number) + " , paste0( hold_names , collapse="+" ) , "+" ,paste0("(" , hold_names, "-1|Case_Number)" , collapse="+" )  ) )

  ri_no_rs <- lmer( ts_form , data=data)
  resid_distribution <- residuals(ri_no_rs )
  
  
  fake_data <- data.frame(Case_Number=data$Case_Number %>% unique %>% rep(each=max_time_in_prediction+1 ), timepoint=seq(from=0, to=max_time_in_prediction, by=1) )
  
  modified_data <-  bs_fun(fake_data$timepoint) %>% as.data.frame %>% set_colnames(hold_names)
  modified_data <- bind_cols(fake_data, modified_data  )

  
  modified_data$predictions<- ri_no_rs %>% predict(newdata=modified_data  ) 
  
  modified_data %<>% arrange(Case_Number, timepoint)
  modified_data %<>% group_by(Case_Number) %>% mutate(coldnow = predictions - first(predictions) + starting_temperature) %>% ungroup
## main ananlyis: smoothed temperatures
## 
  ## average change
  ## probability < 36.0 ever
  ## probability change < -1 ever
  ## probability < 36.0 60
  ## probability change < -1 ever
  average_loss <- modified_data %>% group_by(Case_Number) %>% 
    summarize(
      delta0=max_loss(predictions), delta_600 = max_loss(predictions[timepoint <= 60] ) 
      , time_to_min_a = timepoint[which.min(predictions)]
      ) %>% 
    summarize(
      time_to_min = mean(time_to_min_a) ,
      delta_ever=mean(delta0), 
      delta_60=mean(delta_600), 
      less_1 = mean(delta0 < -1) , 
      less_1_60=mean(delta_600 < -1) , 
      hypothermia_ever= delta0 %>% add(starting_temperature) %>% is_less_than(hypothermia_threshold) %>% mean, 
      hypothermia_ever_low= delta0 %>% add(starting_temperature) %>% is_less_than(hypothermia_threshold2) %>% mean, 
      hypothermia_60= delta_600 %>% add(starting_temperature) %>% is_less_than(hypothermia_threshold) %>% mean,
      hypothermia_60_low= delta_600 %>% add(starting_temperature) %>% is_less_than(hypothermia_threshold2) %>% mean,
      sd_time_to_min = sd(time_to_min_a) ,
      sd_delta0 = sd(delta0), 
      sd_delta_60 = sd(delta_600)
      ) %>% 
      unlist
  
  average_loss <- c( average_loss, measured_very_low=   modified_data %>% group_by(Case_Number) %>% summarize(very_low = any(predictions < hypothermia_threshold2, na.rm=TRUE) ) %>% summarize(mean(very_low)) %>% unlist %>% as.numeric %>% unname )

  ## sensitivity analysis: only actually measured timepoints
  exclude_unmeasured <- modified_data %>% semi_join(data %>% filter(!is.na(temper)), by=c("Case_Number", "timepoint")) %>% group_by(Case_Number) %>% 
    summarize(
      delta0=max_loss(predictions), delta_600 = max_loss(predictions[timepoint <= 60] ) 
      , time_to_min = timepoint[which.min(predictions)]
      ) %>% 
    summarize(
      time_to_min = mean(time_to_min) ,
      delta_ever=mean(delta0), 
      delta_60=mean(delta_600), 
      less_1 = mean(delta0 < -1) , 
      less_1_60=mean(delta_600 < -1) , 
      hypothermia_ever= delta0 %>% add(starting_temperature) %>% is_less_than(hypothermia_threshold) %>% mean, 
      hypothermia_ever_low= delta0 %>% add(starting_temperature) %>% is_less_than(hypothermia_threshold2) %>% mean, 
      hypothermia_60= delta_600 %>% add(starting_temperature) %>% is_less_than(hypothermia_threshold) %>% mean,
      hypothermia_60_low= delta_600 %>% add(starting_temperature) %>% is_less_than(hypothermia_threshold2) %>% mean
) %>% 
      unlist 
  names( exclude_unmeasured) %<>% paste0("_ex")
      
  average_loss <- c(average_loss,exclude_unmeasured)
  
  ## curve on smoothed mean for CI
  mean_curve <- modified_data %>%filter((timepoint %% 10  == 0)) %>% group_by(timepoint) %>% summarize(smooth_temp=mean(predictions ) ) %>% arrange(timepoint)
  
  
  if(is.null(indicies)) {
    return(list(c(average_loss ), modified_data,ri_no_rs,mean_curve ) )
  } else {
  ## get bootstrap variability of average overall_trend
    fake_data <- data.frame(Case_Number=modified_data$Case_Number[1], timepoint=seq(from=0, to=max_time_in_prediction, by=1) )
    fake_data %<>% bind_cols(fake_data$timepoint %>% bs_fun %>%  as.data.frame %>% set_colnames( colnames(modified_data) %>% keep( ~grepl(.x, pattern="^ts")) )  )
    fake_data$predictions <- predict(ri_no_rs , newdata=fake_data, re.form=NA )
    temp1 <- fake_data$predictions
    temp2 <- mean_curve$smooth_temp
    names(temp1) <- paste0("p_", fake_data$timepoint)
    names(temp2) <- paste0("m_", mean_curve$timepoint)
    return(c(average_loss , temp1 , temp2 ))
  }
  
}

set.seed(101)
re_result <- main_data %>% my.boot.fun

boot_out_re <- boot(data=main_data , R=ifelse(run_large, 1000, 100), statistic=my.boot.fun, sim = "ordinary" , ncpus=4L, parallel="multicore")
if(run_large) {
  saveRDS(boot_out_re, "re_boot.RDS")
  re_result %>% saveRDS("re_result.RDS")
  re_result[[2]] %>% select(Case_Number, timepoint, predictions) %>% write_csv("smoothed_temp.csv") 
}

ci_holder <- NULL

used_index <- names(boot_out_re$t0 ) %>% grep(pattern="^[mps]d?_", invert=T)

for(index in  used_index) {
  if(run_large) {
    ci_holder <- rbind( ci_holder , boot.ci(boot_out_re, index =index , type=c("bca") )$bca[ ,4:5] )
  } else {
    ci_holder <- rbind( ci_holder , boot.ci(boot_out_re, index =index , type=c("perc") )$perc[ ,4:5] )
  }
}

ci_holder <- cbind(re_result[[1]][used_index], ci_holder)
colnames(ci_holder) <- c("est", "lower.ci", "upper.ci")
rownames(ci_holder) <- names(re_result[[1]])[used_index]

if(run_large) {
  ci_holder %>% round(2) %>% write.csv("main_outs.csv")
}

measured_index <- names(boot_out_re$t0 ) %>% grep(pattern="^m_")

ci_measured_holder <- NULL
for(index in  measured_index) {
  if(run_large) {
    ci_measured_holder <- rbind( ci_measured_holder , boot.ci(boot_out_re, index =index , type=c("bca") )$bca[ ,4:5] %>% round(1))
  } else {
    ci_measured_holder <- rbind( ci_measured_holder , boot.ci(boot_out_re, index =index , type=c("perc") )$perc[ ,4:5] %>% round(1))
  }
}

ci_measured_holder <- cbind(boot_out_re$t0[measured_index] , ci_measured_holder) %>% round(1) %>% as.data.frame
colnames(ci_measured_holder) <- c("est", "lower.ci", "upper.ci")
ci_measured_holder$timepoint <- names(boot_out_re$t0 ) %>% grep(pattern="^m_", value=T) %>% sub(pattern="^m_", replacement="") %>% as.numeric




table2 <- data.frame( ci_holder)
table2 <- table2[rownames(table2) %>% grep(pattern="_ex$", invert=T),]

table2$est0 <- NA
table2[ -c(1:3), "est0" ] <-  round(table2[ -c(1:3), "est" ] * nrow( main_data )  )
table2[ -c(1:3), "est1" ] <-  paste0( "(", round(table2[ -c(1:3), "est" ] * 100 ) ,"%)" )



table2[ c(1), "est0" ] <- round(table2[ c(1), "est" ] )
table2[ c(1), "est1" ] <- paste0( "(", re_result[[1]] [ re_result[[1]]  %>% names %>% grep(pattern="sd_t") ] %>% unname %>% round ,")" )

table2[ c(2:3), "est0" ] <- sprintf(round(table2[ c(2:3), "est" ] , 1 ), fmt="%1.1f")
table2[ c(2:3), "est1" ] <- paste0( "(",  re_result[[1]] [ re_result[[1]] %>% names %>% grep(pattern="sd_d") ] %>% sprintf(fmt="%1.1f") ,")" )

table2$est <- paste(table2$est0, table2$est1)

table2[-c(1:3) , "ci"] <-  paste0( table2[-c(1:3) , "lower.ci"] %>% multiply_by(100) %>% round , " to " , table2[-c(1:3) , "upper.ci"] %>% multiply_by(100) %>% round , "%")

table2[c(1) , "ci"] <-  paste0( table2[c(1) , "lower.ci"] %>% round , " to " , table2[c(1) , "upper.ci"]  %>% round )

table2[c(2:3) , "ci"] <-  paste0( table2[c(2:3) , "lower.ci"] %>% round(1) %>% sprintf(fmt="%1.1f") , " to " , table2[c(2:3) , "upper.ci"]  %>% round(1) %>% sprintf(fmt="%1.1f") )



## reorder - hard coded
table2<- table2[ c( "less_1", "hypothermia_ever" , "hypothermia_ever_low", "delta_ever" , "less_1_60", "hypothermia_60" , "hypothermia_60_low", "delta_60"  ),  ]


table2$rname <- c("Temp decrease ≥ 1°C", "Minimum Temp < 36°C", "Minimum Temp < 35°C",  "Maximum temperature decrease (°C)" , "Temp decrease ≥ 1°C at 60 min" , "Minimum Temp < 36°C at 60 min", "Minimum Temp < 35°C at 60 min" , "Maximum temp decrease (°C) at 60 min"  )

table2$lower.ci <- NULL
table2$upper.ci <- NULL
table2$est0 <- NULL
table2$est1 <- NULL

table2 %>% write_csv("table2.csv")

# fake_data <- data.frame(Case_Number=re_result[[2]]$Case_Number[1], timepoint=seq(from=0, to=max_time_in_prediction, by=1) )
# fake_data %<>% bind_cols(fake_data$timepoint %>% bs_fun %>%  as.data.frame %>% set_colnames( colnames(re_result[[2]]) %>% keep( ~grepl(.x, pattern="^ts")) )  )
# 
# fake_data$predictions <- predict(re_result[[3]] , newdata=fake_data, re.form=NA )

# fake_data %>% select(timepoint, predictions) %>% write_csv("overall_trend.csv") 

## censor smoothed estimates by the observed window

real_data <- functional_longify(main_data)

real_data <- left_join(real_data, re_result[[2]], by=c("Case_Number", "timepoint") )

## a plot of the smoothing effect
jpeg("scatter_smoothering.jpg")
plot(real_data$temper, real_data$predictions, xlab="measured temp", ylab="smoothed temp estimate" )
abline(a=0, b=1, col="red")
dev.off()

## alpha lines

jpeg("skipped.jpg")

real_data %>% group_by(timepoint) %>% summarize(smooth_temp=mean(coldnow ) ) %>% plot(type="l", lwd=4, ylim=c(35,38), col='red', xlab="time since SA", ylab="smoothed temperature")

real_data  %>%arrange(Case_Number, timepoint) %>% group_by(Case_Number) %>% group_walk( ~lines(.x$timepoint, .x$coldnow, lwd=1) )

dev.off()

mycol <- rgb(0, 0, 0 , maxColorValue = 255, alpha = 100, names = "black50")
mycol2 <- rgb(255, 0, 0 , maxColorValue = 255, alpha = 80, names = "red50")

jpeg("overlay_variation2.jpeg")

real_data %>% group_by(timepoint) %>% summarize(smooth_temp=mean(predictions ) ) %>% plot(type="l", lwd=4, ylim=c(35,38), col='red', xlab="time since SA", ylab="smoothed temperature")

real_data  %>%arrange(Case_Number, timepoint) %>%filter(is.finite(temper) ) %>% filter(Case_Number %in% seq(100)) %>% group_by(Case_Number) %>% group_walk( ~lines(.x$timepoint, .x$predictions, lwd=1, col=mycol) )

polygon( c(ci_measured_holder$timepoint, rev(ci_measured_holder$timepoint)) , c( ci_measured_holder$lower.ci , rev(ci_measured_holder$upper.ci) ) , col=mycol2, border=NA)


dev.off()

real_data %>% group_by(Case_Number) %>% summarize( severe_hypo=any(predictions < 35.0 ), delta_1=max_loss(predictions) < -1 ) -> hypothermia_indicators
hypothermia_indicators %>% write_csv("hypothermia_indicators.csv")


curve_holder <- NULL
for(index in seq(from=ncol(boot_out_re$t) - max_time_in_prediction, to = ncol(boot_out_re$t)) ) {
curve_holder <- rbind(curve_holder, boot.ci(boot_out_re, index =index , type=c("bca") )$bca[ ,4:5] ) }

curve_holder2 <- cbind(seq(nrow(curve_holder))-1,boot_out_re$t0[seq(from=ncol(boot_out_re$t) - max_time_in_prediction , to =
 ncol(boot_out_re$t))] ,curve_holder)

curve_holder2 %<>% set_colnames(c("time", "est", "lowerci", "upperci"))
curve_holder2 %>% write.csv(file="curve_with_ci.csv", row.names=FALSE)




library("tableone")
hypothermia_indicators$Case_Number %<>% as.character %<>% as.numeric %<>% as.integer
main_data$Case_Number %<>% as.character %<>% as.numeric %<>% as.integer
main_data %<>% left_join(hypothermia_indicators, by="Case_Number")

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


tab3 <- CreateTableOne(vars = c(cont_vars, cat_vars) , strata = "delta_1" , data = main_data   , factorVars = cat_vars,  includeNA = FALSE, argsExact=list(simulate.p.values=TRUE), addOverall = TRUE)

tab3 %>% print(contDigits=0, printToggle=FALSE, nonnormal=TRUE, smd=TRUE) %>% (kableExtra::kbl) %>% kableExtra::save_kable(file="table1.html")
tab3 %>% print(contDigits=0, printToggle=FALSE, nonnormal=TRUE, smd=TRUE) %>% write.csv("table1.csv")



