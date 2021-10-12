# load and install required packages
if(!require(pacman)) install.packages('pacman', dependencies = T)

p_load(tidyverse,
       naniar,
       janitor,
       readxl,
       knitr,
       kableExtra,
       expss,
       tableone,
       sjPlot,
       data.table,
       mada,
       metafor,
       lme4,
       influence.ME,
       plotrix,
       DescTools
)

#---------- Data import -----------------

#' Set working directory to where the "MA_Dataset.xlsx" file is located
#' MA_Dataset: Excel file containing all the information for all the articles

db <- readxl::read_xlsx('MA_Dataset.xlsx') 
db <- db %>%
  replace_with_na_all(condition = ~.x == 'N/A') %>%
  clean_names() %>%
  as.data.frame()

# check missing

gg_miss_var(db)

#--------Confidence interval---------------

#' Calculate CI for PPV and NPV
#' Clopper-Pearson's formula is used when sens/spec are very close to 0 or 1
#' Mercaldo's logit formula is used otherwise

data <- db %>%
  select(model_id, tp,fp,tn,fn) %>%
  na.omit()

model_id <- data$model_id
tp <- data$tp
tn <- data$tn
fp <- data$fp
fn <- data$fn
n1 <- tp+fp
n0 <- tn+fn
sens <- tp/(tp+fn)

spec <- tn/(tn+fp)
type <- ifelse(sens == 1 | spec == 1, 'clopper', 'logit')

ppv_clopper <- BinomCI(x=tp,n=n1,conf.level = 0.95, 
                       sides = 'two.sided', 
                       method = 'clopper-pearson')[,2:3]*100

npv_clopper <- BinomCI(x=tn,n=n0,conf.level = 0.95, 
                       sides = 'two.sided', 
                       method = 'clopper-pearson')[,2:3]*100

ppv_logit <- BinomCI(x=tp,n=n1,conf.level = 0.95, 
                     sides = 'two.sided', 
                     method = 'logit')[,2:3]*100

npv_logit <- BinomCI(x=tn,n=n0,conf.level = 0.95, 
                     sides = 'two.sided', 
                     method = 'logit')[,2:3]*100

ci <- cbind(model_id, type, 
            ppv_clopper, npv_clopper, 
            ppv_logit, npv_logit) 

colnames(ci) <- c('model_id', 'type',
                  'ppv_clopper_lower', 'ppv_clopper_upper',
                  'npv_clopper_lower', 'npv_clopper_upper',
                  'ppv_logit_lower', 'ppv_logit_upper',
                  'npv_logit_lower', 'npv_logit_upper')

rownames(ci) <- NULL

ci <- ci %>%
  as.data.frame() %>%
  mutate(ppv_lower = ifelse(type=='clopper', ppv_clopper_lower, ppv_logit_lower),
         ppv_upper = ifelse(type=='clopper', ppv_clopper_upper, ppv_logit_upper),
         npv_lower = ifelse(type=='clopper', npv_clopper_lower, npv_logit_lower),
         npv_upper = ifelse(type=='clopper', npv_clopper_upper, npv_logit_upper)) %>%
  select(model_id, ppv_lower,ppv_upper, npv_lower, npv_upper) %>%
  mutate_at(.vars = c('ppv_lower','ppv_upper', 'npv_lower', 'npv_upper'), as.numeric) %>%
  mutate_at(.vars = c('ppv_lower','ppv_upper', 'npv_lower', 'npv_upper'), round,2) %>%
  mutate(ppv_ci = paste(ppv_lower, '% -', ppv_upper, '%', sep=''),
         npv_ci = paste(npv_lower, '% -', npv_upper, '%', sep=''))


#' Calculate CI for PLR and NLR
#' Altman's formula is used

altman <- db %>%
  select(model_id, tp,fp,tn,fn) %>%
  # Add continuity correction - only rows with zero values
  mutate(cc = ifelse(fp==0 | fn==0, 0.5,0),
         tp = ifelse(cc==0,tp,tp+cc),
         fp = ifelse(cc==0,fp,fp+cc),
         fn = ifelse(cc==0,fn,fn+cc),
         tn = ifelse(cc==0,tn,tn+cc), 
         n1 = tp+fn, #cases
         n0 = tn+fp, # controls
         sens = tp/n1,
         spec = tn/n0) %>%
  mutate(plr = sens/(1-spec), #plr
         plr_se = sqrt(1/tp - 1/n1 + 1/fp - 1/n0),
         plr_lower = round(exp(log(plr) - 1.96 * plr_se),2),
         plr_upper = round(exp(log(plr) + 1.96 * plr_se),2),
         plr_ci = ifelse(plr_lower!='NA',paste(plr_lower, plr_upper, sep = ' - '),'NA'),
         nlr = (1-sens)/spec, #nlr
         nlr_se = sqrt(1/tn - 1/n0 + 1/fn - 1/n1),
         nlr_lower = round(exp(log(nlr) - 1.96 * nlr_se),2),
         nlr_upper = round(exp(log(nlr) + 1.96 * nlr_se),2),
         nlr_ci = ifelse(nlr_lower!='NA',paste(nlr_lower, nlr_upper, sep = ' - '), 'NA')) %>%
  select(model_id, plr, nlr, plr_lower, plr_upper, nlr_lower, nlr_upper, plr_ci, nlr_ci) 

data <- db %>% 
  select(model_id, ppv, npv) %>%
  left_join(ci, by='model_id') %>%
  left_join(altman, by='model_id')

# Export the CI table
# write.table(data,'confidence_interval.txt')

# ------- Data description ------------

data <- db %>%
  select(study, tp, roc_qpoint, with_stage_iv) %>%
  distinct(study, .keep_all=T) %>% 
  mutate(performance = ifelse(is.na(tp),'No','Yes'),
         with_stage_iv = ifelse(with_stage_iv==0, 'No','Yes')) %>%
  apply_labels(performance = 'Performance Data',
               roc_qpoint = 'ROC Data',
               with_stage_iv = 'Stage IV Cases') %>%
  select(performance, roc_qpoint, with_stage_iv) 

CreateTableOne(data = data) %>%
  print(varLabels=TRUE, showAllLevels = T) %>%
  kableone(caption = 'Data information',
           label = knitr::opts_current$get('label')) %>%
  kable_classic(full_width=F) %>%
  kable_styling(bootstrap_options = 'striped', font_size = 16, html_font = 'Cambria', position = 'left') %>%
  row_spec(row = 0, bold = TRUE, font_size = 18, extra_css = c('border-bottom: 1px solid')) %>%
  row_spec(row = 1, bold = TRUE, extra_css = 'border-bottom: 1px solid')



data <- db %>%
  select(study, year, country, avg_median_age_of_cases,
         total_samples, cases_number, control_number) %>%
  mutate_at(.vars = c('year', 'country'), as.factor) %>%
  mutate_at(.vars = 'avg_median_age_of_cases', as.numeric) %>%
  group_by(study) %>%
  mutate(m = max(total_samples)) %>%
  filter(total_samples == m) %>%
  distinct(study, .keep_all=T)

# Year of publication
table(data$year)
ggplot(data, aes(x=year)) +
  geom_bar(stat = 'count', fill='dodgerblue') +
  xlab('Year of publication')

# Country
table(data$country)
ggplot(data, aes(y=reorder(country, table(country)[country]))) +
  geom_bar(stat = 'count', fill='dodgerblue') +
  ylab('Country')

# Mean/median age of cases
round(summary(data$avg_median_age_of_cases)[c("Mean", "NA's")],2)

# Sample size cases + controls (with benign)
round(summary(data$total_samples)[c("Mean", "Min.", "Max.")],2)
# Sample size cases+controls (without benign)
round(summary(data$cases_number+data$control_number)[c("Mean", "Min.", "Max.")],2)
# Number of cases
round(summary(data$cases_number)[c("Mean", "Min.", "Max.")],2)
# Number of controls
round(summary(data$control_number)[c("Mean", "Min.", "Max.")],2)
# Number of benign
round(summary(data$total_samples-(data$cases_number+data$control_number))[c("Mean", "Min.", "Max.")],2)

# Total number of samples (with benign)
sum(data$total_samples)
# Total number of controls
sum(data$control_number)
# Total number of cases
sum(data$cases_number)
# Total number of benign cases
sum(data$total_samples-(data$cases_number+data$control_number))
# Total number of cases+controls (without benign)
sum((data$cases_number+data$control_number))

# Stage
data <- db %>%
  select(study,stage_0_i_ii_percent,stage_iii_percent, stage_iv, with_stage_iv) %>%
  mutate_at(.vars = c('stage_0_i_ii_percent','stage_iii_percent', 'stage_iv'), as.numeric) %>%
  distinct(study,.keep_all=T) %>%
  apply_labels(stage_0_i_ii_percent = 'Stage 0-I-II %',
               stage_iii_percent = 'Stage III %',
               stage_iv = 'Stage IV %')

group <- c('Overall', 'With Stage IV', 'Without Stage IV')
n <- c(dim(data)[1], table(data$with_stage_iv)['1'], table(data$with_stage_iv)['0'])

all <- data %>%
  summarise(across(.cols = c('stage_0_i_ii_percent','stage_iii_percent', 'stage_iv'),
                   list(mean = ~mean(.x,na.rm=T), miss = ~sum(is.na(.x))))) 

with <- data %>%
  filter(with_stage_iv==1) %>%
  summarise(across(.cols = c('stage_0_i_ii_percent','stage_iii_percent', 'stage_iv'),
                   list(mean = ~mean(.x,na.rm=T), miss = ~sum(is.na(.x)))))

without <- data %>%
  filter(with_stage_iv==0) %>%
  summarise(across(.cols = c('stage_0_i_ii_percent','stage_iii_percent', 'stage_iv'),
                   list(mean = ~mean(.x,na.rm=T), miss = ~sum(is.na(.x)))))


data <- cbind(group, n, rbind(all,with,without)) 
data %>%
  kbl(caption = 'Mean percentage of stages 0-I-II, III and IV of BC cases included in the studies',
      label = knitr::opts_current$get('label'), digits = 2, row.names = F,
      col.names = c('', 'n', 'Mean', 'Missing', 'Mean', 'Missing', 'Mean', 'Missing')) %>%
  kable_classic(full_width=F) %>%
  kable_styling(bootstrap_options = 'striped', font_size = 16, html_font = 'Cambria', position = 'left') %>%
  row_spec(row = 0, bold = TRUE, font_size = 18, extra_css = c('border-bottom: 1px solid')) %>%
  add_header_above(c(' '=2, 'Stage 0-I-II %' = 2, 'Stage III %' = 2, 'Stage IV %' = 2), bold = T)


# Complete data
data <- db %>% 
  select(study, model_id, tp,fp,fn,tn) %>%
  rename(TP = 'tp',
         FP = 'fp',
         FN = 'fn',
         TN = 'tn') %>%
  na.omit() 

madad(data, level=0.95)

# Preferred models - Most important model of each study
data_pref <- db %>% 
  select(study, model_id, tp,fp,fn,tn, preferred_model) %>%
  filter(preferred_model=='YES') %>%
  rename(TP = 'tp',
         FP = 'fp',
         FN = 'fn',
         TN = 'tn') %>%
  na.omit() 

# Tiff-figures - skip if no output is desired
# tiff(file = "Forest plot.tiff", width = 6400, height = 3900, units = "px", res = 400)
# par(mfrow = c(1,2))

mada::forest(madad(data_pref, level=0.95), type = "sens", 
             snames = data_pref$study, xlab = 'Sensitivity',
             main = NULL
             )

mada::forest(madad(data_pref, level=0.95), type = "spec", 
             snames = data_pref$study, 
             xlab = 'Specificity',
             main = NULL
             )

# dev.off()

rs <- rowSums(data_pref[c('TP','FP','FN','TN')])
weights <- 4 * rs / max(rs)
crosshair(data_pref, xlim = c(0, 0.6), ylim = c(0.4, 1),
          col = 1:39, lwd=weights)

ROCellipse(data_pref, pch=16, col='dodgerblue', cex=0.8)

# -------------Glmer models------------------

# Functions

#------------Logit--------------

logit <- function(x) {
  y = log(x/(1-x))
  return(y)
}

#------------Inverse logit---------

inv.logit <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}

#--------------- Complete data -------------------

# Preparation of data
data_complete <- db %>%
  select(study, model_id, tp, fp, fn, tn) %>%
  na.omit() 

# Reshaping data in long format for glmer input
reshape.data <- function(data){

  groupD <- data.frame(study=data$study, 
                       wellclassified=data$tp, 
                       misclassified=data$fn, 
                       group="disease",
                       model = 1:dim(data)[1],
                       model_id=data$model_id)

  groupH <- data.frame(study=data$study, 
                       wellclassified=data$tn, 
                       misclassified=data$fp, 
                       group="healthy",
                       model = 1:dim(data)[1],
                       model_id=data$model_id)
  
  data_long <- bind_rows(groupD, groupH) %>%
    as.data.frame()

return(data_long)

}

data_complete_long <- reshape.data(data_complete)

#' Model
#' Nested model: multiple models (if available) per study taken into account
m <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
           data=data_complete_long, 
           family = binomial(link='logit'),           
           nAGQ=0)

# Glmer output
tab.glmer <- function(m){
  
  est <- m@beta
  conf <- confint(m, method = 'Wald')[c('groupdisease', 'grouphealthy'),]
  res <- cbind(est, conf) %>%
    as.data.frame() %>%
    mutate_all(inv.logit) %>%
    mutate_all(round,2) %>%
    mutate_all(format, nsmall=2) %>%
    mutate(CI = paste('[', `2.5 %`, ' - ', `97.5 %`, ']', sep='' )) %>%
    select(est, CI)
  
  std_model <- attr(VarCorr(m)$model, 'stddev')
  corr_model <- attr(VarCorr(m)$model, 'correlation')[1,2]
  n_model <- length(unique(factor(m@flist$model)))
  
  if(is.null(VarCorr(m)$study)) {
    
    tab <- cbind(res, std_model, corr_model, n_model) %>%
      mutate_at(.vars = c('std_model', 'corr_model'), round, 2) %>%
      mutate_at(.vars = c('std_model', 'corr_model'), format, nsmall=2) 
    rownames(tab) <-  c('Sensitivity', 'Specificity')
    tab %>%
      kbl(col.names = c('Estimates', 'CI', 'Std.Dev.', 'Corr', 'n' ), digits = 2) %>%
      kable_classic(full_width = F) %>%
      kable_styling(bootstrap_options = 'striped', 
                    font_size = 16, html_font = 'Cambria', position = 'left') %>%
      row_spec(0, bold = T, extra_css = 'border-bottom: 1px solid') %>%
      add_header_above(c(' '=3, 'Model'=3), bold = T) %>%
      add_header_above(c(' '=1, 'Fixed Effects'=2, 'Random Effects' = 3), bold = T) %>%
      column_spec(1,bold = T) 
    
  } else {
    std_study <- attr(VarCorr(m)$study, 'stddev')
    corr_study <- attr(VarCorr(m)$study, 'correlation')[1,2]
    n_study <- length(unique(factor(m@flist$study)))
    
    tab <- cbind(res, std_model, corr_model, n_model, std_study, corr_study, n_study) %>%
      mutate_at(.vars = c('std_model', 'corr_model'), round, 2) %>%
      mutate_at(.vars = c('std_model', 'corr_model'), format, nsmall=2) 
    rownames(tab) <-  c('Sensitivity', 'Specificity')
    tab %>%
      kbl(col.names = c('Estimates', 'CI', 'Std.Dev.', 'Corr', 'n', 'Std.Dev.', 'Corr', 'n' ), 
          digits = 2 ) %>%
      kable_classic(full_width = F) %>%
      kable_styling(bootstrap_options = 'striped', font_size = 16,
                    html_font = 'Cambria', position = 'left') %>%
      row_spec(0, bold = T, extra_css = 'border-bottom: 1px solid') %>%
      add_header_above(c(' '=3, 'Model'=3, 'Study'=3), bold = T) %>%
      add_header_above(c(' '=1, 'Fixed Effects'=2, 'Random Effects' = 6), bold = T) %>%
      column_spec(1,bold = T) 
    
  }
  
}

#' Results
#' Estimated pooled sensitivity and pooled specificity
tab.glmer(m)

#---------------Preferred models------------------

# Preparation of data
data_pref <- db %>%
  select(study, model_id, tp, fp, tn, fn, preferred_model) %>%
  filter(preferred_model == 'YES') %>%
  na.omit()

data_pref_long <- reshape.data(data_pref)

#' Model
#' One model per study - most important model of each study
mp <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
            data=data_pref_long, 
            family = binomial(link='logit'))

#' Results
#' Estimated pooled sensitivity and pooled specificity
tab.glmer(mp)

#------------------- SROC functions---------------------------

# Extract results from glmer output
fit.glmer <- function(model, data) {
  
  fit <- list(coefficients = fixef(model),
              Psi = VarCorr(model)$model,
              Omega = VarCorr(model)$study,
              alphasens = 1,
              alphafpr = 1,
              logLik = logLik(model),
              freqdata = data[c('tp','fp','tn','fn')],
              vcov = matrix(vcov(model)@x,nrow=2),
              sens = data$tp/(data$tp+data$fn),
              fpr = data$fp/(data$fp+data$tn))
  
  # Change the sign, so that logit FPR is obtained
  fit$coefficients[2] <- - fit$coefficients[2] 
  # Change the sign 
  fit$Psi[1,2] <- fit$Psi[2,1] <- - fit$Psi[2,1]
  
  if(!is.null(fit$Omega)) {
    fit$Omega[1,2] <- fit$Omega[2,1] <- - fit$Omega[2,1]
  }
  
  fit$logLik <- NA
  attr(fit$logLik, "df") <- 5
  names(fit$freqdata) <- c('TP', 'FP', 'TN', 'FN')
  colnames(fit$vcov) <- c('groupdisease', 'grouphealthy')
  
  return(fit)
  
}

# Functions from mada
# Logit transformation
trafo <- function(alpha, x){return(talpha(alpha)$linkfun(x))}
# Inverse logit transformation
inv.trafo <- function(alpha, x){return(talpha(alpha)$linkinv(x))}

# Modified mada functions for glmer output
calc_hsroc_coef <- function(fit){
  coef <- fit$coefficients
  coef <- as.numeric(coef)
  
  Psi <- fit$Psi  
  ran.sd <- sqrt(diag(Psi))
  
  if(attr(fit$logLik,"df")== 5){
    # HSROC parameters (formulae from Harbord et al. 2008)
    Theta <- 0.5*(sqrt(ran.sd[2]/ran.sd[1])*coef[1] + sqrt(ran.sd[1]/ran.sd[2])*coef[2]) 
    Lambda <- sqrt(ran.sd[2]/ran.sd[1])*coef[1] - sqrt(ran.sd[1]/ran.sd[2])*coef[2] 
    sigma2theta <- 0.5*(ran.sd[1]*ran.sd[2] + Psi[1,2]) 
    sigma2alpha <- 2*(ran.sd[1]*ran.sd[2] - Psi[1,2])  
    beta <- log(ran.sd[2]/ran.sd[1])             
    coef_hsroc <- list(Theta = Theta,
                       Lambda = Lambda,
                       beta = beta,
                       sigma2theta = sigma2theta,
                       sigma2alpha = sigma2alpha)
    coef_hsroc <- lapply(coef_hsroc, function(x){attr(x, "names") <- NULL; x})
  }else{
    coef_hsroc = NULL
  }
  if(is.null(coef_hsroc)){
    warning("Can only compute coefficients for SROC curves without covariates. Returning NULL.")
  }
  return(coef_hsroc)
}

# SROC function for glmer output
sroc.glmer <- function(fit, fpr = 1:99/100, type = "ruttergatsonis",
                       return_function = FALSE, ...){
  stopifnot(is.logical(return_function))
  stopifnot(type %in% c("ruttergatsonis", "naive"))
  if(type == "ruttergatsonis"){
    coef_hsroc <- calc_hsroc_coef(fit)
    Lambda <- coef_hsroc$Lambda    
    Beta <- coef_hsroc$beta
    f <- function(x){
      return(inv.trafo(fit$alphasens, (Lambda*exp(-Beta/2) + exp(-Beta)*trafo(fit$alphafpr, x))))
    }
    sens <- f(fpr)
    if(!return_function){
      return(cbind(fpr, sens))
    }else{
      return(f)
    }
  }
}


# ROC ellipse function for glmer output
ROC.ellipse2 <- function(fit, conf.level = 0.95, pch = 1, add = TRUE, 
                         predict = TRUE, predlty = 3, predlwd = 1, predcol = 1, ...)
{
  alpha.sens <- fit$alphasens
  alpha.fpr <- fit$alphafpr
  mu <- fit$coefficients
  Sigma <- fit$vcov
  
  vcov <- fit$vcov
  Psi <- fit$Psi
  Omega <- fit$Omega

  talphaellipse <- ellipse(Sigma, centre = mu, level = conf.level)
  ROCellipse <- matrix(0, ncol = 2, nrow = nrow(talphaellipse))
  ROCellipse[,1] <- inv.trafo(alpha.fpr, talphaellipse[,2])
  ROCellipse[,2] <- inv.trafo(alpha.sens, talphaellipse[,1])
  
  if(predict) {
    if(is.null(Omega)) {
      Sigma_pred <- Psi + vcov
    } else {
      Sigma_pred <- Psi + Omega + vcov
    }
    talphaellipse_pred <- ellipse(Sigma_pred, centre = mu, level = conf.level)
    predellipse <- matrix(0, ncol = 2, nrow = nrow(talphaellipse_pred))
    predellipse[,1] <- inv.trafo(alpha.fpr, talphaellipse_pred[,2])
    predellipse[,2] <- inv.trafo(alpha.sens, talphaellipse_pred[,1])

  }
  
  if(add){
    lines(ROCellipse, ...)
    points(inv.trafo(alpha.fpr, mu[2]), 
           inv.trafo(alpha.sens, mu[1]), pch = pch, ...)
    lines(predellipse, lty = predlty, lwd = predlwd, col=predcol) 
   return(invisible(NULL))
  }
  if(!add){
    return(list(ROCellipse = ROCellipse, 
                fprsens = matrix(c(inv.trafo(alpha.fpr, mu[2]), 
                                   inv.trafo(alpha.sens, mu[1])),nrow = 1)))
  }
  
  }



ROCellipse.glmer <- function(x, level = 0.95, add = FALSE, pch = 1,
                             predict = FALSE, predlty = 3, predlwd = 1, predcol = 1, ...){
  ROC.ellipse2(x, conf.level = level, add = add, pch = pch,  
              predict = predict,  predlty = predlty, predlwd = predlwd, predcol = predcol, ...)
  
}

# Plot function for glmer output
plot.glmer <- function(x, extrapolate = TRUE, plotsumm = TRUE, level = 0.95, 
                       ylim = c(0,1), xlim = c(0,1), pch = 1, 
                       sroclty = 1, sroclwd = 1, 
                       predict = FALSE, 
                       predlty = 3, predlwd = 1, predcol = 1,
                       type = "ruttergatsonis",
                       ...)
{
  plot(c(2,2), ylim = ylim, xlim = xlim, 
       xlab = "False Positive Rate", ylab = "Sensitivity", ...)
  if(length(x$coefficients) == 2){
    FP <- x$freqdata$FP
    negatives <- FP + x$freqdata$TN
    FPR <- FP/negatives
    
    if(extrapolate){bound = c(0,1)}
    if(!extrapolate){bound = c(min(FPR), max(FPR))}
    srocmat <- sroc.glmer(x)
    lines(srocmat[cut(srocmat[,1],bound, "withinbound") == "withinbound",], 
          lty = sroclty, lwd = sroclwd)
  }else{
    warning("Not plotting any SROC for meta-regression")
  }
  if(plotsumm){
    ROCellipse.glmer(x, level = level, add = TRUE, pch = pch, 
                     predict = predict, predlty = predlty, predlwd = predlwd, predcol = predcol, ...)
  }
  return(invisible(NULL))
}

# Plot SROC only on observed data
sroc.glmer.range <- function(fit) {
  
  min_s <- min(fit$sens)
  max_f <- max(fit$fpr)
  
  s <- sroc.glmer(fit) %>%
    as.data.frame() %>%
    filter(fpr <= max_f & sens >= min_s)
}

# Standard plot for subgroups 
plot.glmer.subgroup <- function(fit1, fit2, main) {

  plot.glmer(fit1, xlim=c(0,1), ylim =c(0,1), pty='s', sroclty=2, predict = T,
             main = main)

  ROCellipse.glmer(fit1, lty=1, pch=16, col='red', add=TRUE, 
                   predict = T, predlty = 3, predlwd = 1, predcol = 'red')
  lines(sroc.glmer(fit1), lty=2, col='red')
  lines(sroc.glmer.range(fit1), lty=1, lwd=2, col='red')
  points(fpr(fit1$freqdata), sens(fit1$freqdata), cex=0.8, col='red',pch = 16)

  ROCellipse.glmer(fit2, lty=1, pch=16, col='blue', add=TRUE, 
                   predict = T, predlty = 3, predlwd = 1, predcol = 'blue')
  lines(sroc.glmer(fit2), lty=2, col='blue')
  lines(sroc.glmer.range(fit2), lty=1, lwd=2, col='blue')
  points(fpr(fit2$freqdata), sens(fit2$freqdata), cex=0.8, pch=16, col='blue')
  return(invisible(NULL))
  
}



#-------------- SROC---------------------

# Complete data
fit_complete <- fit.glmer(m, data_complete)

# Tiff-figures - skip if no output is desired
# tiff(file = "SROCs.tiff", width = 3600, height = 1900, units = "px", res = 400)
# par(mfrow=c(1,2))

plot.glmer(fit_complete, type = 'ruttergatsonis', 
           xlim = c(0,1), ylim = c(0,1), sroclwd = 1, sroclty = 2, predict = T,
           main = "A")
lines(sroc.glmer.range(fit_complete), lwd = 2, lty = 1)
points(fpr(fit_complete$freqdata), sens(fit_complete$freqdata), pty = 's', pch= 16, cex=0.8,
       col= factor(data_complete$study))


# Preferred model
fit_pref <- fit.glmer(mp, data_pref)

plot.glmer(fit_pref, type = 'ruttergatsonis', 
           xlim = c(0,1), ylim = c(0,1), sroclwd = 1, sroclty = 2, predict = T,
           main = "B")
lines(sroc.glmer.range(fit_pref), lwd = 2, lty = 1)
points(fpr(fit_pref$freqdata), sens(fit_pref$freqdata), pty = 's', pch= 16, cex=0.8, col = 'blue')

# dev.off()
# ------------- Subgroup analysis: all models-----------------

#----1 Plasma vs Serum-----

# Preparation of data
data <- db %>% 
  select(study, model_id, source, tp,fp,fn,tn) %>%
  na.omit() 

table(data$source)

# Plasma
plasma <- data %>%
  filter(source == 'Plasma')

plasma_long <- reshape.data(plasma)

m_plasma <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                      data=plasma_long, 
                      family = binomial(link='logit'))

tab.glmer(m_plasma)

fit_plasma <- fit.glmer(m_plasma, plasma)

# Serum
serum <- data %>%
  filter(source == 'Serum')

serum_long <- reshape.data(serum)

m_serum <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                 data=serum_long, 
                 family = binomial(link='logit'),
                 nAGQ=0)

tab.glmer(m_serum)

fit_serum <- fit.glmer(m_serum, serum)

# SROC

# tiff(file = "Subgroups_all.tiff", width = 3700, height = 3700, units = "px", res = 400)
# par(mfrow=c(2,2))

plot.glmer.subgroup(fit_plasma, fit_serum, 
                    main = "A")
legend("bottomright", c("Plasma","Serum"), pch=16, lty=1, col = c('red', 'blue'))

#----2 Single vs Multiple miRNA panels-----

# Preparation of data
data <- db %>% 
  select(study,model_id, mi_rna_panel, tp,fp,fn,tn) %>%
  na.omit()

table(data$mi_rna_panel)

# Single miRNA panel
single <- data %>%
  filter(mi_rna_panel == 'Single')

single_long <- reshape.data(single)

m_single <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                  data=single_long, 
                  family = binomial(link='logit'))

tab.glmer(m_single)

fit_single <- fit.glmer(m_single, single)

# Multiple miRNA panel
multiple <- data %>%
  filter(mi_rna_panel == 'Multiple')

multiple_long <- reshape.data(multiple)

m_multiple <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                 data=multiple_long, 
                 family = binomial(link='logit'),
                 nAGQ=0)

tab.glmer(m_multiple)

fit_multiple <- fit.glmer(m_multiple, multiple)

# SROC

plot.glmer.subgroup(fit_single, fit_serum, 
                    main = "B")
legend("bottomright", c("Single miRNA","Multiple miRNAs"), pch=16, lty=1, col = c('red', 'blue'))

#----3 Exogenous vs Endogenous normalizers----

# Preparation of data
data <- db %>% 
  select(study,model_id, normalizer_method, tp,fp,fn,tn) %>%
  na.omit()

table(data$normalizer_method)

# Endogenous normalizers
endogenous <- data %>%
  filter(normalizer_method == 'endogenous')

endogenous_long <- reshape.data(single)

m_endogenous <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                  data=endogenous_long, 
                  family = binomial(link='logit'))

tab.glmer(m_endogenous)

fit_endogenous <- fit.glmer(m_endogenous, endogenous)

# Exogenous normalizers
exogenous <- data %>%
  filter(normalizer_method == 'exogenous')

exogenous_long <- reshape.data(exogenous)

m_exogenous <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                    data=exogenous_long, 
                    family = binomial(link='logit'),
                    nAGQ=0)

tab.glmer(m_exogenous)

fit_exogenous <- fit.glmer(m_exogenous, exogenous)

# SROC

plot.glmer.subgroup(fit_endogenous, fit_exogenous, 
                    main = "C")
legend("bottomright", c("Endogenous normalizer","Exogenous normalizer"), pch=16, lty=1, col = c('red', 'blue'))


#----4 With stage IV cases vs without stage IV cases----

# Preparation of data
data <- db %>% 
  select(study,model_id, with_stage_iv, tp,fp,fn,tn) %>%
  na.omit()

table(data$with_stage_iv)

# Without stage IV cases
without <- data %>%
  filter(with_stage_iv == 0)

without_long <- reshape.data(without)

m_without <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                  data=without_long, 
                  family = binomial(link='logit'))

tab.glmer(m_without)

fit_without <- fit.glmer(m_without, without)

# With stage IV cases
with <- data %>%
  filter(with_stage_iv == 1)

with_long <- reshape.data(with)

m_with <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                    data=with_long, 
                    family = binomial(link='logit'),
                    nAGQ=0)

tab.glmer(m_with)

fit_with <- fit.glmer(m_with, with)

# SROC

plot.glmer.subgroup(fit_without, fit_with, 
                    main = "D")
legend("bottomright", c("without stage IV cases","with stage IV cases"), pch=16, lty=1, col = c('red', 'blue'))

# dev.off()

#----5 Year of publication: 2017-----
# Split (<=2017 v >2017) 

# Preparation of data
data <- db %>% 
  select(study,model_id, year, tp,fp,fn,tn) %>%
  mutate(year2017 = ifelse(year<=2017, 'before', 'after')) %>%
  na.omit()

table(data$year2017)

# Published before and in 2017

before2017 <- data %>%
  filter(year<=2017)

before2017_long <- reshape.data(before2017)

m_before2017 <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                data=before2017_long, 
                family = binomial(link='logit'),
                nAGQ=0)

tab.glmer(m_before2017)

fit_before2017 <- fit.glmer(m_before2017, before2017)

# Published after 2017

after2017 <- data %>%
  filter(year>2017)

after2017_long <- reshape.data(after2017)

m_after2017 <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                      data=after2017_long, 
                      family = binomial(link='logit'),
                      nAGQ=0)

tab.glmer(m_after2017)

fit_after2017 <- fit.glmer(m_after2017, after2017)


# SROC

plot.glmer.subgroup(fit_before2017, fit_after2017,
                    main = "Comparison of studies' models \n published before and after 2017")
legend("bottomright", c("before 2017","after 2017"), pch=16, lty=1, col = c('red', 'blue'))

#----5 Year of publication: 2014-----
# Split (<=2014 v >2014) 

# prep data
data <- db %>% 
  select(study,model_id, year, tp,fp,fn,tn) %>%
  mutate(year2014 = ifelse(year<=2014,'before', 'after')) %>%
  na.omit()

table(data$year2014)

# Published before and in 2014

before2014 <- data %>%
  filter(year<=2014)

before2014_long <- reshape.data(before2014)

m_before2014 <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                      data=before2014_long, 
                      family = binomial(link='logit'),
                      nAGQ=0)

tab.glmer(m_before2014)

fit_before2014 <- fit.glmer(m_before2014, before2014)

# Published after 2014

after2014 <- data %>%
  filter(year>2014)

after2014_long <- reshape.data(after2014)

m_after2014 <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
                     data=after2014_long, 
                     family = binomial(link='logit'),
                     nAGQ=0)

tab.glmer(m_after2014)

fit_after2014 <- fit.glmer(m_after2014, after2014)


# SROC

plot.glmer.subgroup(fit_before2014, fit_after2014, 
                    main = "Comparison of studies' models \n published before and after 2014")
legend("bottomright", c("before 2014","after 2014"), pch=16, lty=1, col = c('red', 'blue'))


#----6 miRNA-21-5p-----
# Studies with models that analysed miRNA-21-5p - all reported models

# Preparation of data
mir21 <- db %>% 
  filter(mi_rna_s == 'miR-21') %>%
  select(study, model_id, tp,fp,fn,tn) %>%
  na.omit() 

# miRNA-21-5p

mir21_long <- reshape.data(mir21)

m_mir21 <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model)+ (group-1|study), 
                 data=mir21_long, 
                 family = binomial(link='logit'), 
                 nAGQ=0)

tab.glmer(m_mir21)

fit_mir21 <- fit.glmer(m_mir21, mir21)

# tiff(file = "miR-21-5p.tiff", width = 3600, height = 1900, units = "px", res = 400)
# par(mfrow=c(1,2))

plot.glmer(fit_mir21, xlim=c(0,1), ylim =c(0,1), pty='s', sroclty = 2, predict = T,
           main = "SROC of miR-21-5p models")
lines(sroc.glmer.range(fit_mir21), lty = 1, lwd=2)
points(fpr(fit_mir21$freqdata), sens(fit_mir21$freqdata), pch = 16,cex=1,
       col = factor(mir21$study))

# If tiff was ran, move to miRNA-21-5 on preferred models to add the second plot

# ------------- Subgroup analysis: Preferred models-----------------

#----1 Plasma vs Serum-----

# Preparation of data
data <- db %>% 
  filter(preferred_model=='YES') %>%
  select(study, model_id, source, tp,fp,fn,tn) %>%
  na.omit() 

table(data$source)

# Plasma
plasma_pref <- data %>%
  filter(source == 'Plasma')

plasma_long_pref <- reshape.data(plasma_pref)

m_plasma_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                       data=plasma_long_pref, 
                       family = binomial(link='logit'))

tab.glmer(m_plasma_pref)

fit_plasma_pref <- fit.glmer(m_plasma_pref, plasma_pref)

# Serum
serum_pref <- data %>%
  filter(source == 'Serum')

serum_long_pref <- reshape.data(serum_pref)

m_serum_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                      data=serum_long_pref, 
                      family = binomial(link='logit'),
                      nAGQ=0)

tab.glmer(m_serum_pref)

fit_serum_pref <- fit.glmer(m_serum_pref, serum_pref)

# SROC
# tiff(file = "Subgroup_pref.tiff", width = 3700, height = 3700, units = "px", res = 400)
# par(mfrow=c(2,2))

plot.glmer.subgroup(fit_plasma_pref, fit_serum_pref,
                    main = "A")
legend("bottomright", c("Plasma","Serum"), pch=16, lty=1, col = c('red', 'blue'))

#----2 Single vs Multiple miRNA panels-----

# Preparation of data
data <- db %>% 
  filter(preferred_model=='YES') %>%
  select(study,model_id, mi_rna_panel, tp,fp,fn,tn) %>%
  na.omit()

table(data$mi_rna_panel)

# Single miRNA panels
single_pref <- data %>%
  filter(mi_rna_panel == 'Single')

single_long_pref <- reshape.data(single_pref)

m_single_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                       data=single_long_pref, 
                       family = binomial(link='logit'))

tab.glmer(m_single_pref)

fit_single_pref <- fit.glmer(m_single_pref, single_pref)

# Multiple miRNA panels
multiple_pref <- data %>%
  filter(mi_rna_panel == 'Multiple')

multiple_long_pref <- reshape.data(multiple_pref)

m_multiple_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                         data=multiple_long_pref, 
                         family = binomial(link='logit'),
                         nAGQ=0)

tab.glmer(m_multiple_pref)

fit_multiple_pref <- fit.glmer(m_multiple_pref, multiple_pref)

# SROC

plot.glmer.subgroup(fit_single_pref, fit_multiple_pref,
                    main = "B")
legend("bottomright", c("Single miRNA","Multiple miRNAs"), pch=16, lty=2, col = c('red', 'blue'))

#----3 Exogenous vs Endogenous normalizers----

# Preparation of data
data <- db %>% 
  filter(preferred_model=='YES') %>%
  select(study,model_id, normalizer_method, tp,fp,fn,tn) %>%
  na.omit()

table(data$normalizer_method)

# Endogenous normalizers
endogenous_pref <- data %>%
  filter(normalizer_method == 'endogenous')

endogenous_long_pref <- reshape.data(single_pref)

m_endogenous_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                           data=endogenous_long_pref, 
                           family = binomial(link='logit'))

tab.glmer(m_endogenous_pref)

fit_endogenous_pref <- fit.glmer(m_endogenous_pref, endogenous_pref)

# Exogenous normalizers
exogenous_pref <- data %>%
  filter(normalizer_method == 'exogenous')

exogenous_long_pref <- reshape.data(exogenous_pref)

m_exogenous_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                          data=exogenous_long_pref, 
                          family = binomial(link='logit'),
                          nAGQ=0)

tab.glmer(m_exogenous_pref)

fit_exogenous_pref <- fit.glmer(m_exogenous_pref, exogenous_pref)

# SROC

plot.glmer.subgroup(fit_endogenous_pref, fit_exogenous_pref,
                    main = "C")
legend("bottomright", c("Endogenous normalizer","Exogenous normalizer"), pch=16, lty=1, col = c('red', 'blue'))

#----4 With stage IV cases vs without stage IV cases----

# Preparation of data
data <- db %>% 
  filter(preferred_model=='YES') %>%
  select(study,model_id, with_stage_iv, tp,fp,fn,tn) %>%
  na.omit()

table(data$with_stage_iv)

# Without stage IV cases
without_pref <- data %>%
  filter(with_stage_iv == 0)

without_long_pref <- reshape.data(without_pref)

m_without_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                        data=without_long_pref, 
                        family = binomial(link='logit'))

tab.glmer(m_without_pref)

fit_without_pref <- fit.glmer(m_without_pref, without_pref)

# With stage IV cases
with_pref <- data %>%
  filter(with_stage_iv == 1)

with_long_pref <- reshape.data(with_pref)

m_with_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                     data=with_long_pref, 
                     family = binomial(link='logit'),
                     nAGQ=0)

tab.glmer(m_with_pref)

fit_with_pref <- fit.glmer(m_with_pref, with_pref)

# SROC

plot.glmer.subgroup(fit_without_pref, fit_with_pref,
                    main = "D")
legend("bottomright", c("without stage IV cases","with stage IV cases"), pch=16, lty=1, col = c('red', 'blue'))

# dev.off()

#----5 Year of publication: 2017-----
# Split (<=2017 v >2017) 

# Preparation of data
data <- db %>% 
  filter(preferred_model=='YES') %>%
  select(study,model_id, year, tp,fp,fn,tn) %>%
  mutate(year2017 = ifelse(year<=2017, 'before', 'after')) %>%
  na.omit()

table(data$year2017)

# Published before and in 2017

before2017_pref <- data %>%
  filter(year<=2017)

before2017_long_pref <- reshape.data(before2017_pref)

m_before2017_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                           data=before2017_long_pref, 
                           family = binomial(link='logit'),
                           nAGQ=0)

tab.glmer(m_before2017_pref)

fit_before2017_pref <- fit.glmer(m_before2017_pref, before2017_pref)

# Published after 2017

after2017_pref <- data %>%
  filter(year>2017)

after2017_long_pref <- reshape.data(after2017_pref)

m_after2017_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                          data=after2017_long_pref, 
                          family = binomial(link='logit'),
                          nAGQ=0)

tab.glmer(m_after2017_pref)

fit_after2017_pref <- fit.glmer(m_after2017_pref, after2017_pref)


# SROC

plot.glmer.subgroup(fit_before2017_pref, fit_after2017_pref,
                    main = "Comparison of studies' preferred models \n published before and after 2017")
legend("bottomright", c("before 2017","after 2017"), pch=16, lty=1, col = c('red', 'blue'))

#----5 Year of publication: 2014-----
# Split (<=2014 v >2014) 

# Preparation of data
data <- db %>% 
  filter(preferred_model=='YES') %>%
  select(study,model_id, year, tp,fp,fn,tn) %>%
  mutate(year2014 = ifelse(year<=2014, 'before', 'after')) %>%
  na.omit()

table(data$year2014)

# Published before and in 2014

before2014_pref <- data %>%
  filter(year<=2014)

before2014_long_pref <- reshape.data(before2014_pref)

m_before2014_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                           data=before2014_long_pref, 
                           family = binomial(link='logit'),
                           nAGQ=0)

tab.glmer(m_before2014_pref)

fit_before2014_pref <- fit.glmer(m_before2014_pref, before2014_pref)

# Published after 2014

after2014_pref <- data %>%
  filter(year>2014)

after2014_long_pref <- reshape.data(after2014_pref)

m_after2014_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                          data=after2014_long_pref, 
                          family = binomial(link='logit'),
                          nAGQ=0)

tab.glmer(m_after2014_pref)

fit_after2014_pref <- fit.glmer(m_after2014_pref, after2014_pref)


# SROC

plot.glmer.subgroup(fit_before2014_pref, fit_after2014_pref,
                    main = "Comparison of studies' preferred model \n published before and after 2014")
legend("bottomright", c("before 2014","after 2014"), pch=16, lty=1, col = c('red', 'blue'))

#----6 miRNA-21-5p-----
# Studies with models that analysed miRNA-21-5p - preferred models

# Preparation of data
mir21_pref <- db %>% 
  filter(mi_rna_s == 'miR-21') %>%
  select(study, model_id, tp,fp,fn,tn) %>%
  filter(model_id != 'Swellam et al. 2019_D') %>%
  na.omit() 


# miRNA-21-5p

mir21_long_pref <- reshape.data(mir21_pref)

m_mir21_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
                 data=mir21_long_pref, 
                 family = binomial(link='logit'), 
                 nAGQ=0)

tab.glmer(m_mir21_pref)

fit_mir21_pref <- fit.glmer(m_mir21_pref, mir21_pref)

plot.glmer(fit_mir21_pref, xlim=c(0,1), ylim =c(0,1), pty='s', sroclty = 2, predict = T,
           main = "SROC of miR-21-5p samples on preferred models")
lines(sroc.glmer.range(fit_mir21_pref), lty=1, lwd=2)
points(fpr(fit_mir21_pref$freqdata), sens(fit_mir21_pref$freqdata), cex=1, pch = 16,
       col = 'blue')

# dev.off()

#------Outlier analysis ------

# Preparation of data
data <- db %>%
  # Performance data with continuity correction 0.1
  select(study, model_id, c_tp, c_fp, c_fn, c_tn) %>%
  rename(tp='c_tp',
         fp='c_fp',
         fn='c_fn',
         tn='c_tn') %>%
  na.omit() %>%
  mutate(sens = tp/(tp+fn),
         spec = tn/(tn+fp),
         fpr = 1-spec,
         # odds ratio
         or = (tp/fp)/(fn/tn),
         # zcore
         z = (or - mean(or))/sd(or))

hist(data$or, breaks=50, xlab = 'Odds Ratio',  ylab = "Models",main = 'Odds Ratio Histogram on all models')
plot(data$or, data$z, xlab = 'Odds Ratio', ylab = 'Z-score', main = 'Identifying outliers')

# Outliers
data$model_id[data$z > 2]

#-------- Influence analysis---------

# Complete data

# Number of models per study
table(data_complete$study)
# Number of studies
length(unique(data_complete$study))
study <- factor(data_complete$study, label=1:37)

influence_complete <- influence(m, group = 'model')
influence_complete_study <- influence(m, group = 'study')
plot.estex(influence_complete_study,
           which='cook', 
           xlab = "Cook's distance", ylab = 'Study')

cook_study <- cooks.distance.estex(influence_complete_study)
cook_study <- as.data.frame(cook_study)
cook_study <- cook_study %>%
  rownames_to_column() %>%
  rename(study='rowname',
         cook = 'V1') %>%
  mutate(z = (cook - mean(cook))/sd(cook)) 

exclude_study <- cook_study %>%
  filter(z > 2)

# Influential study
exclude_study$study

cook <- cooks.distance.estex(influence_complete)
cook <- as.data.frame(cook)
cook <- cook %>%
  rownames_to_column() %>%
  rename(model='rowname',
         cook = 'V1') %>%
  mutate_at(.vars='model', as.numeric) %>%
  left_join(data_complete_long[c('model', 'study','model_id')], by='model') %>%
  mutate(z = (cook - mean(cook))/sd(cook)) %>%
  distinct() 

g1 <- ggplot(cook, aes(x=cook, y=model)) +
  geom_point(col=study)+
  xlab("Cook's distance")+
  ylab(" ") +
  labs(title='B')+
  theme_classic()+
  theme(panel.grid.major.x = element_line(size = .5, color='grey'),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank())

exclude <- cook %>%
  filter(z > 2)

# Influencial models
exclude$model_id

data <- data_complete_long %>%
  filter(model %in% exclude$model == F)

# Model without the influential models
m_ex <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model) + (group-1|study), 
            data=data, 
            family = binomial(link='logit'),           
            nAGQ=0)

tab.glmer(m_ex)


# Influence analysis on preferred models 

influence_pref <- influence(mp, group = 'model')
cook_pref <- cooks.distance.estex(influence_pref)
cook_pref <- as.data.frame(cook_pref)
cook_pref <- cook_pref %>%
  rownames_to_column() %>%
  rename(model='rowname',
         cook = 'V1') %>%
  mutate_at(.vars='model', as.numeric) %>%
  left_join(data_pref_long[c('model', 'study', 'model_id')], by='model') %>%
  mutate(z = (cook-mean(cook))/sd(cook)) %>%
  distinct()

g2 <- ggplot(cook_pref, aes(x=cook, y=study)) +
  geom_point(col='dodgerblue1')+
  xlab("Cook's distance")+
  ylab(" ") +
  labs(title = 'A')+
  theme_classic()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5, color = 'grey'),
        axis.title.y = element_blank())

# Tiff-figures - skip if no output is desired
# tiff(file = "Influence.tiff", width = 4700, height = 2100, units = "px", res = 400)

cowplot::plot_grid(g2,g1)

# dev.off()

exclude <- cook_pref %>%
  filter(z > 2)

# Influential preferred models (studies)
exclude$study

data <- data_pref_long %>%
  filter(model %in% exclude$model==F)

# Model without influential studies
mp_ex <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 + (group-1|model), 
            data=data, 
            family = binomial(link='logit'),           
            nAGQ=0)

tab.glmer(mp_ex)


# --------- Preference - ratios --------
# Based on case/control ratio and 
# predicted positive (TP+FP) and predicted negative (TN+FN) ratio

# 5 groups - cutpoints
data <- db %>%
  select(study, model_id, 
         tp,fp,tn,fn, tp_fp, tn_fn,
         control_number, cases_number) %>%
  mutate(arm_ratio = cases_number/control_number,
         cat_arm_ratio = cut(arm_ratio, breaks = c(0,0.4,0.8,1.2,1.6, Inf)),
         screen_ratio = tp_fp/tn_fn,
         cat_screen_ratio = cut(screen_ratio, breaks = c(0,0.4,0.8,1.2,1.6, Inf))) %>%
  na.omit() %>%
  mutate(sens = tp/(tp+fn),
         spec = tn/(tn+fp),
         fpr= 1-spec)

mean_arm <- data %>%
  group_by(cat_arm_ratio) %>%
  summarise(m_sens = mean(sens), 
            m_fpr = mean(fpr),
            sd_sens = sd(sens),
            sd_fpr = sd(fpr),
            n=n()) %>%
  mutate(a = 1.96*(sd_fpr/n),
         b = 1.96*(sd_sens/n))

mean_screen <- data %>%
  group_by(cat_screen_ratio) %>%
  summarise(m_sens = mean(sens), 
            m_fpr = mean(fpr),
            sd_sens = sd(sens),
            sd_fpr = sd(fpr),
            n=n()) %>%
  mutate(a = 1.96*(sd_fpr/n),
         b = 1.96*(sd_sens/n))

# Tiff-figures - skip if no output is desired

# tiff(file = "Pref_rat_5c.tiff", width = 4800, height = 1800, units = "px", res = 400)
# par(mfrow = c(1,3))

# tiff(file = "Pref_rat_5c_2.tiff", width = 4000, height = 2300, units = "px", res = 300)
# par(mfrow = c(1,2))

# Diagnostic accuracy and case-control-ratio

plot(data$fpr, data$sens, col=data$cat_arm_ratio,
     cex=1, pch=16, xlim=c(0,1), ylim=c(0,1),
     xlab='False Positive Rate', ylab = 'Sensitivity',
     main = "A")
points(mean_arm$m_fpr, mean_arm$m_sens, col=mean_arm$cat_arm_ratio, cex=2, pch=18)
plotrix::draw.ellipse(mean_arm$m_fpr, 
                      mean_arm$m_sens,
                      a = mean_arm$a,
                      b = mean_arm$b,
                      border = mean_arm$cat_arm_ratio, lwd=2)
legend("bottomright",  
       #cex = 1, y.intersp = 0.5, text.width=0.3,
       legend = c("design with sig. fewer cases", 'design with fewer cases', "balanced design","design with more cases", "design with sig. more cases"),
       pch = 19,
       col = 1:5) 

# Diagnostic accuracy and \nratio of positive and negative screens

plot(data$fpr, data$sens, col=data$cat_screen_ratio, 
     cex=1, pch=16, xlim=c(0,1), ylim=c(0,1),
     xlab='False Positive Rate', ylab = 'Sensitivity',
     main = "B")
points(mean_screen$m_fpr, mean_screen$m_sens, col=mean_screen$cat_screen_ratio, cex=2, pch=18)
plotrix::draw.ellipse(mean_screen$m_fpr, 
                      mean_screen$m_sens, 
                      a = mean_screen$a, 
                      b = mean_screen$b,
                      border = mean_screen$cat_screen_ratio, lwd=2)
legend("bottomright",  
       #cex = 1, y.intersp = 0.5, text.width=0.4,
       legend = c("sig. fewer total positives than negatives", "fewer total positives than negatives", "balanced positives and negatives", "more positives than negatives", "sig. more positives than negatives"),
       pch = 19,
       col = 1:5)

# Comparison of predicted screens and case/control ratio

plot(jitter(log(data$screen_ratio)) ~ jitter(log(data$arm_ratio)), 
     ylab = "jittered log pos-neg-ratio screens",
     xlab = "jittered log case-control-ratio",
     pch = 19,
     col = data$cat_screen_ratio, 
     main = "C")
abline(v=0, h = 0, lty = 2, col = "grey50")
legend("topleft", 
       #cex = 1, y.intersp = 0.3, text.width=0.8,
       legend = c("positives << negatives", "positives < negatives","positives = negatives", "positives > negatives", "positives >> negatives"),
       pch = 19,
       col = 1:5)

# dev.off()

# 3 groups - cutpoints
data <- db %>%
  select(study, model_id, 
         tp,fp,tn,fn,tp_fp, tn_fn, 
         control_number, cases_number) %>%
  mutate(arm_ratio = cases_number/control_number,
         cat_arm_ratio = cut(arm_ratio, breaks = c(0,0.7,1.3, Inf)),
         screen_ratio = tp_fp/tn_fn,
         cat_screen_ratio = cut(screen_ratio, breaks = c(0,0.7,1.3, Inf))) %>%
  na.omit() %>%
  mutate(sens = tp/(tp+fn),
         spec = tn/(tn+fp),
         fpr= 1-spec)

mean_arm <- data %>%
  group_by(cat_arm_ratio) %>%
  summarise(m_sens = mean(sens), 
            m_fpr = mean(fpr),
            sd_sens = sd(sens),
            sd_fpr = sd(fpr),
            n=n()) %>%
  mutate(a = 1.96*(sd_fpr/n),
         b = 1.96*(sd_sens/n))

mean_screen <- data %>%
  group_by(cat_screen_ratio) %>%
  summarise(m_sens = mean(sens), 
            m_fpr = mean(fpr),
            sd_sens = sd(sens),
            sd_fpr = sd(fpr),
            n=n()) %>%
  mutate(a = 1.96*(sd_fpr/n),
         b = 1.96*(sd_sens/n))

# Tiff-figures - skip if no output is desired

# tiff(file = "Pref_rat_3c.tiff", width = 4800, height = 1800, units = "px", res = 400)
# par(mfrow = c(1,3))

# tiff(file = "Pref_rat_3c_2.tiff", width = 4000, height = 2300, units = "px", res = 400)
# par(mfrow = c(1,2))

# Diagnostic accuracy and case-control-ratio

plot(data$fpr, data$sens, col=data$cat_arm_ratio,
      main = "A",
      cex=1, pch=16, xlim=c(0,1), ylim=c(0,1),
      xlab='False Positive Rate', ylab = 'Sensitivity')
 points(mean_arm$m_fpr, mean_arm$m_sens, col=mean_arm$cat_arm_ratio, cex=2, pch=18)
 plotrix::draw.ellipse(mean_arm$m_fpr, 
                       mean_arm$m_sens, 
                       a = mean_arm$a, 
                       b = mean_arm$b,
                       border = mean_arm$cat_arm_ratio, lwd=2)
 legend("bottomright",  
        #cex = 1, y.intersp = 0.5, text.width=0.3,
        legend = c('design with fewer cases', "balanced design","design with more cases"),
        pch = 19,
        col = 1:3)

# Diagnostic accuracy and ratio of positive and negative screens
 
 plot(data$fpr, data$sens, col=data$cat_screen_ratio, 
      cex=1, pch=16, xlim=c(0,1), ylim=c(0,1),
      xlab='False Positive Rate', ylab = 'Sensitivity',
      main = "B")
 points(mean_screen$m_fpr, mean_screen$m_sens, col=mean_screen$cat_screen_ratio, cex=2, pch=18)
 plotrix::draw.ellipse(mean_screen$m_fpr, 
                       mean_screen$m_sens, 
                       a = mean_screen$a, 
                       b = mean_screen$b,
                       border = mean_screen$cat_screen_ratio, lwd=2)
 legend("bottomright",  
        #cex = 1, y.intersp = 0.5, text.width=0.4,
        legend = c("fewer total positives than negatives", "balanced positives and negatives", "more positives than negatives"),
        pch = 19,
        col = 1:3)
 
# Comparison of predicted positive (TP+FP) and negative
# screens (TN+FN) and case/control ratio
 
plot(jitter(log(data$screen_ratio)) ~ jitter(log(data$arm_ratio)), 
     ylab = "jittered log pos-neg-ratio screens",
     xlab = "jittered log case-control-ratio",
     pch = 19,
     col = data$cat_screen_ratio, 
     main = "C")
abline(v=0, h = 0, lty = 2, col = "grey50")
legend("topleft", 
       #cex = 1, y.intersp = 0.3, text.width=0.8,
       legend = c("positives < negatives","positives = negatives", "positives > negatives"),
       pch = 19,
       col = 1:3)

# dev.off()

#------- Preference - Alpha & C1  -----------

# C1 - author's percieved cost of not detecting a BC patient
# The shape parameter Alpha quantifies the (a)symmetry of the study level ROC curve.
# C1 and Alpha method, graph with the performance and imabalance of proportions

data_wide <- db %>%
  select(study, model_id, picture_name, 
         control_number, cases_number, benign_number,
         starts_with('roc_sens'), starts_with('roc_1_spec'))

data_long <- melt(setDT(data_wide), 
                  measure = patterns('roc_sens', 'roc_1_spec'), 
                  value.name = c('sens', 'fpr'),
                  variable.name = 'point')

data <- data_long %>%
  mutate_at(.vars = c('sens','fpr'), as.numeric) %>%
  mutate(spec = 1 - fpr) %>%
  # Calculate total sample size (cases and controls)
  #' These studies (models) include benign in their control cohort 
  #' Check if they were accounted for in the previous analyses too (such as glmer)
  mutate(ngesamt = ifelse(model_id %in% c('Swellam et al. 2019_2_A', 
                                           'Swellam et al. 2019_2_B',
                                           'Swellam et al. 2019_2_C',
                                           'Swellam et al. 2019_2_D',
                                           'Swellam et al. 2019_2_E',
                                           'Swellam et al. 2019_2_F',
                                           'Swellam et al. 2019_2_G',
                                           
                                           'Swellam et al. 2021_A',
                                           'Swellam et al. 2019_A',
                                           'Swellam et al. 2019_B',
                                           'Swellam et al. 2019_C',
                                           'Swellam et al. 2019_D',
                                           'Swellam et al. 2019_E',
                                           'Swellam et al. 2019_F',
                                           'Swellam et al. 2019_G',
                                           
                                           'Fang_et_al_2019_B'),                                            
                          control_number + cases_number + benign_number,
                          control_number + cases_number),
         npos = cases_number, 
         nneg = ifelse(model_id %in% c('Swellam et al. 2019_2_A',
                                        'Swellam et al. 2019_2_B',
                                        'Swellam et al. 2019_2_C',
                                        'Swellam et al. 2019_2_D',
                                        'Swellam et al. 2019_2_E',
                                        'Swellam et al. 2019_2_F',
                                        'Swellam et al. 2019_2_G',
                                        
                                        'Swellam et al. 2021_A',
                                        'Swellam et al. 2019_A',
                                        'Swellam et al. 2019_B',
                                        'Swellam et al. 2019_C',
                                        'Swellam et al. 2019_D',
                                        'Swellam et al. 2019_E',
                                        'Swellam et al. 2019_F',
                                        'Swellam et al. 2019_G',
                                        
                                        'Fang_et_al_2019_B'),
                       control_number + benign_number,
                       control_number),
         # calculate TP,FN,FP,TN for every pair of (sens, 1-spec) using totals
         TP = round(npos*sens),
         FN = round(npos*(1-sens)),
         FP = round(nneg*(1-spec)), 
         TN = round(nneg*spec)) %>%
  # add continuity correction 0.1 only for rows with zero values
  mutate(cc = ifelse(FP==0 | FN==0, 0.1,0),
         TP=ifelse(cc==0, TP, TP+cc),
         TN=ifelse(cc==0, TN, TN+cc),
         FP=ifelse(cc==0, FP, FP+cc),
         FN=ifelse(cc==0, FN, FN+cc)) %>%
  # approximate 1 and 0 values for (sens, spec) to apply talpha
  mutate_at(.vars = c('sens', 'spec'), 
            function(x){
              replace(x, x<=0, 0.0001)
              replace(x, x>=1, 0.9999)
            }) %>%
  select(-c(cases_number, control_number, benign_number)) %>%
  na.omit() %>%
  # id number for studies for the following loop
  arrange(model_id) %>%
  group_by(point) %>%
  mutate(study = row_number())

# Loop on studies

# Lists of df, one for each study
study <- list() 
points_study <- list() 
alphamin_study <- list()
# n = number of studies
n <- length(unique(data$study))

for(i in 1:n){
  
  # study[[i]] contains the coordinates of 3 points extracted from ROC of study i
  study[[i]] <- data %>%
    filter(study==i) 
  
  # loop on points of the study
  
  points <- list()
  # m = number of points extracted from ROC = 3
  m <- length(study[[i]]$point)
  
  for(j in 1:m) {
    # points[[j]] contains the jth point of study i 
    points[[j]] <- study[[i]] %>%
      filter(point==j) %>%
      # for each point generate 201 alpha values
      slice(rep(1:n(), each = 201)) %>% 
      mutate(alpha = 0:200/100) %>%
      # for each alpha calculate theta = talpha(p) - talpha(q) and variance
      mutate(theta = talpha(alpha)$linkfun(sens)-talpha(alpha)$linkfun(1-spec),
             vartp = (((2-alpha) * FN + alpha * TP)^2)/(TP*FN*npos),
             vartq = (((2-alpha) * TN + alpha * FP)^2)/(FP*TN*nneg),
             vartheta = vartp + vartq)
  }
  
  # group together all the points info (depending on alpha) of the study i
  for (j in 1:(m-1)) {
    points[[j+1]] <- bind_rows(points[[j+1]], points[[j]])
    
  }
  
  # points_study[[i]] contains theta and var(theta) for each alpha for each point of study i
  points_study[[i]] <- points[[m]] %>%
    # for each alpha calculate Q and mean(theta)
    group_by(alpha) %>%
    mutate(thetabar = mean(theta),
           Q = sum((theta-thetabar)^2/vartheta)) %>%
    # find alpha for which Q is minimum
    group_by(point) %>%
    mutate(alphamin = alpha[which.min(Q)])
  
  # alphamin_study[[i]] contains only the row with minimum Q
  alphamin_study[[i]] <- points_study[[i]] %>%
    filter(alpha == alphamin) %>%
    filter(point == 1) %>%
    select(point, picture_name, model_id, study, thetabar, Q, alphamin)
  
}    

# Group together minimum Q for each study
temp <- alphamin_study
for(i in 1:(n-1)) {
  temp[[i+1]] <- bind_rows(temp[[i+1]], temp[[i]])
}    

alphamin_db <- temp[[n]] %>%
  arrange(study)
alphamin_db$point <- NULL

remove(temp, points, study, alphamin_study)

# db with cost c1 information
cost_db <- alphamin_db %>%
  left_join(db[c('model_id', 
                 'sensitivity', 
                 'specificity', 
                 'accuracy')], by='model_id') %>%
  na.omit()  %>%
  mutate_at(.vars = c('sensitivity', 'specificity'), as.numeric) %>%
  # approximate 1 and 0 values for (sens,spec) to apply talpha
  mutate_at(.vars = c('sensitivity', 'specificity'), 
            function(x){
              replace(x, x<=0, 0.0001)
              replace(x, x>=1, 0.9999)
            })


# Fix values
q <- 1 - cost_db$specificity
theta <- cost_db$accuracy
alpha <- cost_db$alphamin

# talpha (talpha_expr to distiguish from talpha in mada)
talpha_expr = expression(alpha*log(x)-(2-alpha)*log(1-x))

# Derivative of talpha
D_talpha_expr = D(talpha_expr, 'x')

# talpha(q)
x = q
talpha_val_q <- eval(eval(talpha_expr, list(alpha = alpha)), list(x = x))

# Inverse of talpha in talpha(q)-theta
inv_talpha_val <- numeric()
for(i in 1:length(alpha)) {
  inv_talpha_val[i] <- talpha(alpha[i])$linkinv(talpha_val_q[i]+theta[i])
}

p <- inv_talpha_val

# Derivative of talpha in p and q
x = p
d_talpha_val_p <- eval(eval(D_talpha_expr, list(alpha = alpha)), list(x = x))

x = q
d_talpha_val_q <- eval(eval(D_talpha_expr, list(alpha = alpha)), list(x = x))

# Add the values to cost db
cost_db$talpha_val_q <- talpha_val_q
cost_db$inv_talpha_val <- inv_talpha_val
cost_db$d_talpha_val_q <- d_talpha_val_q
cost_db$d_talpha_val_p <- d_talpha_val_p

# Fixed prevalence 
cost_db$prevalence <- 0.02

# C1 cost

preference <- cost_db %>%
  left_join(db[c('model_id', 'study','preferred_model')],by='model_id') %>%
  rename(study_num = 'study.x',
         study = 'study.y') %>%
  mutate_at(.vars = 'study', as.factor) %>%
  mutate(c1_prev = (1-prevalence)/prevalence * (d_talpha_val_p/d_talpha_val_q),
         c1 =  (d_talpha_val_p/d_talpha_val_q),
         z_c1 = (c1-1)/sd(c1),
         z_alphamin = (alphamin-mean(alphamin))/sd(alphamin),
         cat_c1 = cut(z_c1,breaks=c(-Inf,-0.8,0.8,Inf)),
         cat_alphamin = cut(z_alphamin,breaks=c(-Inf,-0.8,0.8,Inf)))

preference_pref <- cost_db %>%
  left_join(db[c('model_id', 'study','preferred_model')],by='model_id') %>%
  rename(study_num = 'study.x',
         study = 'study.y') %>%
  mutate_at(.vars = 'study', as.factor) %>%
  filter(preferred_model == 'YES') %>%
  mutate(c1_prev = (1-prevalence)/prevalence * (d_talpha_val_p/d_talpha_val_q),
         c1 =  (d_talpha_val_p/d_talpha_val_q),
         z_c1 = (c1-1)/sd(c1),
         z_alphamin = (alphamin-mean(alphamin))/sd(alphamin),
         cat_c1 = cut(z_c1,breaks=c(-Inf,-0.8,0.8,Inf)),
         cat_alphamin = cut(z_alphamin,breaks=c(-Inf,-0.8,0.8,Inf)))

# All reported models

# Alpha method
summary(preference$alphamin)
table(preference$cat_alphamin)

g3 <- ggplot(preference, aes(x=alphamin,y=log(sensitivity/specificity), col=study)) +
  geom_smooth(aes(colour=NA),
              se=F,
              colour = 'dimgrey') +
  stat_smooth(method = "loess", 
              colour = "black", 
              geom = "ribbon", 
              fill = NA,
              linetype = 'dashed') +
  geom_point() +
  ggtitle("A")+
  theme_classic()+
  theme(legend.position = 'none')+
  xlab('Alpha for Q(min)')

summary(mgcv::gam(log(sensitivity/specificity)~s(alphamin), 
                  data=preference))

spec_preference_alpha <- preference %>%
  filter(z_alphamin< -0.8)

sens_preference_alpha <- preference %>%
  filter(z_alphamin>0.8)


# C1 method
summary(preference$c1)
table(preference$cat_c1)

g4 <- ggplot(preference, aes(x=c1,y=log(sensitivity/specificity), col=study)) +
  geom_smooth(aes(colour=NA),
              se=F,
              colour = 'dimgrey') +
  stat_smooth(method = "loess", 
              colour = "black", 
              geom = "ribbon", 
              fill = NA,
              linetype = 'dashed')+
  geom_point() +
  ggtitle("B")+
  theme_classic()+
  theme(legend.position = 'none')+
  xlab('Relative perceived cost of misdiagnosis')

summary(mgcv::gam(log(sensitivity/specificity)~s(c1), 
                  data=preference))

spec_preference_c1 <- preference %>%
  filter(z_c1< -0.8)
  
sens_preference_c1 <-preference %>%
  filter(z_c1> 0.8)

table(preference$cat_alphamin, preference$cat_c1)

# Common studies
preference %>%
  filter(cat_alphamin==cat_c1) %>%
  mutate(cat = factor(cat_alphamin, labels = c('pref spec', 'no sig pref'))) %>%
  select(cat,model_id) %>%
  arrange(cat) %>%
  as.data.frame()

# Tiff-figures - skip if no output is desired
# tiff(file = "Pref_alpha_C1.tiff", width = 4200, height = 2200, units = "px", res = 400)


cowplot::plot_grid(g3,g4)

# dev.off()

# Preferred models

# Alpha method

summary(preference_pref$alphamin)
table(preference_pref$cat_alphamin)

ggplot(preference_pref, aes(x=alphamin,y=log(sensitivity/specificity))) +
  geom_smooth() +
  geom_point()+
  theme_classic()+
  theme(legend.position = 'none')+
  xlab('Alpha for Q(min)')

summary(mgcv::gam(log(sensitivity/specificity)~s(alphamin), 
                  data=preference_pref))


# C1 method
summary(preference_pref$c1)
table(preference_pref$cat_c1)

ggplot(preference_pref, aes(x=c1,y=log(sensitivity/specificity)), col=study) +
  geom_smooth() +
  geom_point()+
  theme_classic()+
  theme(legend.position = 'none')+
  xlab('Relative perceived cost of misdiagnosis')

summary(mgcv::gam(log(sensitivity/specificity)~s(c1), 
                  data=preference_pref))

table(preference_pref$cat_alphamin, preference_pref$cat_c1)

# Common studies
preference_pref %>%
  mutate(cat = factor(cat_alphamin, labels = c('pref spec','no sig pref', 'pref sens'))) %>%
  filter(cat_alphamin==cat_c1) %>%
  select(cat,model_id) %>%
  arrange(cat) %>%
  as.data.frame()

# ----- Publication bias---------

data <- db %>%
  select(study, model_id, tp,fp,fn,tn) %>%
  na.omit() 

model_id <- data$model_id
tp <- data$tp
fp <- data$fp
fn <- data$fn
tn <- data$tn
study <- data$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi


results <- bind_cols(model_id, effect_size, sample_variance) 
names(results) <- c('model_id', 'effect_size', 'sample_variance')
results <- results %>%
  left_join(db[c('model_id', 'study','number_of_study')], by='model_id')

# Check for missing values

gg_miss_var(results)

funnel(x=results$effect_size, vi = results$sample_variance, yaxis = "sei", 
       level=c(90, 95, 99), back = "white",
       shade=c("gray95", "gray55", "gray75"), 
       col = results$number_of_study, 
       refline=0, pch = 16, legend = F,)
legend("topright", c(expression("p"<"0.90"),expression("0.90<p"<="0.95"), expression("0.95<p"<="0.99")),
       pch=19, col = c("gray95", "gray55", "gray75"),
       y.intersp = 0.7, text.width=2.1, cex=1)

# Egger's test

test.egger = rma.mv(effect_size,sample_variance, mod = sample_variance, random = ~1|study)
summary(test.egger)

# Trim fill

rma <- rma(effect_size, sample_variance,
           method = 'REML', slab = study)

(tf <- trimfill(rma))

lab <- str_sub(tf$slab[1:105], start = 1, end = -3)
lab.fill <- str_sub(tf$slab[106:137], start = 1, end = 6)
lab <- append(lab, lab.fill)

lab <- factor(lab, labels = 1:39)

# Tiff-figures - skip if no output is desired
# tiff(file = "Publication_bias.tiff", width = 3700, height = 2700, units = "px", res = 400)

funnel(tf$yi, tf$vi, yaxis = "sei", 
       level=c(90, 95, 99), back = "white",
       shade=c("gray95", "gray55", "gray75"), 
       col = c(lab[1:105], rep('grey40',32)),
       refline=0, pch = 16, legend = T)

# dev.off()

#------Q points Univariate analysis-------

tab_dor <- function(mq, name) {
  
  t <- cbind(name, paste(round(mq$beta,2), ' [', round(mq$ci.lb,2),'-' ,round(mq$ci.ub,2), ']', sep = ''),
             paste(round(mq$QE,2),'  ', ifelse(mq$QEp<0.001, ' <0.001', round(mq$QEp,2)), sep='')) %>%
    as.data.frame() 
  rownames(t) <- NULL
  
  t %>% kbl(col.names = c('Subgroup', 'Pooled DOR', "Cochran's Q")) %>%
    kable_classic(full_width = F) %>%
    kable_styling(#bootstrap_options = 'striped', 
      font_size = 16, html_font = 'Cambria', position = 'left')  %>%
    row_spec(0, bold = T, extra_css = 'border-bottom: 1px solid')
}

#-----1.1 metafor DOR: All reported models--------

data <- db %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn,qc_dor) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

model_id <- data$model_id
tp <- data$tp
fp <- data$fp
fn <- data$fn
tn <- data$tn
study <- data$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

results <- bind_cols(model_id, effect_size, sample_variance) 
names(results) <- c('model_id', 'effect_size', 'sample_variance')
results <- results %>%
  left_join(db[c('model_id', 'study','number_of_study')], by='model_id')

mq = rma.mv(effect_size, sample_variance, 
             random = ~1|study)
summary(mq)
metafor::forest(mq, cex = 0.5, cex.lab = 1, annotate=F, slab=NA, xlim=c(-2,12))

# -------1.2 metafor DOR: Preferred models--------

data <- db %>%
  filter(preferred_model=='YES') %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn,qc_dor) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

model_id <- data$model_id
tp <- data$tp
fp <- data$fp
fn <- data$fn
tn <- data$tn
study <- data$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

results <- bind_cols(model_id, effect_size, sample_variance) 
names(results) <- c('model_id', 'effect_size', 'sample_variance')
results <- results %>%
  left_join(db[c('model_id', 'study','number_of_study')], by='model_id')

mq_pref <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(mq_pref)
metafor::forest(mq_pref, cex = 0.5, cex.lab = 1, annotate=F, slab=NA, xlim=c(-2,12))


#----2 Univariate subgroup analyses with metafor: All reported models-----

#----2.1 Plasma/Serum-------

data <- db %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn, source) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

table(data$source)

# Plasma

plasma_q <- data %>%
  filter(source=='Plasma')

model_id <- plasma_q$model_id
tp <- plasma_q$tp
fp <- plasma_q$fp
fn <- plasma_q$fn
tn <- plasma_q$tn
study <- plasma_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_plasma_q <- rma.mv(effect_size, sample_variance, random = ~1|study,
                     method = 'REML', slab = study)
summary(m_plasma_q)
metafor::forest(m_plasma_q, cex = 0.5, cex.lab = 1, 
                annotate=F, slab=NA, xlim=c(-2,12))

# Serum

serum_q <- data %>%
  filter(source=='Serum')

model_id <- serum_q$model_id
tp <- serum_q$tp
fp <- serum_q$fp
fn <- serum_q$fn
tn <- serum_q$tn
study <- serum_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_serum_q <- rma.mv(effect_size, sample_variance, random = ~1|study,
                    method = 'REML', slab = study)
summary(m_serum_q)
metafor::forest(m_serum_q, cex = 0.5, cex.lab = 1, 
                annotate=F, slab=NA, xlim=c(-2,12))

#-----2.2 single/multiple miRNA panels--------

data <- db %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn, mi_rna_panel) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

table(data$mi_rna_panel)

# Single miRNA panel

single_q <- data %>%
  filter(mi_rna_panel=='Single')

model_id <- single_q$model_id
tp <- single_q$tp
fp <- single_q$fp
fn <- single_q$fn
tn <- single_q$tn
study <- single_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_single_q <- rma.mv(effect_size, sample_variance, random = ~1|study,
                     method = 'REML', slab = study)
summary(m_single_q)
metafor::forest(m_single_q, cex = 0.5, cex.lab = 1, 
                annotate=F, slab=NA, xlim=c(-2,12))

# Multiple miRNA panel

multiple_q <- data %>%
  filter(mi_rna_panel=='Multiple')

model_id <- multiple_q$model_id
tp <- multiple_q$tp
fp <- multiple_q$fp
fn <- multiple_q$fn
tn <- multiple_q$tn
study <- multiple_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_multiple_q <- rma.mv(effect_size, sample_variance, random = ~1|study,
                       method = 'REML', slab = study)
summary(m_multiple_q)
metafor::forest(m_multiple_q, cex = 0.5, cex.lab = 1, 
                annotate=F, slab=NA, xlim=c(-2,12))

#-----2.3 endogenous/exogenous normalizers--------

data <- db %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn, normalizer_method) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

table(data$normalizer_method)

# Endogenous normalizers

endogenous_q <- data %>%
  filter(normalizer_method=='endogenous')

model_id <- endogenous_q$model_id
tp <- endogenous_q$tp
fp <- endogenous_q$fp
fn <- endogenous_q$fn
tn <- endogenous_q$tn
study <- endogenous_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_endogenous_q <- rma.mv(effect_size, sample_variance, random = ~1|study,
                         method = 'REML', slab = study)
summary(m_endogenous_q)
metafor::forest(m_endogenous_q, cex = 0.5, cex.lab = 1, 
                annotate=F, slab=NA, xlim=c(-2,12))

# Exogenous normalizers

exogenous_q <- data %>%
  filter(normalizer_method=='exogenous')

model_id <- exogenous_q$model_id
tp <- exogenous_q$tp
fp <- exogenous_q$fp
fn <- exogenous_q$fn
tn <- exogenous_q$tn
study <- exogenous_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_exogenous_q <- rma.mv(effect_size, sample_variance, random = ~1|study,
                        method = 'REML', slab = study)
summary(m_exogenous_q)
metafor::forest(m_exogenous_q, cex = 0.5, cex.lab = 1, 
                annotate=F, slab=NA, xlim=c(-2,12))

#-----2.4 with or without stage IV cases--------

data <- db %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn, with_stage_iv) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

table(data$with_stage_iv)

# with stage IV cases

with_q <- data %>%
  filter(with_stage_iv==1)

model_id <- with_q$model_id
tp <- with_q$tp
fp <- with_q$fp
fn <- with_q$fn
tn <- with_q$tn
study <- with_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_with_q <- rma.mv(effect_size, sample_variance, random = ~1|study,
                   method = 'REML', slab = study)
summary(m_with_q)
metafor::forest(m_with_q, cex = 0.5, cex.lab = 1, 
                annotate=F, slab=NA, xlim=c(-2,12))

# without stage iV cases

without_q <- data %>%
  filter(with_stage_iv==0)

model_id <- without_q$model_id
tp <- without_q$tp
fp <- without_q$fp
fn <- without_q$fn
tn <- without_q$tn
study <- without_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_without_q <- rma.mv(effect_size, sample_variance, random = ~1|study,
                      method = 'REML', slab = study)
summary(m_without_q)
metafor::forest(m_without_q, cex = 0.5, cex.lab = 1, 
                annotate=F, slab=NA, xlim=c(-2,12))


#----3 Univariate subgroup analysis with metafor: Preferred models-----

#----3.1 Plasma/Serum-------

data <- db %>%
  filter(preferred_model=='YES') %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn, source) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

table(data$source)

# Plasma

plasma_q <- data %>%
  filter(source=='Plasma')

model_id <- plasma_q$model_id
tp <- plasma_q$tp
fp <- plasma_q$fp
fn <- plasma_q$fn
tn <- plasma_q$tn
study <- plasma_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_plasma_q <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(m_plasma_q)
metafor::forest(m_plasma_q)

# Serum

serum_q <- data %>%
  filter(source=='Serum')

model_id <- serum_q$model_id
tp <- serum_q$tp
fp <- serum_q$fp
fn <- serum_q$fn
tn <- serum_q$tn
study <- serum_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_serum_q <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(m_serum_q)
metafor::forest(m_serum_q)

#-----3.2 Single/Multiple miRNA panel--------

data <- db %>%
  filter(preferred_model=='YES') %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn, mi_rna_panel) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

table(data$mi_rna_panel)

# Single miRNA panel

single_q <- data %>%
  filter(mi_rna_panel=='Single')

model_id <- single_q$model_id
tp <- single_q$tp
fp <- single_q$fp
fn <- single_q$fn
tn <- single_q$tn
study <- single_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_single_q <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(m_single_q)
metafor::forest(m_single_q)

# Multiple miRNA panel

multiple_q <- data %>%
  filter(mi_rna_panel=='Multiple')

model_id <- multiple_q$model_id
tp <- multiple_q$tp
fp <- multiple_q$fp
fn <- multiple_q$fn
tn <- multiple_q$tn
study <- multiple_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_multiple_q <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(m_multiple_q)
metafor::forest(m_multiple_q)

#-----3.3 Endogenous/Exogenous normalizers--------

data <- db %>%
  filter(preferred_model=='YES') %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn, normalizer_method) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

table(data$normalizer_method)

# Endogenous normalizers

endogenous_q <- data %>%
  filter(normalizer_method=='endogenous')

model_id <- endogenous_q$model_id
tp <- endogenous_q$tp
fp <- endogenous_q$fp
fn <- endogenous_q$fn
tn <- endogenous_q$tn
study <- endogenous_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_endogenous_q <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(m_endogenous_q)
metafor::forest(m_endogenous_q)

# Exogenous normalizers

exogenous_q <- data %>%
  filter(normalizer_method=='exogenous')

model_id <- exogenous_q$model_id
tp <- exogenous_q$tp
fp <- exogenous_q$fp
fn <- exogenous_q$fn
tn <- exogenous_q$tn
study <- exogenous_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_exogenous_q <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(m_exogenous_q)
metafor::forest(m_exogenous_q)

#-----3.4 With or without stage IV cases--------

data <- db %>%
  filter(preferred_model=='YES') %>%
  select(study, model_id, q_tp,q_fp,q_fn,q_tn, with_stage_iv) %>%
  rename(tp = 'q_tp',
         fp = 'q_fp',
         fn = 'q_fn',
         tn = 'q_tn') %>%
  na.omit() 

table(data$with_stage_iv)

# With stage IV cases

with_q <- data %>%
  filter(with_stage_iv==1)

model_id <- with_q$model_id
tp <- with_q$tp
fp <- with_q$fp
fn <- with_q$fn
tn <- with_q$tn
study <- with_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_with_q <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(m_with_q)
metafor::forest(m_with_q)

# Without stage IV cases

without_q <- data %>%
  filter(with_stage_iv==0)

model_id <- without_q$model_id
tp <- without_q$tp
fp <- without_q$fp
fn <- without_q$fn
tn <- without_q$tn
study <- without_q$study

effect_size <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$yi
sample_variance <- escalc(measure = "OR", ai = tp, ci = fp, bi = fn, di = tn)$vi

m_without_q <- rma.uni(effect_size, sample_variance, method = 'REML', slab = study)
summary(m_without_q)
metafor::forest(m_without_q)

#------Fixed effects analysis-------

#----All reported models-------

data_adjusted <- db %>%
  select(model_id,
         source, mi_rna_panel, 
         normalizer_method, with_stage_iv) %>%
  filter(source %in% c('Plasma', 'Serum'),
         mi_rna_panel %in% c('Single', 'Multiple'),
         normalizer_method %in% c('endogenous','exogenous')) %>%
  na.omit() %>%
  mutate_at(.vars = 'with_stage_iv', as.factor)

data_adjusted_long <- data_complete_long %>%
  inner_join(data_adjusted, by='model_id')

#' Model
#' Nested model: more models per study taken into account
m_adj <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 
               + source 
               + mi_rna_panel
               + normalizer_method
               + with_stage_iv
               + (group-1|model) + (group-1|study), 
           data=data_adjusted_long, 
           family = binomial(link='logit'),           
           nAGQ=0)
tab_model(m_adj)
# Details on model with fixed effects - comparison between binary groups in fixed effects
summary(m_adj)

#-----Preferred models-------

data_adjusted <- db %>%
  select(model_id,
         source, mi_rna_panel, 
         normalizer_method, with_stage_iv) %>%
  filter(source %in% c('Plasma', 'Serum'),
         mi_rna_panel %in% c('Single', 'Multiple'),
         normalizer_method %in% c('endogenous','exogenous')) %>%
  na.omit() %>%
  mutate_at(.vars = 'with_stage_iv', as.factor)

data_adjusted_pref_long <- data_pref_long %>%
  inner_join(data_adjusted, by='model_id')

#' Model
#' Nested model: multiple models per study taken into account
m_adj_pref <- glmer(formula = cbind(wellclassified, misclassified) ~ group-1 
               + source 
               + mi_rna_panel
               + normalizer_method
               + with_stage_iv
               + (group-1|model) + (group-1|study), 
               data=data_adjusted_pref_long, 
               family = binomial(link='logit'),           
               nAGQ=0)

# Details on model with fixed effects - comparison between binary groups in fixed effects
summary(m_adj_pref)

