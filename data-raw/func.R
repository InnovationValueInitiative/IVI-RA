# STANDARD ERROR GIVEN NORMDAL DISTRIBUTION AND BOUNDS FROM CONFIDENE INTERVAL -
se_normal <- function(l, u, p = .975){
  return((u - l)/(2 * qnorm(p)))
}

# FIND BEST PARAMETRIC FIT TO SURVIVAL CURVE -----------------------------------
# Loop over distributions fitting parametric survival model
parametric_fit <- function(dists, data){
  fit <- vector(length(dists), mode = "list"); names(fit) <- dists
  for (i in 1:length(dists)){
    fit[[i]] <- flexsurvreg(Surv(t, event) ~ 1,
                            data = data, dist = dists[i])
  }
  return(fit)
}

# Estimate either AIC or BIC
IC <- function(mod, type){
  ret <- if (type == "AIC") AIC(mod) else {BIC(mod)}
  return(ret)
}

# Matrix containing AIC/BIC for each distribution
ic_mat <- function(fits, modnames){
  ic <- matrix(NA, nrow = length(fits), ncol = 2)
  rownames(ic) <- modnames
  colnames(ic) <- ic.type <- c("AIC", "BIC")
  for (i in 1:ncol(ic)){
    for (m in 1:length(modnames)){
      ic[m, i] <- IC(fits[[m]], ic.type[i])
    }
  }
  return(ic)
}

# SAMPLE FROM CUMULATIVE HAZARD ------------------------------------------------
#Sample survival data from a distribution given by a cumulative hazard 
#
# Args: 
#   time Vector containing time points to sample from.
#   Haz The cumulative hazard at those time points
#   size The size of the sample
#   replace Sample with or without replacement? (default=TRUE, with replacement) 
#
#  Return:
#   A random sample from Haz
cumhaz_sample <- function(time, Haz, size = 1, replace = TRUE){
  p <- diff(c(0, 1 - exp(-Haz)))
  p <- c(p, exp(-Haz[length(Haz)])) # add probability of sampling time=Inf 
  return(sample(c(time, Inf), size = size, prob = p, replace = replace))
}

# FIT SPLINES TO RECONSTRUCTED SURVIVAL DATA -----------------------------------
surv_splines <- function(data, time, knots = 4){
  km <- survfit(Surv(t, event) ~1, data = data)
  spline <- flexsurvspline(Surv(t, event) ~ 1, data = data, k = knots,
                               scale = "hazard")
  spline.surv <- summary(spline, t = time, type = "survival", ci = FALSE)[[1]]
  pdat <- data.frame(t = c(km$time, time),
                         surv = c(km$surv, spline.surv$est),
                         lab = rep(c("KM", "Spline"), 
                                   times = c(length(km$surv), length(time))))
  p <- ggplot(pdat, aes(x = t, y = surv, col = lab)) + geom_line() + xlab("Time") +
    ylab("Survival") + scale_colour_discrete("") +
    theme(legend.position = "bottom")
  spline.haz <- summary(spline, t = time[-1], type = "hazard", ci = FALSE)[[1]]
  return(list(p = p, time = time, haz =  spline.haz[, "est"], 
              surv = spline.surv[, "est"], km = km))
}

# KAPLAN MEIER PLOT ------------------------------------------------------------
km_plot <- function(x, timelab = "Years", convert_years = TRUE){
  dat <- data.table(time = x$time, surv = x$surv)
  if (convert_years == TRUE){
    dat$time <- dat$time/12
  }
  p <- ggplot(dat, aes(x = time, y = surv)) + geom_line() + xlab(timelab) +
    ylab("Survival")
  return(p)
}

# SURVIAL DATA TO COMPARE PARAMETRIC FITS WITH NON-PARAMETRIC ------------------
# compare km to survival data
survfit_data <- function(km, fits, modnames, long = TRUE){
  surv <- data.table(time =  km$time, surv = km$surv,
                     mod = "Kaplan-Meier")
  for (i in 1:length(modnames)){
    surv.tmp <- summary(fits[[i]], type = "surv", t = km$time)[[1]]
    colnames(surv.tmp)[1:2] <- c("time", "surv")
    surv.tmp$mod <- modnames[i]
    surv <- rbind(surv, surv.tmp[, c("time", "surv", "mod")])
  }
  surv$cumhaz <- -log(surv$surv)
  if (long == TRUE){
    km.indx <- which(surv$mod == "Kaplan-Meier")
    par.indx <- which(surv$mod != "Kaplan-Meier")
    surv <- surv[c(rep(km.indx, times = length(modnames)), par.indx)]
    surv[, type := ifelse(mod == "Kaplan-Meier", mod, "Parametric")]
    km.indx.new <- which(surv$mod == "Kaplan-Meier")
    surv[km.indx.new, mod := rep(modnames, each = length(km.indx))]
  }
  return(surv)
}

# compare spline to survival data
survfit_data_spline <- function(time, spline, fits, modnames){
  surv <- data.frame(time =  time, surv = spline,
                     mod = "spline")
  for (i in 1:length(modnames)){
    surv.tmp <- summary(fits[[i]], type = "surv", t = time)[[1]]
    colnames(surv.tmp)[1:2] <- c("time", "surv")
    surv.tmp$mod <- modnames[i]
    surv <- rbind(surv, surv.tmp[, c("time", "surv", "mod")])
  }
  surv$cumhaz <- -log(surv$surv)
  return(surv)
}

# PARAMETRIC FIT DIAGNOSTICS ---------------------------------------------------
parametric_fits_diag <- function(km, parfits, modlabs){
  survfit.data <- survfit_data(km, parfits, modlabs)
  p.survfits.gg <- ggplot(survfit.data[mod == "Generalized gamma"], 
                                 aes(x = time/12, y = surv, col = type)) + geom_line() +
    xlab("Years") + ylab("Survival") + scale_colour_discrete("") +
    theme(legend.position="bottom")
  p.survfits <- ggplot(survfit.data, 
                              aes(x = time/12, y = surv, col = type)) + geom_line() +
    facet_wrap("mod", scales = "free_x") + xlab("Years") + ylab("Survival") + 
    scale_colour_discrete("") + theme(legend.position="bottom")
  p.cumhazfits <- ggplot(survfit.data, 
                                aes(x = time/12, y = cumhaz, col = type)) + geom_line() +
    facet_wrap("mod", scales = "free_x") + xlab("Years") + ylab("Survival") + 
    scale_colour_discrete("") + theme(legend.position="bottom")
  return(list(data = survfit.data, p.survfits.gg = p.survfits.gg, p.survfits = p.survfits,
              p.cumhazfits = p.cumhazfits))
}

# COMPARE TREATMENT DURATION FOR MODERATE AND GOOD EULAR RESPONDERS ------------
survfit_by_eular <- function(parfits.m, parfits.g, t, dist = "gengamma"){
  survfit.data <-  rbind(data.table(summary(parfits.m[[dist]], t = t)[[1]],
                             eular = "Moderate"),
                         data.table(summary(parfits.g[[dist]], t = t)[[1]],
                                   eular = "Good"))
  p.survfit <- ggplot(survfit.data, aes(x = time/12, y = est, col = eular)) + 
    geom_line() +  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = eular, linetype = NA), alpha = 0.5) +
    xlab("Year since initiating treatment") + ylab("Proportion surviving") +
    theme(legend.position = "bottom") +
    scale_fill_discrete(name = "Eular response") +
    scale_colour_discrete(guide = FALSE)
  return(list(data = survfit.data, p = p.survfit))
}

# SAMPLE SURVIVAL DATA FROM CUMULATIVE HAZARD ----------------------------------
# Sample survival data from a distribution given by a cumulative hazard 
#
# Args:
#   time: Vector containing time points to sample from.
#   Haz: The cumulative hazard at those time points.
#   size: The size of the sample.
#   replace: Sample with or without replacement? (default=TRUE, with replacement)
#
# Returns:
#   A random sample (a vector) of survival probabilties from Haz.
cumhaz_sample <- function(time, Haz, size = 1, replace = TRUE){
  p <- diff(c(0, 1 - exp(-Haz)))
  p <- c(p, exp(-Haz[length(Haz)])) # add probability of sampling time=Inf 
  return(sample(c(time, Inf), size = size, prob = p, replace = replace))
}

# EXTRACT FLEXSURVREG MODEL PARAMETERS -----------------------------------------
# Extract parameters from a survival model fit using flexsurvreg
#
# Args:
#   x: A model fit using flexsurvreg.
#
# Returns:
#   A list to be used in sample_pars.
flexsurvreg_pars <- function(x){
  d <- x$dlist
  basepar.loc.indx <- which(d$pars == d$location)
  basepar.anc1.indx <- which(d$pars != d$location)[1]
  basepar.anc2.indx <-  which(d$pars != d$location)[2]
  loc.indx <- c(x$basepars[basepar.loc.indx], x$covpars[x$mx[[d$location]]])
  anc1.indx <- c(x$basepars[basepar.anc1.indx], x$covpars[x$mx[[d$pars[basepar.anc1.indx]]]])
  anc2.indx <- c(x$basepars[basepar.anc2.indx], x$covpars[x$mx[[d$pars[basepar.anc2.indx]]]])
  l <- list(est = x$coef, vcov = x$cov, loc.index = loc.indx, 
           anc1.index = anc1.indx, anc2.index = anc2.indx)
  return(l)
}

# RECONSTRUCT INDIVIDUAL PATIENT DATA ------------------------------------------
# Reconstruct individual-level patient data from digitized curve.
#
# Args: 
#   sc_data: Dataframe of survival curve data. Should be two columns. First column is 
#            is time and second column is survival probability.
#   risk_date: Data frame of two columns representing number at risk. First column is
#              time and second column is number at risk. These correspond to the times and 
#              number at risk frequently reported under Kaplan Meier curves in journal articles.
#   arm_id Choose id number for arm.
#   tot_events Total number of events (e.g. deaths) at end of study time. Default is NULL.
#   
# Returns:
#    Dataframe with three columns. Fist column, t, is time to event. Second column, event,
#    is a binary indicator = 1 if an event occured and 0 if an individual was censored. Third
#    column, arm, is the id for the trial arm. 
reconstruct_ipd <- function(sc_data, risk_date, arm_id, tot_events = NULL){
  
  # interval columns
  sc_data <- cbind(seq(nrow(sc_data)), sc_data)
  risk_date <- cbind(seq(nrow(risk_date)), risk_date)
  
  # name dataset columns
  colnames(sc_data) <- c("k", "tk", "sk")
  colnames(risk_date) <- c("i", "trisk", "nrisk")
  
  # coordinates
  t.S <- sc_data[, "tk"]
  S <- sc_data[, "sk"]
  
  # number at risk
  risk_date$lower <- NA
  risk_date$upper <- NA
  if (nrow(risk_date) > 1){
    for (i in 1:(nrow(risk_date) - 1)){
      risk_date$upper[i] <- max(sc_data[sc_data$tk < risk_date$trisk[i + 1], 1])
    }
  }
  risk_date$upper[nrow(risk_date)] <- sc_data[nrow(sc_data), "k"]
  risk_date$lower <- c(1, risk_date$upper[-length(risk_date$upper)] + 1)
  n.int <- nrow(risk_date)
  if (risk_date[n.int, "lower"] >= risk_date[n.int, "upper"]){
    risk_date[n.int, "lower"] <- risk_date[n.int, "lower"] - 1
    risk_date[n.int - 1, "upper"] <- risk_date[n.int - 1, "upper"] - 1
  }
  
  t.risk <- risk_date[, "trisk"]
  n.risk <- risk_date[, "nrisk"]
  lower <- risk_date[, "lower"]
  upper <- risk_date[, "upper"]
  n.t <- upper[n.int]
  
  # initialise vectors
  arm <- rep(arm_id, n.risk[1])
  n.censor <- rep(0, (n.int - 1))
  n.hat <- rep(n.risk[1] + 1, n.t)
  cen <- rep(0, n.t)
  d <- rep(0, n.t)
  KM.hat<-rep(1,n.t)
  last.i <- rep(1, n.int)
  sumdL <- 0
  
  if (n.int > 1){
    #Time intervals 1,...,(n.int-1)
    for (i in 1:(n.int-1)){
      #First approximation of no. censored on interval i
      n.censor[i]<- round(n.risk[i]*S[lower[i+1]]/S[lower[i]]- n.risk[i+1])
      #Adjust tot. no. censored until n.hat = n.risk at start of interval (i+1)
      while((n.hat[lower[i+1]]>n.risk[i+1])||((n.hat[lower[i+1]]<n.risk[i+1])&&(n.censor[i]>0))){
        if (n.censor[i]<=0){
          cen[lower[i]:upper[i]]<-0
          n.censor[i]<-0
        }
        if (n.censor[i]>0){
          cen.t<-rep(0,n.censor[i])
          for (j in 1:n.censor[i]){
            cen.t[j]<- t.S[lower[i]] +
              j*(t.S[lower[(i+1)]]-t.S[lower[i]])/(n.censor[i]+1)
          }
          #Distribute censored observations evenly over time. Find no. censored on each time interval.
          cen[lower[i]:upper[i]]<-hist(cen.t,breaks=t.S[lower[i]:lower[(i+1)]],
                                       plot=F)$counts
        }
        #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
        n.hat[lower[i]]<-n.risk[i]
        last<-last.i[i]
        for (k in lower[i]:upper[i]){
          if (i==1 & k==lower[i]){
            d[k]<-0
            KM.hat[k]<-1
          }
          else {
            d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
            KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          }
          n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
          if (d[k] != 0) last<-k
        }
        n.censor[i]<- n.censor[i]+(n.hat[lower[i+1]]-n.risk[i+1])
      }
      if (n.hat[lower[i+1]]<n.risk[i+1]) n.risk[i+1]<-n.hat[lower[i+1]]
      last.i[(i+1)]<-last
    }
  }
  #Time interval n.int.
  if (n.int > 1){
    #Assume same censor rate as average over previous time intervals.
    n.censor[n.int]<- min(round(sum(n.censor[1:(n.int-1)])*(t.S[upper[n.int]]-
                                                              t.S[lower[n.int]])/(t.S[upper[(n.int-1)]]-t.S[lower[1]])), n.risk[n.int])
  }
  if (n.int == 1){n.censor[n.int] <- 0}
  if (n.censor[n.int] <= 0){
    cen[lower[n.int]:(upper[n.int] - 1)] <- 0
    n.censor[n.int] <- 0
  }
  if (n.censor[n.int]>0){
    cen.t<-rep(0,n.censor[n.int])
    for (j in 1:n.censor[n.int]){
      cen.t[j]<- t.S[lower[n.int]] +
        j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
    }
    cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                             plot=F)$counts
  }
  #Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
  n.hat[lower[n.int]]<-n.risk[n.int]
  last<-last.i[n.int]
  for (k in lower[n.int]:upper[n.int]){
    if(KM.hat[last] !=0){
      d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))} else {d[k]<-0}
    KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
    n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
    #No. at risk cannot be negative
    if (n.hat[k+1] < 0) {
      n.hat[k+1]<-0
      cen[k]<-n.hat[k] - d[k]
    }
    if (d[k] != 0) last<-k
  }
  #If total no. of events reported, adjust no. censored so that total no. of events agrees.
  if (!is.null(tot_events)){
    if (n.int>1){
      sumdL<-sum(d[1:upper[(n.int-1)]])
      #If total no. events already too big, then set events and censoring = 0 on all further time intervals
      if (sumdL >= tot_events){
        d[lower[n.int]:upper[n.int]]<- rep(0,(upper[n.int]-lower[n.int]+1))
        cen[lower[n.int]:(upper[n.int]-1)]<- rep(0,(upper[n.int]-lower[n.int]))
        n.hat[(lower[n.int]+1):(upper[n.int]+1)]<- rep(n.risk[n.int],(upper[n.int]+1-lower[n.int]))
      }
    }
    #Otherwise adjust no. censored to give correct total no. events
    if ((sumdL < tot_events)|| (n.int==1)){
      sumd<-sum(d[1:upper[n.int]])
      while ((sumd > tot_events)||((sumd< tot_events)&&(n.censor[n.int]>0))){
        n.censor[n.int]<- n.censor[n.int] + (sumd - tot_events)
        if (n.censor[n.int]<=0){
          cen[lower[n.int]:(upper[n.int]-1)]<-0
          n.censor[n.int]<-0
        }
        if (n.censor[n.int]>0){
          cen.t<-rep(0,n.censor[n.int])
          for (j in 1:n.censor[n.int]){
            cen.t[j]<- t.S[lower[n.int]] +
              j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
          }
          cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                                   plot=F)$counts
        }
        n.hat[lower[n.int]]<-n.risk[n.int]
        last<-last.i[n.int]
        for (k in lower[n.int]:upper[n.int]){
          d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
          KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          if (k != upper[n.int]){
            n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
            #No. at risk cannot be negative
            if (n.hat[k+1] < 0) {
              n.hat[k+1]<-0
              cen[k]<-n.hat[k] - d[k]
            }
          }
          if (d[k] != 0) last<-k
        }
        sumd<- sum(d[1:upper[n.int]])
      }
    }
  }
  #write.table(matrix(c(t.S,n.hat[1:n.t],d,cen),ncol=4,byrow=F),paste(path,KMdatafile,sep=""),sep="\t")
  ### Now form IPD ###
  #Initialise vectors
  t.IPD<-rep(t.S[n.t],n.risk[1])
  event.IPD<-rep(0,n.risk[1])
  #Write event time and event indicator (=1) for each event, as separate row in t.IPD and event.IPD
  k=1
  for (j in 1:n.t){
    if(d[j]!=0){
      t.IPD[k:(k+d[j]-1)]<- rep(t.S[j],d[j])
      event.IPD[k:(k+d[j]-1)]<- rep(1,d[j])
      k<-k+d[j]
    }
  }
  #Write censor time and event indicator (=0) for each censor, as separate row in t.IPD and event.IPD
  for (j in 1:(n.t-1)){
    if(cen[j]!=0){
      t.IPD[k:(k+cen[j]-1)]<- rep(((t.S[j]+t.S[j+1])/2),cen[j])
      event.IPD[k:(k+cen[j]-1)]<- rep(0,cen[j])
      k<-k+cen[j]
    }
  }
  
  IPD <- matrix(c(t.IPD, event.IPD, arm), ncol = 3, byrow = F)
  IPD <- data.frame(t = t.IPD, event = event.IPD, arm = arm)
  return(IPD)
}

