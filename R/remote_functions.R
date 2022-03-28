# Installing uninstalled packages and loading packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


# Reading in .sav data
read_in <- function(df = files, directory) {
  for (fname in df) {
    df_list[[fname]] <- haven::read_sav(paste0(directory ,fname), encoding = NULL, user_na = FALSE, col_select = NULL,skip = 0, n_max = Inf, .name_repair = "unique")
  }
  names(df_list) <- paste0("", gsub(".sav","",names(df_list)))
  ff <- df_list
}

# Creating nicely formatted correlation table
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower", "none"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    ")))
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  else if(removeTriangle[1]=="none"){
    Rnew <- as.matrix(Rnew)
    Rnew <- as.data.frame(Rnew)
  }
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex")
  }
}


corstars_no_stars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower", "none"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    ")))
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R,sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  else if(removeTriangle[1]=="none"){
    Rnew <- as.matrix(Rnew)
    Rnew <- as.data.frame(Rnew)
  }
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex")
  }
}


# Getting F and corresponding p value from lm object
Regressionp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)}


# Pasting together two character vectors; if corresponding element in either vector is NA, return only element from
# vector with !NA. If both elements are NA, return NA.
paste3 <- function(...,sep="") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}
val<- paste3(c("a","b", "c", NA), c("A","B", NA, NA))
val

# Saving indirect effect of mediation model
indirectsaved <- function(comp_df, random) {
  d <- comp_df[random, ] #rndomize by row
  apath <- lm(NA. ~ finstrain_center, data = d)
  bpath <- lm(satis ~ NA., data = d)
  indirect <- apath$coefficients[2]* bpath$coefficients[2]
  return(indirect)
}

# Creating time since business foundation in days
timediff <- function() {
  year_n <- as.numeric(names(attr(obs_df$T1year_1,"labels")[match(obs_df$T1year_1,attr(obs_df$T1year_1,"labels"))]))
  month_n <- as.numeric(obs_df$T1month_1 %>% remove_all_labels(.))
  time <- as.yearmon(paste(year_n, month_n), "%Y %m") %>% as_date(.)
  record <- as.yearmon(obs_df$T1RecordedDate, "%Y %m") %>% as_date(.)
  difftime <- difftime(record, time, UTC,
                       units = c("days"))
  return(difftime)
}

# Creating attention fail output T1. This may then be used to remove attention fails.
attention_fail_t1 <- function() {
  attention_f <- obs_df[(obs_df$T1chal_3 != 4) | (obs_df$T1threat_4 != 4) | (obs_df$T1satis_4 != 2),]
  ID_vals <- data.frame(table(attention_f$T1ID))
  Rows_fails <- attention_f$T1ID %in% ID_vals[ID_vals$Freq > 0,1]
  Att_fails <- attention_f[Rows_fails,]
  return(Att_fails)
}


# Creating attention fail output T2. This may then be used to remove attention fails.
attention_fail_t2 <- function() {
  attention_f <- t2[(t2$T2_lifesat_4 != 2) | (t2$T2_vita_5 != 3) | (t2$T2_entresat_4 != 6),]
  ID_vals <- data.frame(table(attention_f$T2_PID))
  Rows_fails <- attention_f$T2_PID %in% ID_vals[ID_vals$Freq > 0,1]
  Att_fails <- attention_f[Rows_fails,]
  return(Att_fails)
}


# Returning summary object for lavaan latent mediation model for models with variables in measurement part loading
# on three items each. With control variables = T1employable + T1employees + T1timebuiss + T1occ.

lavaan_sem <- function(IV, mediator, DV, controls = F) {
  if (controls == T) {
    controlinput <- paste('+ T1employable + T1employees + T1timebuiss + T1occ', '\n')
  } else {
    controlinput <- paste("", "\n")
  }

  measurement_c <- function() {
    measurement_part <- function(...) {
      params <- list(...)
      stopifnot(length(params)%%2==0)
      lefts = params[seq(1,length(params), by=2)]
      rights = params[seq(2,length(params), by=2)]
      rights <- Map(paste, rights, collapse="+")
      paste(paste0(lefts, " =~", rights), collapse="\n", "\n")
    }

    if(controls == F){
      meas <- measurement_part(paste0(IV), c(paste0(IV,"_1"), paste0(IV, "_2"), paste0(IV,"_3")),
                               paste0(mediator), c(paste0(mediator,"_1"), paste0(mediator, "_2"), paste0(mediator,"_3")),
                               paste0(DV), c(paste0(DV,"_1"), paste0(DV, "_2"), paste0(DV,"_3")))
    } else {
      meas <- measurement_part(paste0(IV), c(paste0(IV,"_1"), paste0(IV, "_2"), paste0(IV,"_3")),
                               paste0(mediator), c(paste0(mediator,"_1"), paste0(mediator, "_2"), paste0(mediator,"_3")),
                               paste0(DV), c(paste0(DV,"_1"), paste0(DV, "_2"), paste0(DV,"_3")),
                               "T1employable", c("T1employable_1", "T1employable_2", "T1employable_3"))

    }
    meas
  }
  model<- paste(measurement_c(),
                DV, "~ c*", IV, controlinput,
                mediator, "~ a*", IV, "\n",
                DV, "~ b*", mediator, "\n",
                "ab := a*b", '\n',
                "total := c + (a*b)")
  fit <- sem(model, data = obs_df)
  summary(fit, fit.measures=TRUE)
}


# center x and y variables in polynomial regression analysis on their grand mean
Center <- function(longer_comp) {

  # center x and y on grand mean across both variables
  grand.M <- mean( c(longer_comp$vitality, longer_comp$learn), na.rm = TRUE )
  longer_comp$x.c <- longer_comp$vitality-grand.M
  longer_comp$y.c <- longer_comp$learn-grand.M


  # compute higher order terms based on the centered predictors
  longer_comp$x2.c <- longer_comp$x.c^2
  longer_comp$xy.c <- longer_comp$x.c*longer_comp$y.c
  longer_comp$y2.c <- longer_comp$y.c^2

  # save means and compute their higher order terms
  longer_comp$x.mean  <- mean( longer_comp$vitality, na.rm = TRUE )
  longer_comp$y.mean  <- mean( longer_comp$learn, na.rm = TRUE )
  longer_comp$x2.mean <- longer_comp$x.mean^2
  longer_comp$xy.mean <- longer_comp$x.mean*longer_comp$y.mean
  longer_comp$y2.mean <- longer_comp$y.mean^2

  longer_comp

}

# Function to compute the intraclass correlation coefficient (ICC)
compute_icc <- function(lmer_object){
  var_dat <- lmer_object %>% VarCorr %>% as.data.frame
  icc <- var_dat$vcov[1]/(var_dat$vcov[1]+var_dat$vcov[2])
  return(icc)
}

# Create standardized coefficients
stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- round(fixef(object)*sdx/sdy, 2)
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc))
}


#  R function multilevel_alpha(), a wrapper for computing
#  the reliability indices discussed in
#  Lai, M. H. C. (2021). Composite reliability of multilevel data:
#      It’s about observed scores and construct meanings.
#      Psychological Methods, 26(1), 90–102.
#      https://doi.org/10.1037/met0000287
#  Copyright (C) 2021 Lai, Hok Chio (Mark)
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

multilevel_alpha <- function(data, id, nsim = 5000, conf_level = .95,
                             se = "robust.huber.white",
                             pa_args = list(fa = "pc"), ...) {
  if (!require(lavaan)) stop("The lavaan package needs to be installed.")
  if (!require(psych)) stop("The psych package needs to be installed.")
  nitem <- ncol(data)
  ynames <- paste0("y", seq_len(nitem))
  colnames(data) <- ynames
  data <- cbind(data, id = id)
  hmean_cs <- 1 / mean(1 / table(id))
  # Alpha
  # Generate syntax for saturated model
  sat_syntax <- (function(y) {
    if (length(y) <= 1) {
      return(NULL)
    }
    paste(
      c(paste(y[1], "~~", paste(y[-1], collapse = " + ")),
        Recall(y[-1])),
      collapse = "\n  "
    )
  })(ynames)
  msat <- paste0("level: 1\n  ", sat_syntax, "\nlevel: 2\n  ", sat_syntax)
  msat_fit <- cfa(msat, data = data, cluster = "id", se = se,
                  test = "none", h1 = FALSE, baseline = FALSE, ...)
  coef_msat <- coef(msat_fit, type = "user")
  vcov_msat <- vcov(msat_fit)
  vw <- names(coef_msat)[
    with(msat_fit@ParTable, which(op == "~~" & lhs == rhs & level == 1))]
  cvw <- names(coef_msat)[
    with(msat_fit@ParTable, which(op == "~~" & lhs != rhs & level == 1))]
  vb <- names(coef_msat)[
    with(msat_fit@ParTable, which(op == "~~" & lhs == rhs & level == 2))]
  cvb <- names(coef_msat)[
    with(msat_fit@ParTable, which(op == "~~" & lhs != rhs & level == 2))]
  Sw <- sum(coef_msat[vw], 2 * coef_msat[cvw])
  Sb <- sum(coef_msat[vb], 2 * coef_msat[cvb])
  alpha_const <- nitem / (nitem - 1)
  alphaw <- alpha_const * sum(2 * coef_msat[cvw]) / Sw
  alpha2l <- alpha_const * sum(2 * coef_msat[c(cvw, cvb)]) / (Sb + Sw)
  alphab <- alpha_const * sum(2 * coef_msat[cvb]) / (Sb + Sw / hmean_cs)
  sim_coef_msat <- MASS::mvrnorm(nsim,
                                 mu = coef_msat[c(vw, vb, cvw, cvb)],
                                 Sigma = vcov_msat[c(vw, vb, cvw, cvb),
                                                   c(vw, vb, cvw, cvb)])
  sim_Sw <- rowSums(cbind(sim_coef_msat[ , vw], 2 * sim_coef_msat[ , cvw]))
  sim_Sb <- rowSums(cbind(sim_coef_msat[ , vb], 2 * sim_coef_msat[ , cvb]))
  sim_alphaw <- alpha_const * rowSums(2 * sim_coef_msat[ , cvw]) / sim_Sw
  sim_alpha2l <- alpha_const * rowSums(2 * sim_coef_msat[ , c(cvw, cvb)]) /
    (sim_Sb + sim_Sw)
  sim_alphab <- alpha_const * rowSums(2 * sim_coef_msat[ , cvb]) /
    (sim_Sb + sim_Sw / hmean_cs)
  sim_alpha_cis <- lapply(list(alpha2l = sim_alpha2l,
                               alphab = sim_alphab,
                               alphaw = sim_alphaw),
                          quantile,
                          probs = .5 + c(- conf_level, conf_level) / 2)
  # Omega
  loading_labels <- paste0("l", seq_len(nitem))
  g_syntax <- paste(loading_labels, "*", ynames,
                    collapse = " + ")
  mcfa <- paste0("level: 1\n  fw =~ ",
                 g_syntax,
                 "\nlevel: 2\n  fb =~ ",
                 g_syntax)
  mcfa_fit <- cfa(mcfa, data = data, cluster = "id", se = se,
                  test = "none", h1 = TRUE, baseline = FALSE, ...)
  mcfa_pt <- partable(mcfa_fit)
  coef_mcfa <- coef(mcfa_fit)
  vcov_mcfa <- vcov(mcfa_fit)
  coef_mcfa[loading_labels]
  ld <- names(coef_mcfa)[with(mcfa_pt, free[which(op == "=~" &
                                                    level == 1)])]
  evw <- names(coef_mcfa)[with(mcfa_pt,
                               free[which(op == "~~" &
                                            lhs == rhs &
                                            lhs != "fw" &
                                            level == 1)])]
  fvw <- names(coef_mcfa)[with(mcfa_pt, free[which(op == "~~" &
                                                     lhs == "fw")])]
  evb <- names(coef_mcfa)[with(mcfa_pt,
                               free[which(op == "~~" &
                                            lhs == rhs & lhs != "fb" &
                                            level == 2)])]
  fvb <- names(coef_mcfa)[with(mcfa_pt, free[which(op == "~~" &
                                                     lhs == "fb")])]
  sumldsq <- sum(1, coef_mcfa[ld])^2
  sumevw <- sum(coef_mcfa[evw])
  sumevb <- sum(coef_mcfa[evb])
  omegaw <- sumldsq * coef_mcfa[[fvw]] /
    (sumldsq * coef_mcfa[[fvw]] + sumevw)
  omega2l <- sum(sumldsq * coef_mcfa[c(fvw, fvb)]) /
    sum(sumldsq * coef_mcfa[c(fvw, fvb)], sumevw, sumevb)
  omegab <- sumldsq * coef_mcfa[[fvb]] /
    (sumldsq * (coef_mcfa[[fvb]] + coef_mcfa[[fvw]] / hmean_cs) +
       sumevb + sumevw / hmean_cs)
  sim_coef_mcfa <- MASS::mvrnorm(nsim,
                                 mu = coef_mcfa[c(ld, fvw, fvb, evw, evb)],
                                 Sigma = vcov_mcfa[c(ld, fvw, fvb, evw, evb),
                                                   c(ld, fvw, fvb, evw, evb)])
  sim_sumldsq <- (1 + rowSums(sim_coef_mcfa[ , ld]))^2
  sim_sumevw <- rowSums(sim_coef_mcfa[ , evw])
  sim_sumevb <- rowSums(sim_coef_mcfa[ , evb])
  sim_omegaw <- sim_sumldsq * sim_coef_mcfa[ , fvw] /
    (sim_sumldsq * sim_coef_mcfa[ , fvw] + sim_sumevw)
  sim_omega2l <- rowSums(sim_sumldsq * sim_coef_mcfa[ , c(fvw, fvb)]) /
    (rowSums(sim_sumldsq * sim_coef_mcfa[ , c(fvw, fvb)]) +
       sim_sumevw + sim_sumevb)
  sim_omegab <- sim_sumldsq * sim_coef_mcfa[ , fvb] /
    (sim_sumldsq * (sim_coef_mcfa[ , fvb] +
                      sim_coef_mcfa[ , fvw] / hmean_cs) +
       sim_sumevb + sim_sumevw / hmean_cs)
  sim_omega_cis <- lapply(list(omega2l = sim_omega2l,
                               omegab = sim_omegab,
                               omegaw = sim_omegaw),
                          quantile,
                          probs = .5 + c(- conf_level, conf_level) / 2)
  # resid_corb <- resid(mcfa_fit, type = "cor")$id$cov
  # diag(resid_corb) <- 1
  # psych::KMO(resid_corb)$MSA
  # psych::fa.parallel(resid_corb, fm = "pa", fa = "fa",
  #                    n.obs = lavTech(msat_fit, "nclusters")[[1]],
  #                    n.iter = 30 * nitem,
  #                    plot = FALSE)$nfact
  # Dimensionality
  corw <- lavTech(msat_fit, what = "cor.ov")$within
  corb <- lavTech(msat_fit, what = "cor.ov")$id
  paw <- do.call(fa.parallel,
                 args = c(list(x = corw,
                               n.obs = nobs(msat_fit) -
                                 lavTech(msat_fit, "nclusters")[[1]],
                               n.iter = 30 * nitem,
                               plot = FALSE), pa_args))
  pab <- do.call(fa.parallel,
                 args = c(list(x = corw,
                               n.obs = lavTech(msat_fit, "nclusters")[[1]],
                               n.iter = 30 * nitem,
                               plot = FALSE), pa_args))
  if (pa_args$fa == "pc") {
    ncompw <- paw$ncomp
    ncompb <- pab$ncomp
  } else if (pa_args$fa == "fa") {
    ncompw <- paw$nfact
    ncompb <- pab$nfact
  }
  list(alpha = c(alpha2l = alpha2l, alphab = alphab, alphaw = alphaw),
       alpha_ci = do.call(rbind, sim_alpha_cis),
       omega = c(omega2l = omega2l, omegab = omegab, omegaw = omegaw),
       omega_ci = do.call(rbind, sim_omega_cis),
       ncomp = c(within = ncompw, between = ncompb))
}

