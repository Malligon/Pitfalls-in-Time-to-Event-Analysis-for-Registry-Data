# Effect of censoring ----

library(survival)
library(stringr)
library(ggsci)

set.seed(37)

N = 1000

shape = 2
scale = 30

T <- rweibull(N, shape, scale) # true times
C <- runif(N, 0, 85) # censored times

status <- T <= C # status is 0 if censored, 1 either, ie if censoring happens before event
Tobs <- pmin(T, C) # observed time, ie min between time of censoring and time of event

xaxis<-seq(0,50,by=0.01)

DB <- data.frame(Tobs = c(Tobs, Tobs[status == 1]), 
                 status = as.numeric(c(status, status[status == 1])),
                 Cens = c(rep(1, length(status)), rep(2, length(status[status == 1]))))

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure1.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)
plot_km_trunc(fit = survfit(Surv(Tobs, status) ~ Cens, data = DB),
              data = DB,
              Trunc = FALSE,
              col = c("blue", "red"),
              strata.names = c("KM", "Naive"),
              xlim = c(-1, 50),
              ylim = c(-0.05, 1.05),
              legend = FALSE,
              lwd = 0.5,
              legend.lwd = 0.5)
lines(xaxis, exp(-(xaxis/scale)^shape), type = "l", lwd = 0.5) # true survival
legend("topright", legend = c("True survival", "Kaplan-Meier estimator", "Naive approach"), col = c("black", "blue", "red"), lty = 1)
dev.off()

# Effect of truncature ----

set.seed(37)

N = 1000

shape = 2
scale = 30

T <- rweibull(N, shape, scale) # true times
C <- runif(N, 0, 85) # censored times

Ttrunc <- runif(N, 0, 50) # Time of truncature

status <- T <= C
Tobs <- pmin(T, C)

Tobs2 <- Tobs[Tobs >= Ttrunc] # only times of event after truncature could be observed
status2 <- status[Tobs >= Ttrunc] # idem for the status of event
Ttrunc2 <- Ttrunc[Tobs >= Ttrunc] # idem for times of truncature

KM <- survfit(Surv(Ttrunc2, Tobs2, status2) ~ 1)
KM_false <- survfit(Surv(Tobs2, status2) ~ 1)

DB <- data.frame(Tobs = c(Tobs2, Tobs2), 
                 status = c(status2, status2),
                 start = c(Ttrunc2, rep(0, length(status2))),
                 trunc = c(rep(1, length(status2)), rep(2, length(status2))))

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure2.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)
plot_km_trunc(survfit(Surv(start, Tobs, status) ~ trunc, data = DB),
              DB,
              Trunc = TRUE,
              col = c("blue", "red"),
              strata.names = c("KM", "Naive"),
              xlim = c(-1, 50),
              ylim = c(-0.05, 1.05),
              legend = FALSE,
              lwd = 0.5,
              legend.lwd = 0.5)
lines(xaxis, exp(-(xaxis/scale)^shape), type = "l",
      lwd = 0.5 ) # true survival
legend("topright", legend = c("True survival", "KM taking into account truncation", "KM ignoring truncation"), col = c("black", "blue", "red"), lty = 1)
dev.off()

# Effect of competing risk ----

set.seed(37)

N = 1000

shape = 2
scale_1 = 30
T1 <- rweibull(N, shape, scale_1) 

shape = 2 
scale_2 = 40
T2 <- rweibull(N, shape, scale_2) # competing risk

T <- pmin(T1, T2)
C <- runif(N, 0, 85)

status <- as.numeric(T2 < T1) + 1 
status[C<T] <- 0 # 0 for censoring, 1 for event, 2 for competing risk

time_obs <- pmin(C, T)

true_surv <- scale_1^shape/(scale_1^shape + scale_2^shape) + exp(-xaxis^shape*(scale_1^shape + scale_2^shape)/(scale_1^shape * scale_2^shape))*(1-scale_1^shape/(scale_1^shape + scale_2^shape)) # True survival for outcome 1

true_surv_2 <- scale_2^shape/(scale_1^shape + scale_2^shape) + exp(-xaxis^shape*(scale_1^shape + scale_2^shape)/(scale_1^shape * scale_2^shape))*(1-scale_2^shape/(scale_1^shape + scale_2^shape)) # True survival for outcome 2

DB <- data.frame(Tobs = time_obs, 
                 status = status)

DB$status_CR <- as.factor(DB$status)

KM_false <- survfit(Surv(time_obs, status==1)~1) # Survival for outcome 1 with censoring on outcome 2
KM_false_2 <- survfit(Surv(time_obs, status==2)~1) # Survival for outcome 2 with censoring on outcome 1

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure3.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)
plot_km_trunc_CR(fit = survfit(Surv(Tobs, status_CR) ~ 1, data = DB),
                 data = DB,
                 Trunc = FALSE,
                 col = c("blue"),
                 xlim = c(-1, 50),
                 ylim = c(-0.05, 1.05),
                 legend = FALSE,
                 risk.table = TRUE,
                 lwd = 0.5)
lines(xaxis, 1-true_surv, col = "black",
      lwd = 0.5) # True survival for outcome 1
lines(xaxis, 1-true_surv_2, col = "black", lty = 2,
      lwd = 0.5) # True survival for outcome 2
lines(KM_false$time, 1-KM_false$surv, col = "red", type = "s",
      lwd = 0.5) # Survival for outcome 1 with censored data
lines(KM_false_2$time, 1-KM_false_2$surv, col = "red", lty = 2, type = "s",
      lwd = 0.5) # survival for outcome 2 with censored data
lines(KM_false$time, 1-KM_false$upper, col = "red", lty = 3, type = "s",
      lwd = 0.5)
lines(KM_false$time, 1-KM_false$lower, col = "red", lty = 3, type = "s",
      lwd = 0.5)
lines(KM_false_2$time, 1-KM_false_2$upper, col = "red", lty = 3, type = "s",
      lwd = 0.5)
lines(KM_false_2$time, 1-KM_false_2$lower, col = "red", lty = 3, type = "s",
      lwd = 0.5)

legend("topleft", legend = c("True survival", "Competing risk method", "KM (naive method)", "", "Event", "Competing risk"), col = c("black", "blue", "red", "white", "black", "black"), lty = c(1, 1, 1, 1, 1, 2), lwd = 0.5)
dev.off()

# Recurrent event ----

# Recurrent event: simulation ----

set.seed(37)

N = 1000
shape = 2
scale = 30

shape_event = 2.5
scale_event = 20

randef = rep(1, N) # when simulating a random effect, to cancel, randef = 1

C <- runif(N, 0, 85)
T <- rweibull(N, shape, scale)

status <- T <= C
Tobs <- pmin(T, C)

E1 <- scale_event*(rexp(N, 1)/randef)^(1/shape_event)
j = 1
nbevent <- rep(NA, N)
recdata <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(recdata) <- c("id", "start", "stop", "status", "terminal")
for (i in 1:N)
{
  stop = E1[i]
  start = 0
  while(stop < Tobs[i])
  {
    recdata[j, ] <- c(i, start, stop, 1, 0)
    start = stop
    stop = scale_event*((stop/scale_event)^(shape_event)+rexp(1, 1)/randef[i])^(1/shape_event)
    j = j+1
  }
  nbevent[i] <- j-1
  recdata[j,] <- c(i, start, Tobs[i], 0, status[i])
  j = j+1
}

# Recurrent event: estimation and plot ----

rec_est <- rec_event(id = "id", start="start", stop="stop", status="status", terminal="terminal", recdata, recdata)

f=function(x){(x/scale_event)^{shape_event-1}*exp(-(x/scale)^shape)*shape_event/scale_event}

# plot

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure4.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)

par(mfrow = c(1, 2), mar = c(6 + 1, 4.5, 4, 3))

plot.rec_event(x = rec_est,
               conf.int = TRUE,
               risk.table = TRUE,
               main = "", 
               xlab="Time", 
               ylab="Mean number of recurrent events",
               col = "blue",
               xlim = c(0, 50),
               lwd = 0.5)
lines(seq(0, 50, by=0.01), sapply(xaxis,function(y){integrate(f,lower=0,upper=y)$value}), type = "l", col = "black", lty = 1,
      lwd = 0.5)
lines(rec_est$expecNP$time, cumsum(rec_est$expecNP$n.event)/N, type="l", col="red", lty=1,
      lwd = 0.5)

legend("topleft",c("Truth","Correct estimator", "Naive estimator"),lwd=0.5,lty=c(1,1)
       ,col=c("black", "blue", "red"),inset=0.05)

plot_km_trunc(survfit(Surv(start, stop, terminal) ~ 1, data = recdata), 
              data = recdata, 
              xlim = c(0, 50), 
              col = "blue", 
              legend = FALSE, 
              Trunc = TRUE,
              lwd = 0.5)
lines(xaxis, exp(-(xaxis/scale)^shape), type = "l", col="black",
      lwd = 0.5)

legend("bottomleft",c("Truth","Correct estimator"),lwd=0.5,lty=c(1,1)
       ,col=c("black", "blue"),inset=0.05)

par(mfrow = c(1, 1))
dev.off()



# CEREDIH: importation ----

DB <- read.csv2("C:/Users/Mika/OneDrive/Desktop/TADbis/TAd epistat/CEREDIH_randomised.csv", dec=",", na.strings=c("NA",""," ",".","--", "N/A"))

# KM start diag ----

DB_cvid <- DB[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID"), ]

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure6A.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)

plot_km_trunc(fit = survfit(Surv(time, status) ~ age_clinical_diagnosis_class, data = DB_cvid),
              data = DB_cvid,
              Trunc = FALSE,
              conf.int = FALSE,
              col = pal_jco(palette = c("default"), alpha = 1)(5),
              strata.names = c("0-4", "5-9", "10-19", "20-39", "40+"),
              xlim = c(-1, 50),
              ylim = c(-0.05, 1.05),
              legend = TRUE,
              legend.x = "bottomleft",
              lwd = 0.5)
dev.off()

# KM start birth ----

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure6B.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)
plot_km_trunc(fit = survfit(Surv(age_clinical_diagnosis, last_news_age, status) ~ cat_dih1ter, data = DB),
              data = DB,
              Trunc = TRUE,
              conf.int = FALSE,
              col = pal_jco(palette = c("default"), alpha = 1)(6),
              strata.names = c("CVID", "non-CVID", "SCID", "CID", "OtherT", "Innate"),
              xlim = c(-1, 100),
              ylim = c(-0.05, 1.05),
              legend = TRUE,
              legend.x = "topright",
              legend.legend = c("CVID", "non-CVID", "SCID", "CID", "Other Tcell def.", "Innate"),
              risk.table = TRUE,
              lwd = 0.5)
dev.off()

# Naive vs true estimate truncation cvid ----

DB_comp <- rbind(DB[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID"), ], 
                 DB[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID"), ])

DB_comp$start <- c(DB$age_clinical_diagnosis[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID")], 
                   rep(0, nrow(DB[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID"), ])))

DB_comp$flag <- c(rep(1, nrow(DB[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID"), ])), 
                  rep(2, nrow(DB[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID"), ])))

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure7.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)

plot_km_trunc(Surv(start, last_news_age, status) ~ flag,
              data = DB_comp,
              Trunc = TRUE,
              conf.int = TRUE,
              col = c("blue", "red"),
              strata.names = c("KM/trunc+", "KM/trunc-"),
              xlim = c(-1, 100),
              ylim = c(-0.05, 1.05),
              legend = TRUE,
              legend.x = "bottomleft",
              risk.table = TRUE,
              lwd = 0.5)

dev.off()

# Competing risks ----

DB$time_cancer_cmprsk <- pmin(DB$first_malignancy_age, DB$last_news_age_with_censored_hsct, na.rm = TRUE)
DB$status_cancer_cmprsk <- ifelse(DB$last_news_age_with_censored_hsct < DB$last_news_age | DB$Survival_STATUS == "DEAD", "Competing risks", 0)
DB$status_cancer_cmprsk <- ifelse(!is.na(DB$first_malignancy_age) & DB$first_malignancy_age <= DB$last_news_age_with_censored_hsct, "Cancer", DB$status_cancer_cmprsk)

KM_cancer_cvid <- survfit(Surv(time_cancer_cmprsk, status_cancer_cmprsk == "Cancer") ~ 1, data = DB[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID"), ])
KM_cancer_noncvid <- survfit(Surv(time_cancer_cmprsk, status_cancer_cmprsk == "Cancer") ~ 1, data = DB[which(DB$cat_dih1ter == "2.Adaptive Immunity B-cell deficiencies non-CVID"), ])
KM_cancer_scid <- survfit(Surv(time_cancer_cmprsk, status_cancer_cmprsk == "Cancer") ~ 1, data = DB[which(DB$cat_dih1ter == "3.Adaptive Immunity T-cell deficiencies SCID"), ])
KM_cancer_cid <- survfit(Surv(time_cancer_cmprsk, status_cancer_cmprsk == "Cancer") ~ 1, data = DB[which(DB$cat_dih1ter == "4.Adaptive Immunity T-cell deficiencies CID"), ])
KM_cancer_otherT <- survfit(Surv(time_cancer_cmprsk, status_cancer_cmprsk == "Cancer") ~ 1, data = DB[which(DB$cat_dih1ter == "5.Adaptive Immunity Other T-cell deficiencies"), ])
KM_cancer_innate <- survfit(Surv(time_cancer_cmprsk, status_cancer_cmprsk == "Cancer") ~ 1, data = DB[which(DB$cat_dih1ter == "6.Innate Immunity Deficiencies"), ])

KM_surv_cvid <- survfit(Surv(age_clinical_diagnosis, time_cancer_cmprsk, status_cancer_cmprsk == "Competing risks")~1, data = DB[which(DB$cat_dih1ter == "1.Adaptive Immunity B-cell deficiencies CVID"), ])
KM_surv_noncvid <- survfit(Surv(age_clinical_diagnosis, time_cancer_cmprsk, status_cancer_cmprsk == "Competing risks")~1, data = DB[which(DB$cat_dih1ter == "2.Adaptive Immunity B-cell deficiencies non-CVID"), ])
KM_surv_scid <- survfit(Surv(age_clinical_diagnosis, time_cancer_cmprsk, status_cancer_cmprsk == "Competing risks")~1, data = DB[which(DB$cat_dih1ter == "3.Adaptive Immunity T-cell deficiencies SCID"), ])
KM_surv_cid <- survfit(Surv(age_clinical_diagnosis, time_cancer_cmprsk, status_cancer_cmprsk == "Competing risks")~1, data = DB[which(DB$cat_dih1ter == "4.Adaptive Immunity T-cell deficiencies CID"), ])
KM_surv_otherT <- survfit(Surv(age_clinical_diagnosis, time_cancer_cmprsk, status_cancer_cmprsk == "Competing risks")~1, data = DB[which(DB$cat_dih1ter == "5.Adaptive Immunity Other T-cell deficiencies"), ])
KM_surv_innate <- survfit(Surv(age_clinical_diagnosis, time_cancer_cmprsk, status_cancer_cmprsk == "Competing risks")~1, data = DB[which(DB$cat_dih1ter == "6.Innate Immunity Deficiencies"), ])

KM_cancer_stepfun_cvid <- stepfun(KM_cancer_cvid$time, c(1, KM_cancer_cvid$surv))
KM_cancer_stepfun_noncvid <- stepfun(KM_cancer_noncvid$time, c(1, KM_cancer_noncvid$surv))
KM_cancer_stepfun_scid <- stepfun(KM_cancer_scid$time, c(1, KM_cancer_scid$surv))
KM_cancer_stepfun_cid <- stepfun(KM_cancer_cid$time, c(1, KM_cancer_cid$surv))
KM_cancer_stepfun_otherT <- stepfun(KM_cancer_otherT$time, c(1, KM_cancer_otherT$surv))
KM_cancer_stepfun_innate <- stepfun(KM_cancer_innate$time, c(1, KM_cancer_innate$surv))

KM_surv_stepfun_cvid <- stepfun(KM_surv_cvid$time, c(1, KM_surv_cvid$surv))
KM_surv_stepfun_noncvid <- stepfun(KM_surv_noncvid$time, c(1, KM_surv_noncvid$surv))
KM_surv_stepfun_scid <- stepfun(KM_surv_scid$time, c(1, KM_surv_scid$surv))
KM_surv_stepfun_cid <- stepfun(KM_surv_cid$time, c(1, KM_surv_cid$surv))
KM_surv_stepfun_otherT <- stepfun(KM_surv_otherT$time, c(1, KM_surv_otherT$surv))
KM_surv_stepfun_innate <- stepfun(KM_surv_innate$time, c(1, KM_surv_innate$surv))

# Warnings are all normal, indeed some patients are at risk to develop a cancer due to the retrospective aspect of this data, while they are not at risk to die. A clean way would be to to create a new database including only patients who have their time_cancer_cmprsk above age_clinical_diagnosis. Results would be te same

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure8A.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)

par(mfrow = c(1, 1), mar = c(6 + 6, 4.5, 4, 3))

plot(cbind(KM_cancer_cvid$time,
           cumsum(c(1, (KM_cancer_stepfun_cvid(KM_cancer_cvid$time)*KM_surv_stepfun_cvid(KM_cancer_cvid$time))[1:(length(KM_cancer_cvid$time)-1)])* KM_cancer_cvid$n.event/KM_cancer_cvid$n.risk)),
     col = pal_jco(palette = c("default"), alpha = 1)(6)[1],
     main = "",
     xlab = list("Time", cex = 1.2),
     ylab = list("Cumulative incidence function", cex = 1.2),
     type = "s",
     ylim = c(-0.05, 1.05),
     xlim = c(-3, 100),
     lwd = 0.5)
lines(cbind(KM_cancer_noncvid$time,
            cumsum(c(1, (KM_cancer_stepfun_noncvid(KM_cancer_noncvid$time)*KM_surv_stepfun_noncvid(KM_cancer_noncvid$time))[1:(length(KM_cancer_noncvid$time)-1)])* KM_cancer_noncvid$n.event/KM_cancer_noncvid$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[2], lwd = 0.5, type = "s")
lines(cbind(KM_cancer_scid$time,
            cumsum(c(1, (KM_cancer_stepfun_scid(KM_cancer_scid$time)*KM_surv_stepfun_scid(KM_cancer_scid$time))[1:(length(KM_cancer_scid$time)-1)])* KM_cancer_scid$n.event/KM_cancer_scid$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[3], lwd = 0.5, type = "s")
lines(cbind(KM_cancer_cid$time,
            cumsum(c(1, (KM_cancer_stepfun_cid(KM_cancer_cid$time)*KM_surv_stepfun_cid(KM_cancer_cid$time))[1:(length(KM_cancer_cid$time)-1)])* KM_cancer_cid$n.event/KM_cancer_cid$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[4], lwd = 0.5, type = "s")
lines(cbind(KM_cancer_otherT$time,
            cumsum(c(1, (KM_cancer_stepfun_otherT(KM_cancer_otherT$time)*KM_surv_stepfun_otherT(KM_cancer_otherT$time))[1:(length(KM_cancer_otherT$time)-1)])* KM_cancer_otherT$n.event/KM_cancer_otherT$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[5], lwd = 0.5, type = "s")
lines(cbind(KM_cancer_innate$time,
            cumsum(c(1, (KM_cancer_stepfun_innate(KM_cancer_innate$time)*KM_surv_stepfun_innate(KM_cancer_innate$time))[1:(length(KM_cancer_innate$time)-1)])* KM_cancer_innate$n.event/KM_cancer_innate$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[6], lwd = 0.5, type = "s")
legend("topleft", legend = c("CVID", "non-CVID", "SCID", "CID", "Other Tcell def.", "Innate"), col = pal_jco(palette = c("default"), alpha = 1)(6), lty = 1, lwd = 0.5)

name_at_risk <- "Risk table"
times.print <- axis(1, labels = FALSE, tick = FALSE)
at.risk <- matrix(c(summary(KM_cancer_cvid, times.print)$n.risk, rep(0, length(times.print)-length(summary(KM_cancer_cvid, times.print)$n.risk)),
                    summary(KM_cancer_noncvid, times.print)$n.risk, rep(0, length(times.print)-length(summary(KM_cancer_noncvid, times.print)$n.risk)),
                    summary(KM_cancer_scid, times.print)$n.risk, rep(0, length(times.print)-length(summary(KM_cancer_scid, times.print)$n.risk)),
                    summary(KM_cancer_cid, times.print)$n.risk, rep(0, length(times.print)-length(summary(KM_cancer_cid, times.print)$n.risk)),
                    summary(KM_cancer_otherT, times.print)$n.risk, rep(0, length(times.print)-length(summary(KM_cancer_otherT, times.print)$n.risk)),
                    summary(KM_cancer_innate, times.print)$n.risk, rep(0, length(times.print)-length(summary(KM_cancer_innate, times.print)$n.risk))), nrow = 6, byrow = TRUE)

mtext(side = 1, at = times.print[2]/2, 
      line = 4, name_at_risk, cex = par("cex"), adj = 1)

for (i in 1:nrow(at.risk))
{
  mtext(side = 1, at = -0.3 * (times.print[2] - times.print[1]) + 
          times.print[1], line = i + 4, c("CVID", "nonCVID", "SCID", "CID", "OtherT", "Innate")[i], cex = par("cex"), 
        adj = 1)
  mtext(side = 1, at = times.print, line = i + 4, at.risk[i, ], cex = par("cex"))
}

dev.off()

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure8B.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)

plot(cbind(KM_surv_cvid$time, 
           cumsum(c(1, (KM_cancer_stepfun_cvid(KM_surv_cvid$time)*KM_surv_stepfun_cvid(KM_surv_cvid$time))[1:(length(KM_surv_cvid$time)-1)])* KM_surv_cvid$n.event/KM_surv_cvid$n.risk)),
     col = pal_jco(palette = c("default"), alpha = 1)(6)[1],
     main = "",
     xlab = list("Time", cex = 1.2),
     ylab = list("Cumulative incidence function", cex = 1.2),
     type = "s",
     lwd = 0.5,
     ylim = c(-0.05, 1.05),
     xlim = c(-3, 100))
lines(cbind(KM_surv_noncvid$time, 
            cumsum(c(1, (KM_cancer_stepfun_noncvid(KM_surv_noncvid$time)*KM_surv_stepfun_noncvid(KM_surv_noncvid$time))[1:(length(KM_surv_noncvid$time)-1)])* KM_surv_noncvid$n.event/KM_surv_noncvid$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[2], lwd = 0.5, type = "s")
lines(cbind(KM_surv_scid$time, 
            cumsum(c(1, (KM_cancer_stepfun_scid(KM_surv_scid$time)*KM_surv_stepfun_scid(KM_surv_scid$time))[1:(length(KM_surv_scid$time)-1)])* KM_surv_scid$n.event/KM_surv_scid$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[3], lwd = 0.5, type = "s")
lines(cbind(KM_surv_cid$time, 
            cumsum(c(1, (KM_cancer_stepfun_cid(KM_surv_cid$time)*KM_surv_stepfun_cid(KM_surv_cid$time))[1:(length(KM_surv_cid$time)-1)])* KM_surv_cid$n.event/KM_surv_cid$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[4], lwd = 0.5, type = "s")
lines(cbind(KM_surv_otherT$time, 
            cumsum(c(1, (KM_cancer_stepfun_otherT(KM_surv_otherT$time)*KM_surv_stepfun_otherT(KM_surv_otherT$time))[1:(length(KM_surv_otherT$time)-1)])* KM_surv_otherT$n.event/KM_surv_otherT$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[5], lwd = 0.5, type = "s")
lines(cbind(KM_surv_innate$time, 
            cumsum(c(1, (KM_cancer_stepfun_innate(KM_surv_innate$time)*KM_surv_stepfun_innate(KM_surv_innate$time))[1:(length(KM_surv_innate$time)-1)])* KM_surv_innate$n.event/KM_surv_innate$n.risk)),
      col = pal_jco(palette = c("default"), alpha = 1)(6)[6], lwd = 0.5, type = "s")
legend("bottomright", legend = c("CVID", "non-CVID", "SCID", "CID", "Other Tcell def.", "Innate"), col = pal_jco(palette = c("default"), alpha = 1)(6), lty = 1, lwd = 0.5)

name_at_risk <- "Risk table"
times.print <- axis(1, labels = FALSE, tick = FALSE)

at.risk <- c()
for (k in levels(as.factor(as.character(DB$cat_dih1ter))))
{
  for (i in times.print)
  {
    at.risk <- c(at.risk, nrow(DB[which(DB$cat_dih1ter == k
                                        & DB$age_clinical_diagnosis <= i 
                                        & DB$time_cancer_cmprsk >= i), ]))
  }
}
at.risk <- matrix(at.risk, nrow = 6, byrow = TRUE)

# mtext(side = 1, at = times.print[2]/2, 
# line = 4, name_at_risk, cex = par("cex"), adj = 1)

for (i in 1:nrow(at.risk))
{
  mtext(side = 1, at = -0.3 * (times.print[2] - times.print[1]) + 
          times.print[1], line = i + 4, c("CVID", "nonCVID", "SCID", "CID", "OtherT", "Innate")[i], cex = par("cex"), 
        adj = 1)
  mtext(side = 1, at = times.print, line = i + 4, at.risk[i, ], cex = par("cex"))
}

dev.off()

# Counting process ----

CP <- as.data.frame(matrix(ncol = 7))

names(CP) <- c("id", "start", "stop", "event", "terminal", "cat_dih", "age_diag")

for (i in 1:nrow(DB))
{
  times <- c()
  id <- c()
  start <- c()
  stop <- c()
  event <- c()
  terminal <- c()
  cat_dih <- c()
  age_diag <- c()
  
  for (j in 1:length(CM_type))
  {
    if (!is.na(DB[i, CM_age[j]]) & DB[i, CM_type[j]] %in% c("MALIGNANCY", "AUTOIMMUNITY"))
    {
      times <- c(times, DB[i, CM_age[j]])
    }
  }
  times[times <= 0] <- 1/365
  if (length(times)>1)
  {
    times <- times[order(times)]
    
    for (j in 1:(length(times)-1))
    {
      if(times[j] > times[j+1] - 0.0001)
      {
        times[j+1] <- times[j] + 1/365
      }
    }
  }
  
  times <- times[times < DB$last_news_age_with_censored_hsct[i]]
  
  id <- rep(DB$PATIENT_ID[i], length(times)+1)
  start <- c(0, times)
  stop <- c(times, DB$last_news_age_with_censored_hsct[i])
  event <- c(rep(1, length(times)), 0)
  terminal <- c(rep(0, length(times)), DB$last_news_age_with_censored_hsct[i] < (DB$last_news_age[i]-0.0001) | DB$status_cens[i] == 1)
  cat_dih <- c(rep(as.character(DB$cat_dih1ter[i]), length(start)))
  age_diag <- c(rep(DB$age_clinical_diagnosis[i], length(start)))
  
  if (length(times)>1)
  {
    if (times[length(times)-1] > (times[length(times)]-0.0001))
    {
      times[length(times)] <- times[length(times)-1] + 1/365
    }
  }
  
  x <- cbind(id, start, stop, event, terminal, cat_dih, age_diag)
  
  CP <- rbind(CP, x)
}

CP <- CP[-1, ]

CP$stop <- as.numeric(CP$stop)
CP$start <- as.numeric(CP$start)
CP$event <- as.numeric(CP$event)
CP$terminal <- as.numeric(CP$terminal)
CP$age_diag <- as.numeric(CP$age_diag)

CP <- CP[which(!CP$id %in% CP$id[which(CP$start >= CP$stop)]), ]

CP2 <- CP[which(CP$stop > CP$age_diag), ]
CP2$start[which(CP2$start < CP2$age_diag)] <- CP2$age_diag[which(CP2$start < CP2$age_diag)]

CP <- CP[which(CP$id %in% CP2$id), ]

# View(CP[which(CP$id %in% 10039), ])

# Recurrent event ----

test <- rec_event(id = "id", 
                  start="start", 
                  stop="stop", 
                  status="event", 
                  terminal="terminal", 
                  recdata=CP[which(CP$cat_dih == "4.Adaptive Immunity T-cell deficiencies CID"), ],
                  recdata2=CP2[which(CP2$cat_dih == "4.Adaptive Immunity T-cell deficiencies CID"), ])

recdata2 <- CP2[which(CP2$cat_dih == "4.Adaptive Immunity T-cell deficiencies CID"), ]

png("O:/CEREDIH/Registre/Mika/CEREDIH/Survie/EPISTAT_figure9.png", width = 3.25, height = 3.25, units = "in", res = 1200, pointsize = 4)

par(mfrow = c(1, 2), mar = c(6 + 1, 4.5, 4, 3))

plot.rec_event(x = test,
               conf.int = TRUE,
               risk.table = TRUE,
               main = "", 
               xlab="Time", 
               ylab="Mean number of recurrent events",
               col = pal_jco(palette = c("default"), alpha = 1)(6)[4],
               xlim = c(0, 100),
               lwd = 0.5)

plot_km_trunc(survfit(Surv(start, stop, terminal) ~ 1, data = recdata2), 
              data = recdata2, 
              main = "",
              xlim = c(0, 100), 
              col = pal_jco(palette = c("default"), alpha = 1)(6)[4], 
              legend = FALSE, 
              Trunc = TRUE,
              risk.table = TRUE,
              conf.int = TRUE,
              lwd = 0.5)
dev.off()

par(mfrow = c(1, 1))
