# plot_km_trunc ----

# Documentation
# fit: must be either a survfit or a survival formula 
# survfit object should be displayed exactly that way : survfit(survival fomula, data = database) // ex : survfit(Surv(start, stop, status) ~ var+var2+var3, data = DB)
# survival formula should be displayed exactly way : Surv() ~ covariates // ex : Surv(stop, status) ~ 1, Surv(start, stop, status) ~ var+var2+var3
# data: dataframe, if fit is a survfit object, must be the same dataframe than in fit
# Trunc: logical value, must be set to TRUE if fit implies a start stop notation, FALSE otherwise
# risk.table: logical value informating if the risk table should be written under plot or no
# all risk.table options are options for risk table settings
# margin : number of row in the final risk table
# times.print: chose the time you want to display n at risk at
# at.name: adjustement of the name of each rows, see mtext for more information
# at.number: adjustement of the all the n at risk, see mtext for more information
# risk.table.cex.name: text size for the name of each rows, see mtext for more information
# risk.table.cex.number: text size for all the n at risk, see mtext for more information
# risk.table.adj.name: alignement for the name of each rows, see mtext for more information
# risk.table.adj.number: alignement for the all the n at risk, see mtext for more information
# strata.names: string vector, contains names of each strata of fit, must have the same length than the number of strata (~ 1 is consider as one unique strata)
# conf.int: Logical value informating weither or not the confidence intervals should be displayed
# main: title for the plot
# xlab: label for the horizontal axis
# ylab: label for the vertical axis
# xlim: numeric vector containing the min and max values for horizontal axis
# ylim: numeric vector containing the min and max values for vertical axis
# col: vector of characters informating the colours of strata, must have the same length than the number of strata (~ 1 is consider as one unique strata)
# legend: set to TRUE by default to display the legend, FALSE otherwise
# all legend.XXXX options are options for legend

plot_km_trunc <- function(fit,
                          data,
                          Trunc = FALSE,
                          strata.names = NULL,
                          conf.int = TRUE,
                          main = "",
                          xlab = "Time",
                          ylab = "Survival",
                          xlim = NULL,
                          ylim = c(-0.05, 1.05),
                          col = NULL,
                          lty = 1,
                          lwd = 1,
                          plot.margin = TRUE,
                          risk.table = TRUE,
                          risk.table.times.print = NULL,
                          risk.table.at.name = NULL,
                          risk.table.at.number = NULL,
                          risk.table.cex.name = par("cex"),
                          risk.table.cex.number = par("cex"),
                          risk.table.adj.name = 1,
                          risk.table.adj.number = NA,
                          legend = TRUE,
                          legend.x = "topright",
                          legend.y = NULL,
                          legend.legend = strata.names,
                          legend.col = col,
                          legend.border = "black",
                          legend.lty = 1,
                          legend.lwd,
                          legend.pch,
                          legend.angle = 45,
                          legend.bty = "0",
                          legend.bg = par("bg"),
                          legend.pt.bg = NA,
                          legend.cex = 1,
                          legend.pt.cex = legend.cex,
                          legend.pt.lwd = legend.lwd,
                          legend.xjust = 0,
                          legend.yjust = 1,
                          legend.x.intersp = 1,
                          legend.y.intersp = 1,
                          legend.adj = c(0, 0.5),
                          legend.text.width = NULL,
                          legend.text.col = par("col"),
                          legend.text.font = NULL,
                          legend.trace = FALSE,
                          legend.plot = TRUE,
                          legend.ncol = 1,
                          legend.horiz = FALSE,
                          legend.title = NULL,
                          legend.inset = 0,
                          legend.xpd,
                          legend.title.col = legend.text.col,
                          legend.title.adj = 0.5,
                          legend.seg.len = 2) {
  
  require(survival)
  require(stringr)
  
  formula <- gsub("[[:space:]]+", " ", paste(deparse(substitute(fit)), collapse = " "))
  
  if(!str_detect(formula, "~")) stop(paste("formula needs to be directly displayed in the function and not in the object '",
                                           formula, "'", sep = ""))
  
  if(any(!(str_detect(unlist(str_split(formula, ","))[2], fixed(")"))))
     & Trunc ==FALSE) stop("Start-stop notation should only be used for left-truncated data")
  
  if(!any(!(str_detect(unlist(str_split(formula, ","))[2], fixed(")"))))
     & Trunc == TRUE) stop("Start-stop notation should be used for left-truncated data")
  
  if(!str_detect(formula, "survfit.")) fit <- survfit(fit, data = data)
  
  varnames <- gsub("\\<survfit\\>|\\<Surv\\>|\\(|\\)| \\?\\(.*\\)", "", formula)
  varnames <- gsub("[[:space:]]+", "", unlist(strsplit(varnames, "[[:space:]]*(\\+|,|~)[[:space:]]*")))
  varnames <- varnames[!str_detect(varnames, "data=")]
  if(any(!setdiff(varnames, "1") %in% names(data))) stop(paste(paste(setdiff(varnames, "1")[!setdiff(varnames, "1") %in% names(data)], 
                                                                     collapse = " and "), 
                                                               " not in ", 
                                                               deparse(substitute(data)), "\n", sep = ""))
  
  n.group <- ifelse("strata" %in% names(fit), length(fit$strata), 1)
  
  if (n.group == 1) varnames <- varnames[-length(varnames)]
  
  if (is.null(col)) col <- rainbow(n.group)
  if (length(col[!is.na(col)]) != n.group) stop("One color per survival curve is necessary.\n")
  
  if(is.null(xlim) & !Trunc) xlim <- c(-1, 10*(floor(1+max(data[, varnames[1]])/10)))
  if(is.null(xlim) & Trunc) xlim <- c(-1, 10*(floor(1+max(data[, varnames[1]], data[, varnames[2]])/10)))
  
  if (plot.margin)
  {
    if(risk.table) par(mar = c(6 + n.group, 4.5, 4, 3))
    if(!risk.table) par(mar = c(6, 4.5, 4, 3))
  }
  
  if(!is.null(strata.names) & length(strata.names) != n.group) stop("strata.names don't fit the number of strata")
  if(is.null(strata.names) & n.group > 1) strata.names <- names(fit$strata)
  if(is.null(strata.names) & n.group == 1) strata.names <- "# at risk"
  
  plot(x = fit,
       conf.int = conf.int,
       main = main,
       xlab = list(xlab, cex = 1.2),
       ylab = list(ylab, cex = 1.2),
       xlim = xlim,
       ylim = ylim,
       # type = "l",
       lwd = lwd,
       col = col,
       lty = lty)
  
  if(n.group > 1 & legend == TRUE)
  {
    legend(x = legend.x,
           y = legend.y,
           legend= legend.legend,
           col = col,
           border = legend.border,
           lty = legend.lty,
           lwd = legend.lwd,
           pch = legend.pch,
           angle = legend.angle,
           bty = legend.bty,
           bg = legend.bg,
           pt.bg = legend.pt.bg,
           cex = legend.cex,
           pt.cex = legend.pt.cex,
           pt.lwd = legend.pt.lwd,
           xjust = legend.xjust,
           yjust = legend.yjust,
           x.intersp = legend.x.intersp,
           y.intersp = legend.y.intersp,
           adj = legend.adj,
           text.width = legend.text.width,
           text.col = legend.text.col,
           text.font = legend.text.font,
           trace = legend.trace,
           plot = legend.plot,
           ncol = legend.ncol,
           horiz = legend.horiz,
           title = legend.title,
           inset = legend.inset,
           xpd = legend.xpd,
           title.col = legend.title.col,
           title.adj = legend.title.adj,
           seg.len = legend.seg.len
    )
  }
  
  if (is.null(risk.table.times.print)) risk.table.times.print <- axis(1, labels = FALSE, tick = FALSE)
  at.risk <- c()
  
  if (n.group == 1)
  {
    if (!Trunc)
    {
      at.risk <- c(summary(fit, risk.table.times.print)$n.risk,
                   rep(0, length(risk.table.times.print)-length(summary(fit, risk.table.times.print)$n.risk)))
    } else {
      for (t in risk.table.times.print)
      {
        at.risk <- c(at.risk, nrow(data[which(data[, varnames[1]] <= t & data[, varnames[2]] >= t), ]))
      }
    }
  } else {
    for (i in 1:n.group)
    {
      if (!Trunc)
      {
        at.risk <- c(at.risk,
                     summary(fit, risk.table.times.print)$n.risk[summary(fit, risk.table.times.print)$strata==names(fit$strata)[i]],
                     rep(0, length(risk.table.times.print)-length(summary(fit, risk.table.times.print)$n.risk[summary(fit, risk.table.times.print)$strata==names(fit$strata)[i]])))
      }
    }
    if (Trunc)
    {
      covariate <- varnames[-c(1, 2, 3)]
      
      strata <- as.data.frame(matrix(ncol = length(covariate), nrow = length(names(table(data[, covariate[length(covariate)]])))))
      names(strata) <- covariate
      strata[, length(covariate)] <- names(table(data[, covariate[length(covariate)]]))
      
      if (length(covariate) > 1)
      {
        for (c in (length(covariate)-1):1)
        {
          x <- data.frame()
          for (i in 1:length(names(table(data[, covariate[c]]))))
          {
            x <- rbind(x, strata)
          }
          strata <- x
          strata[, c] <- rep(names(table(data[, covariate[c]])), each = nrow(strata)/length(names(table(data[, covariate[c]]))))
        }
      }
      
      for (s in 1:n.group)
      {
        for (t in risk.table.times.print)
        {
          DB_strata <- data
          for (c in 1:length(covariate))
          {
            DB_strata <- DB_strata[which(DB_strata[, covariate[c]] == strata[s, c]), ]
          }
          at.risk <- c(at.risk, nrow(DB_strata[which(DB_strata[, varnames[1]] <= t & DB_strata[, varnames[2]] >= t), ]))
        }
      }
    }
  }
  at.risk <- matrix(at.risk, nrow = n.group, byrow = TRUE)
  
  
  if(risk.table)
  {
    if(is.null(risk.table.at.name)) risk.table.at.name <- -0.3 * (risk.table.times.print[2] - risk.table.times.print[1]) + risk.table.times.print[1]
    if(is.null(risk.table.at.number)) risk.table.at.number <- risk.table.times.print
    for (i in 1:n.group)
    {
      mtext(side = 1, at = risk.table.at.name, line = i + 4, strata.names[i], cex = risk.table.cex.name, 
            adj = 1)
      mtext(side = 1, at = risk.table.at.number, line = i + 4, at.risk[i, ], cex = risk.table.cex.number)
    }
  }
}

# competing risks ----

plot_km_trunc_CR <- function(fit,
                             fit2 = NULL,
                             data,
                             data2 = NULL,
                             Trunc = FALSE,
                             Trunc2 = NULL,
                             strata.names = NULL,
                             strata.names2 = NULL,
                             conf.int = TRUE,
                             main = "",
                             xlab = "Time",
                             ylab = "Cumulative incidence function",
                             xlim = NULL,
                             ylim = c(-0.05, 1.05),
                             col = NULL,
                             lwd = 1,
                             plot.margin = TRUE,
                             risk.table = TRUE,
                             risk.table.times.print = NULL,
                             risk.table.at.name = NULL,
                             risk.table.at.number = NULL,
                             risk.table.cex.name = par("cex"),
                             risk.table.cex.number = par("cex"),
                             risk.table.adj.name = 1,
                             risk.table.adj.number = NA,
                             legend = TRUE,
                             legend.x = "topright",
                             legend.y = NULL,
                             legend.legend = strata.names,
                             legend.col = col,
                             legend.border = "black",
                             legend.lty = 1,
                             legend.lwd,
                             legend.pch,
                             legend.angle = 45,
                             legend.bty = "0",
                             legend.bg = par("bg"),
                             legend.pt.bg = NA,
                             legend.cex = 1,
                             legend.pt.cex = legend.cex,
                             legend.pt.lwd = legend.lwd,
                             legend.xjust = 0,
                             legend.yjust = 1,
                             legend.x.intersp = 1,
                             legend.y.intersp = 1,
                             legend.adj = c(0, 0.5),
                             legend.text.width = NULL,
                             legend.text.col = par("col"),
                             legend.text.font = NULL,
                             legend.trace = FALSE,
                             legend.plot = TRUE,
                             legend.ncol = 1,
                             legend.horiz = FALSE,
                             legend.title = NULL,
                             legend.inset = 0,
                             legend.xpd,
                             legend.title.col = legend.text.col,
                             legend.title.adj = 0.5,
                             legend.seg.len = 2) {
  require(survival)
  require(stringr)
  
  formula <- gsub("[[:space:]]+", " ", paste(deparse(substitute(fit)), collapse = " "))
  
  if(!str_detect(formula, "~")) stop(paste("formula in fit needs to be directly displayed in the function and not in the object '",
                                           formula, "'", sep = ""))
  
  if(any(!(str_detect(unlist(str_split(formula, ","))[2], fixed(")"))))
     & Trunc ==FALSE) stop("Start-stop notation should only be used for left-truncated data (1st model)")
  
  if(!any(!(str_detect(unlist(str_split(formula, ","))[2], fixed(")"))))
     & Trunc == TRUE) stop("Start-stop notation should be used for left-truncated data (1st model)")
  
  if(!str_detect(formula, "survfit.")) fit <- survfit(fit, data = data)
  
  varnames <- gsub("\\<survfit\\>|\\<Surv\\>|\\(|\\)| \\?\\(.*\\)", "", formula)
  varnames <- gsub("[[:space:]]+", "", unlist(strsplit(varnames, "[[:space:]]*(\\+|,|~)[[:space:]]*")))
  varnames <- varnames[!str_detect(varnames, "data=") & !str_detect(varnames, "id=")]
  if(any(!setdiff(varnames, "1") %in% names(data))) stop(paste(paste(setdiff(varnames, "1")[!setdiff(varnames, "1") %in% names(data)], 
                                                                     collapse = " and "), 
                                                               " not in ", 
                                                               deparse(substitute(data)), "\n", sep = ""))
  
  n.group <- ifelse("strata" %in% names(fit), length(fit$strata), 1)
  
  if (n.group == 1) varnames <- varnames[-length(varnames)]
  
  if(is.null(fit2) & is.null(col)) col <- rainbow(n.group)
  if(is.null(fit2) & length(col[!is.na(col)]) != n.group) stop("One color per survival curve is necessary.\n")
  
  if(is.null(fit2) & is.null(xlim) & !Trunc) xlim <- c(-1, 10*(floor(1+max(data[, varnames[1]])/10)))
  if(is.null(fit2) & is.null(xlim) & Trunc) xlim <- c(-1, 10*(floor(1+max(data[, varnames[1]], data[, varnames[2]])/10)))
  
  if (plot.margin)
  {
    if (is.null(fit2) & risk.table) par(mar = c(6 + n.group, 4.5, 4, 3))
    if (is.null(fit2) & !risk.table) par(mar = c(6, 4.5, 4, 3))
  }
  
  if(is.null(fit2) & !is.null(strata.names) & length(strata.names) != n.group) stop("strata.names don't fit the number of strata")
  if(is.null(fit2) & is.null(strata.names) & n.group > 1) strata.names <- names(fit$strata)
  if(is.null(fit2) & is.null(strata.names) & n.group == 1) strata.names <- "# at risk"
  
  if (!is.null(fit2))
  {
    if (is.null(data2)) data2 <- data
    formula2 <- gsub("[[:space:]]+", " ", paste(deparse(substitute(fit2)), collapse = " "))
    
    if(!str_detect(formula2, "~")) stop(paste("formula in fit2 needs to be directly displayed in the function and not in the object '",
                                              formula2, "'", sep = ""))
    
    if(any(!(str_detect(unlist(str_split(formula2, ","))[2], fixed(")"))))
       & Trunc2 == FALSE) stop("Start-stop notation should only be used for left-truncated data (2nd model)")
    
    if(!any(!(str_detect(unlist(str_split(formula2, ","))[2], fixed(")"))))
       & Trunc2 == TRUE) stop("Start-stop notation should be used for left-truncated data (2nd model)")
    
    if(!str_detect(formula2, "survfit.")) fit2 <- survfit(fit2, data = data2)
    
    varnames2 <- gsub("\\<survfit\\>|\\<Surv\\>|\\(|\\)| \\?\\(.*\\)", "", formula2)
    varnames2 <- gsub("[[:space:]]+", "", unlist(strsplit(varnames2, "[[:space:]]*(\\+|,|~)[[:space:]]*")))
    varnames2 <- varnames2[!str_detect(varnames2, "data=")]
    if(any(!setdiff(varnames2, "1") %in% names(data2))) stop(paste(paste(setdiff(varnames2, "1")[!setdiff(varnames2, "1") %in% names(data2)], 
                                                                         collapse = " and "), 
                                                                   " not in ", 
                                                                   deparse(substitute(data2)), "\n", sep = ""))
    if ("strata" %in% c(names(fit), names(fit2))) stop("Double model method doesn't handle covariates")
    
    n.group2 <- 1
    
    varnames2 <- varnames2[-length(varnames2)]
    
    if(is.null(col)) col <- rainbow(1)
    if(length(col[!is.na(col)]) != 1) stop("One color is enough for double model method.\n")
    
    if(is.null(xlim) & !Trunc) xlim <- c(-1, 10*(floor(1+max(data[, varnames[1]], data2[, varnames2[1]])/10)))
    if(is.null(xlim) & Trunc) xlim <- c(-1, 10*(floor(1+max(data[, varnames[1]], data[, varnames[2]], data2[, varnames2[1]], data2[, varnames2[2]])/10)))
    
    if (plot.margin)
    {
      if (risk.table) par(mar = c(6 + 2, 4.5, 4, 3))
      if (!risk.table) par(mar = c(6, 4.5, 4, 3))
    }
    
    
    if( !is.null(strata.names) & length(strata.names) != 1) stop("One strata.name maximum is allowed for double model method")
    if( !is.null(strata.names2) & length(strata.names2) != 1) stop("One strata.name2 maximum is allowed for double model method")
    if(is.null(strata.names)) strata.names <- "Outcome 1"
    if(is.null(strata.names2)) strata.names2 <- "Outcome 2"
  }
  
  if (is.null(fit2))
  {
    if(Trunc) if(!is.factor(data[, varnames[3]])) stop("Status covariate should be a factor for competing risks")
    if(!Trunc) if(!is.factor(data[, varnames[2]])) stop("Status covariate should be a factor for competing risks")
  }
  
  if (!is.null(fit2))
  {
    stepfun_fit <- stepfun(fit$time, c(1, fit$surv))
    stepfun_fit2 <- stepfun(fit2$time, c(1, fit2$surv))
  }
  
  
  if (is.null(fit2))
  {
    plot(x = fit,
         conf.int = FALSE,
         main = main,
         xlab = list(xlab, cex = 1.2),
         ylab = list(ylab, cex = 1.2),
         xlim = xlim,
         ylim = ylim,
         # type = "l",
         col = rep(col, n.group),
         lty = rep(c(1, 2), each = n.group),
         lwd = lwd)
  } else {
    plot(cbind(fit$time[fit$time <= xlim[2]],
               cumsum(c(1, (stepfun_fit(fit$time)*stepfun_fit2(fit$time))[1:(length(fit$time)-1)])* fit$n.event/fit$n.risk)[fit$time <= xlim[2]]),
         col = col,
         main = main,
         type = "s",
         xlab = list(xlab, cex = 1.2),
         ylab = list(ylab, cex = 1.2),
         xlim = xlim,
         ylim = ylim,
         lty = 1,
         lwd = lwd)
    lines(cbind(fit2$time[fit$time <= xlim[2]], 
                cumsum(c(1, (stepfun_fit(fit2$time)*stepfun_fit2(fit2$time))[1:(length(fit2$time)-1)])* fit2$n.event/fit2$n.risk)[fit$time <= xlim[2]]),
          col = col, lty = 2, type = "s")
  }
  
  if (conf.int & is.null(fit2))
  {
    times <- fit$time[fit$time <= xlim[2]]
    upper2 <- fit$upper[, 2][fit$time <= xlim[2]]
    lower2 <- fit$lower[, 2][fit$time <= xlim[2]]
    upper3 <- fit$upper[, 3][fit$time <= xlim[2]]
    lower3 <- fit$lower[, 3][fit$time <= xlim[2]]
    indic <- c(1, which(times[1:(length(times) - 1)] > times[2:length(times)])+1, length(times)+1)
    
    for (i in 1:n.group)
    {
      lines(times[indic[i]:(indic[i+1]-1)], 
            upper2[indic[i]:(indic[i+1]-1)], 
            col = col[i], lty = 3, type = "s", lwd = lwd)
      lines(times[indic[i]:(indic[i+1]-1)], 
            lower2[indic[i]:(indic[i+1]-1)], 
            col = col[i], lty = 3, type = "s", lwd = lwd)
      
      lines(times[indic[i]:(indic[i+1]-1)], 
            upper3[indic[i]:(indic[i+1]-1)], 
            col = col[i], lty = 3, type = "s", lwd = lwd)
      lines(times[indic[i]:(indic[i+1]-1)], 
            lower3[indic[i]:(indic[i+1]-1)], 
            col = col[i], lty = 3, type = "s", lwd = lwd)
    }
  }
  
  if (is.null(fit2) & risk.table)
  {
    if (is.null(risk.table.times.print)) risk.table.times.print <- axis(1, labels = FALSE, tick = FALSE)
    at.risk <- c()
    
    if (n.group == 1)
    {
      if (!Trunc)
      {
        at.risk <- c(summary(fit, risk.table.times.print)$n.risk[, 1],
                     rep(0, length(risk.table.times.print)-length(summary(fit, risk.table.times.print)$n.risk[, 1])))
      } else {
        for (t in risk.table.times.print)
        {
          at.risk <- c(at.risk, nrow(data[which(data[, varnames[1]] <= t & data[, varnames[2]] >= t), ]))
        }
      }
    } else {
      for (i in 1:n.group)
      {
        if (!Trunc)
        {
          at.risk <- c(at.risk,
                       summary(fit, risk.table.times.print)$n.risk[, 1][summary(fit, risk.table.times.print)$strata==names(fit$strata)[i]],
                       rep(0, length(risk.table.times.print)-length(summary(fit, risk.table.times.print)$n.risk[, 1][summary(fit, risk.table.times.print)$strata==names(fit$strata)[i]])))
        }
      }
      if (Trunc)
      {
        covariate <- varnames[-c(1, 2, 3)]
        
        strata <- as.data.frame(matrix(ncol = length(covariate), nrow = length(names(table(data[, covariate[length(covariate)]])))))
        names(strata) <- covariate
        strata[, length(covariate)] <- names(table(data[, covariate[length(covariate)]]))
        
        if (length(covariate) > 1)
        {
          for (c in (length(covariate)-1):1)
          {
            x <- data.frame()
            for (i in 1:length(names(table(data[, covariate[c]]))))
            {
              x <- rbind(x, strata)
            }
            strata <- x
            strata[, c] <- rep(names(table(data[, covariate[c]])), each = nrow(strata)/length(names(table(data[, covariate[c]]))))
          }
        }
        
        for (s in 1:n.group)
        {
          for (t in risk.table.times.print)
          {
            DB_strata <- data
            for (c in 1:length(covariate))
            {
              DB_strata <- DB_strata[which(DB_strata[, covariate[c]] == strata[s, c]), ]
            }
            at.risk <- c(at.risk, nrow(DB_strata[which(DB_strata[, varnames[1]] <= t & DB_strata[, varnames[2]] >= t), ]))
          }
        }
      }
    }
    at.risk <- matrix(at.risk, nrow = n.group, byrow = TRUE)
    
    if(is.null(risk.table.at.name)) risk.table.at.name <- -0.3 * (risk.table.times.print[2] - risk.table.times.print[1]) + risk.table.times.print[1]
    if(is.null(risk.table.at.number)) risk.table.at.number <- risk.table.times.print
    for (i in 1:n.group)
    {
      mtext(side = 1, at = risk.table.at.name, line = i + 4, strata.names[i], cex = risk.table.cex.name, 
            adj = 1)
      mtext(side = 1, at = risk.table.at.number, line = i + 4, at.risk[i, ], cex = risk.table.cex.number)
    }
  }
  
  if (!is.null(fit2) & risk.table)
  {
    if (is.null(risk.table.times.print)) risk.table.times.print <- axis(1, labels = FALSE, tick = FALSE)
    at.risk <- c()
    if (Trunc)
    {
      for (t in risk.table.times.print)
      {
        at.risk <- c(at.risk, nrow(data[which(data[, varnames[1]] <= t & data[, varnames[2]] >= t), ]))
      }
      
    } else {
      at.risk <- c(at.risk,
                   summary(fit, risk.table.times.print)$n.risk,
                   rep(0, length(risk.table.times.print)-length(summary(fit, risk.table.times.print)$n.risk)))
    }
    
    if (!is.null(fit2) & risk.table)
    {
      if (Trunc2)
      {
        for (t in risk.table.times.print)
        {
          at.risk <- c(at.risk, nrow(data2[which(data2[, varnames2[1]] <= t & data2[, varnames2[2]] >= t), ]))
        }
      } else {
        at.risk <- c(at.risk,
                     summary(fit2, risk.table.times.print)$n.risk,
                     rep(0, length(risk.table.times.print)-length(summary(fit2, risk.table.times.print)$n.risk)))
      }
    }
    
    at.risk <- matrix(at.risk, nrow = 2, byrow = TRUE)
    
    if(is.null(risk.table.at.name)) risk.table.at.name <- -0.3 * (risk.table.times.print[2] - risk.table.times.print[1]) + risk.table.times.print[1]
    if(is.null(risk.table.at.number)) risk.table.at.number <- risk.table.times.print
    for (i in 1:2)
    {
      mtext(side = 1, at = risk.table.at.name, line = i + 4, c(strata.names, strata.names2)[2], cex = risk.table.cex.name, 
            adj = 1)
      mtext(side = 1, at = risk.table.at.number, line = i + 4, at.risk[i, ], cex = risk.table.cex.number)
    }
  }
  
  if (conf.int & !is.null(fit2))
  {
  }
}

# rec_event ----

rec_event <- function(id, # id of patient 
                      start, # start of observation
                      stop, # stop of observation
                      status, # status of event at the end of observation
                      terminal, # status of competing risk at the end of observation
                      recdata, # database for event cox model
                      recdata2, # database for competing risk cox model. Set to NULL if you want recdata to be used here instead.
                      conf_int=TRUE){ # set it to TRUE if you want estimations of confidence interval to be calculated.
  
  recdata <- recdata[, c(id, start, stop, status)]
  names(recdata) <- c("id", "start", "stop", "status")
  
  recdata2 <- recdata2[, c(id, start, stop, terminal)]
  names(recdata2) <- c("id", "start", "stop", "terminal")
  
  NPfit <- coxph(Surv(start,stop,status) ~ 1, data=recdata)
  expecNP <- survfit(NPfit,type="aalen")
  
  NPfit <- coxph(Surv(start,stop,terminal) ~ 1, data=recdata2)
  survT <- survfit(NPfit)
  
  surv <- stepfun(survT$time, c(1, survT$surv))
  
  N <- length(unique(recdata$id))
  
  if(conf_int)
  {
    Newdata <- survSplit(recdata, cut=expecNP$time, end="stop", event="status", start="start")
    lambda=cbind(expecNP$time, cumsum(expecNP$n.event/expecNP$n.risk))
    residRec=matrix(rep(0, length(expecNP$time)*N), ncol=length(expecNP$time))
    j <- 1
    for (i in unique(recdata$id))
    { 
      indiv <- Newdata[which(Newdata$id == i), ]
      lambdai <- matrix(lambda[which(lambda[,1]<=max(indiv$stop)), ], ncol = 2) 
      d <- dim(lambdai)[1] #\hat M_i(t) 
      
      if (nrow(indiv)>1)
      {
        if(indiv$stop[nrow(indiv)] < (indiv$stop[nrow(indiv)-1]+0.00001))
        {
          indiv$status[nrow(indiv)-1] <- indiv$status[nrow(indiv)]
          indiv <- indiv[-nrow(indiv), ]
        }
      }
      
      indiv <- indiv[which(indiv$star + 0.00000001 < indiv$stop), ]
      
      residRec[j, ] <- c(cumsum(indiv$status) - lambdai[, 2],
                         rep(max(cumsum(indiv$status)) - max(lambdai[,2]),
                             length(expecNP$time) - d))
      j <- j + 1
    }
    
    event <- stepfun(survT$time, c(0, cumsum(survT$n.event/survT$n.risk)))
    
    NewdataD <- survSplit(recdata2, cut = expecNP$time, end = "stop", event = "terminal", start = "start")
    lambdaT <- cbind(expecNP$time, event(expecNP$time)) 
    resid <- matrix(rep(0, length(expecNP$time)*N), ncol=length(expecNP$time)) 
    j <- 1
    for (i in unique(recdata2$id)) 
    { 
      indiv <- NewdataD[which(NewdataD$id==i), ]
      lambdaiT <- matrix(lambdaT[which(lambdaT[, 1] <= max(indiv$stop)), ], ncol = 2) 
      dT <- dim(lambdaiT)[1] #\hat M^D_i(t) 
      
      if (nrow(indiv)>1)
      {
        if(indiv$stop[nrow(indiv)] < (indiv$stop[nrow(indiv)-1]+0.00001))
        {
          indiv$terminal[nrow(indiv)-1] <- indiv$terminal[nrow(indiv)]
          indiv <- indiv[-nrow(indiv), ]
        }
      }
      
      indiv <- indiv[which(indiv$star + 0.00000001 < indiv$stop), ]
      
      resid[j, ] <- c(rep(0, length(which(lambdaiT[, 1] < min(indiv$stop)))), indiv$terminal - lambdaiT[which(lambdaiT[, 1] %in% indiv$stop), 2],
                      rep(max(indiv$terminal) - max(lambdaiT[, 2]),
                          length(expecNP$time) - dT))
      j <- j + 1
    }
    
    Shat <- matrix(rep(surv(expecNP$time), N), byrow = TRUE, ncol = length(expecNP$time))
    Ybar <- matrix(rep(expecNP$n.risk, N), byrow = TRUE, ncol = length(expecNP$time))
    muT <- cumsum(surv(expecNP$time)*expecNP$n.event/expecNP$n.risk)
    mu <- matrix(rep(muT, N), byrow = TRUE, ncol = length(expecNP$time)) 
    dMD <- t(apply(cbind(rep(0, N), resid), 1, diff))
    dM <- t(apply(cbind(rep(0, N), residRec), 1, diff))
    psi2 <- N*(t(apply(Shat*dM/Ybar, 1, cumsum)) - mu*t(apply(dMD/Ybar, 1, cumsum)) + t(apply(mu*dMD/Ybar, 1, cumsum)))^2
    sigma2 <- apply(psi2, 2, sum)
    
    CIleft <- muT*exp((1/sqrt(N))*qnorm(0.975)*sqrt(sigma2)/muT) 
    CIright <- muT*exp((-1/sqrt(N))*qnorm(0.975)*sqrt(sigma2)/muT)
    
    lower <- stepfun(survT$time, c(1, survT$lower))
    upper <- stepfun(survT$time, c(1, survT$upper))
    
    result <- list(expecNP, survT, CIleft, CIright, expecNP$time, cumsum(surv(expecNP$time)*expecNP$n.event/expecNP$n.risk), 1-surv(expecNP$time), 1-lower(expecNP$time), 1-upper(expecNP$time))
    names(result) <- c("expecNP", "survT", "CIleft", "CIright", "time", "mean_number", "cif_surv", "CIleft_surv", "CIright_surv")
  }
  class(result) <- "rec_event"
  return(result)
}

plot.rec_event <- function(x,
                           conf.int = TRUE,
                           risk.table = TRUE,
                           risk.table.margin = TRUE,
                           risk.table.times.print = NULL,
                           risk.table.name = NULL,
                           risk.table.at.name = NULL,
                           risk.table.at.number = NULL,
                           risk.table.cex.name = par("cex"),
                           risk.table.cex.number = par("cex"),
                           risk.table.adj.name = 1,
                           risk.table.adj.number = NA,
                           main = "", 
                           xlab="Time", 
                           ylab="Mean number of recurrent events",
                           col, lwd, ...) {
  
  if (risk.table & risk.table.margin)
  {
    if (!is.null(risk.table.margin))
    {
      par(mar = c(6 + risk.table.margin, 4.5, 4, 3))
    } else {
      par(mar = c(6 + n.group, 4.5, 4, 3))
    }
  } else {
    par(mar = c(6, 4.5, 4, 3))
  }
  
  plot(x$time,
       x$mean_number,
       type="s",
       xlab = list(xlab, cex = 1.2),
       ylab = list(ylab, cex = 1.2),
       main=main,
       lty=1,
       col = col,
       lwd = lwd,
       ...)
  if (conf.int)
  {
    lines(x$time, x$CIleft, type="s", col=col, lty=2, lwd = lwd)
    lines(x$time, x$CIright, type="s", col=col, lty=2, lwd = lwd)
  }
  
  if (is.null(risk.table.times.print)) risk.table.times.print <- axis(1, labels = FALSE, tick = FALSE)
  
  at.risk <- matrix(c(summary(x$expecNP, risk.table.times.print)$n.risk,
                      rep(0, length(risk.table.times.print)-length(summary(x$expecNP, risk.table.times.print)$n.risk))), 
                    nrow = 1, byrow = TRUE)
  if(risk.table)
  {
    if(is.null(risk.table.at.name)) risk.table.at.name <- -0.3 * (risk.table.times.print[2] - risk.table.times.print[1]) + risk.table.times.print[1]
    if(is.null(risk.table.at.number)) risk.table.at.number <- risk.table.times.print
    if(is.null(risk.table.name)) risk.table.name <- "# at risk"
    mtext(side = 1, at = risk.table.at.name, line = 5, risk.table.name, cex = risk.table.cex.name, 
          adj = 1)
    mtext(side = 1, at = risk.table.at.number, line = 5, at.risk, cex = risk.table.cex.number)
  }
}