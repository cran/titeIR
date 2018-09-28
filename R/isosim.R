#' Simulate TITE-IR designs
#'
#' Simulates trials based on the TITE-IR design.
#'
#' @param PI A vector of true toxicity probabilities at each dose
#' @param target Target DLT rate
#' @param n Sample size of the trial
#' @param nsim Number of trial replicates
#' @param obswin The observation window with respect to which the MTD is defined
#' @param rate Patient arrival rate: expected number of arrivals per observation window
#' @param safety The safety factor to prevent overly aggressive escalation
#' @param accrual Specify the accrual distribution. Can be either \code{"poisson"} or \code{"fixed"}. Partial strings are also acceptable.
#' @param restrict If \code{TRUE}, do not allow escalation immediately after a toxic outcome (require coherent escalation)
#'
#' @return Object of type \code{isotite} which provides results from TITE-IR simulations
#' @aliases titeIR
#' @seealso \code{\link{isotitedose}} for dose recommendation
#' @importFrom Iso pava
#' @examples
#' isotitesim(PI = c(0.05, 0.10, 0.20, 0.30, 0.50, 0.70),
#'  target = 1/3, n = 24, nsim = 10, obswin = 6, rate = 12)
#' @export
isotitesim <- function(PI,
                       target,
                       n,
                       nsim,
                       obswin = 1,
                       rate = 1,
                       safety = 0.05,
                       accrual = "poisson",
                       restrict = TRUE){
  doselevels <- length(PI)
  totalTime <- rep(0, nsim)
  earlyStop <- rep(0, nsim)
  MTD <- rep(0, nsim)
  totpatients = matrix(0, nsim, ncol=doselevels)
  simassignments <- matrix(0, nrow=nsim, ncol=n)
  simtoxes <- matrix(0, nrow=nsim, ncol=n)
  allTimes <- matrix(0, nrow=nsim, ncol=n)

  accmatch <- pmatch(accrual, c("fixed", "poisson"))

  if(is.na(accmatch)) stop("accrual must be one of 'fixed' or 'poisson'")


  for(k in 1:nsim){

    DLTs <- rep(0, doselevels)
    initDLTs <- DLTs

    initpatients <- rep(0, doselevels)

    cohortCount <- 1
    currdose = 1
    if(accmatch == 2){
      recTime <- cumsum(stats::rexp(n, rate / obswin))
    }else if(accmatch == 1){
      recTime <- (1:n)*obswin/rate
    }else{
      stop("accrual must be either 'fixed' or 'poisson'")
    }
    recTime <- recTime
    currTime <- 0
    assignment <- rep(0, n)

    results <- rep(0, n)
    eventTime <- stats::runif(n, 0, obswin)

    assignment[1] <- 1
    results[1] <- stats::rbinom(1, 1, PI[1])

    if(results[1] == 0){
      eventTime[1] <- obswin
    }
    for(i in 2:n){
      tmpDLT <- initDLTs
      patients <- initpatients
      tmpDenom <- initDLTs
      tmpadj <- initDLTs
      currassign <- tmpDLT
      ### Sum over all previous patients, get the numerator for the pseudo-proportions

        currTime <- recTime[i]
        for(j in 1:(i - 1)){
          if((currTime >= recTime[j] + eventTime[j])){
            patients[assignment[j]] <- patients[assignment[j]] + 1
            tmpDLT[assignment[j]] <- tmpDLT[assignment[j]] + results[j]
            tmpDenom[assignment[j]] <- sum(assignment[j] == assignment)
          }else if(currTime < recTime[j] + eventTime[j]){
            tmpDLT[assignment[j]] <- tmpDLT[assignment[j]] + (safety+target)*(obswin - (currTime - recTime[j])) / obswin
            tmpDenom[assignment[j]] <- sum(assignment[j] == assignment)
          }
        }
        for(j in 1:length(initDLTs)){
          currassign[j] <- sum(assignment==j)
        }


        probvec <- rep(0,doselevels)
        # Placeholder probability if prior dose hasn't been tried yet
        for(j in 2:doselevels){
          if(currassign[j-1] == 0) probvec[j] <- 1
        }
        # get raw pseudo-probabilities
        probvec[currassign > 0] <- tmpDLT[currassign > 0]/currassign[currassign > 0]

        tmpProbs <- pava(probvec, w=currassign) #Calculate tentative estimates with isotonic regression


        # Keep coherence as per Cheung, Chappell
        coherent <- TRUE
        curRes <- results[assignment==currdose]

        if(restrict == TRUE){
          lastFull <- patients[currdose]
          if(lastFull==0){
            coherent <- TRUE
          }else{
            if(curRes[lastFull] == 1) coherent <- FALSE
          }
        }else{
          coherent <- TRUE
        }
        #coherent <- TRUE

        if((tmpProbs[currdose] < target)){
          if(currdose < doselevels){
            # If next dose is closer to target, and we have 3 or more assigned to current dose, then coherence is satisfied, then escalate
            if(((target - tmpProbs[currdose]) >= (tmpProbs[currdose+1] - target)) & (sum(assignment == currdose) > 2) & (coherent)) currdose <- currdose + 1
          }
        }else if(tmpProbs[currdose] >= target){
          if(currdose > 1){
            # If previous dose is closer to target, and we have at least 3 assigned to current dose, then de-escalate
            if(((target - tmpProbs[currdose - 1]) < (tmpProbs[currdose] - target)) & (sum(assignment == currdose) > 2)) currdose <- currdose - 1
          }
        }

        assignment[i] <- currdose
        results[i] <- stats::rbinom(1, 1, PI[currdose])

        if(results[i] == 0){
          eventTime[i] <- obswin
        }

    }

    DLTs <- initDLTs
    patients <- initpatients
    for(i in 1:n){
      DLTs[assignment[i]] <- DLTs[assignment[i]] + results[i]
      patients[assignment[i]] <- patients[assignment[i]] + 1
    }

    probvec <- DLTs/patients
    probvec[is.nan(probvec)] <- 1

    finalprob <- pava(probvec, w=patients)


    if(any(finalprob > target)){
      MTD[k] <- max(which.min(finalprob <= target) - 1, 1)
    }else{
      MTD[k] <- doselevels
    }

    totpatients[k,] <- patients
    totalTime[k] <- max(recTime[assignment > 0]) + obswin
    simassignments[k,] <- assignment
    if(earlyStop[k] == 1){
      simtoxes[k,] <- results*as.numeric(currTime > (recTime + eventTime))
    }else{
      simtoxes[k,] <- results
    }
    allTimes[k,] <- recTime

    #if(earlyStop[k] == 1) break
  }
  retlist <- list(Patients = totpatients,
                  Time = totalTime,
                  MTD = MTD,
                  Assignments = simassignments,
                  Toxicities = simtoxes,
                  target = target,
                  allTimes = allTimes,
                  PI = PI  )
  class(retlist) <- "isotite"
  return(retlist)
}


#' @export
summary.isotite <- function(object, ...){
  if(any(object$PI > object$target)){
    MTD = which.min(object$PI <= object$target) - 1
  }else{
    MTD = length(object$PI)
  }
  meanMTD <- mean(object$MTD == MTD)
  meanTox <- mean(apply(object$Toxicities, 1, sum))
  numAssigned <- mean(apply(object$Assignments > 0, 1, sum))
  duration <- mean(object$Time)
  belowMTD <- 100*mean(apply((object$Assignments < MTD) & (object$Assignments > 0), 1, sum) / apply(object$Assignments > 0, 1, sum))
  atMTD <- 100*mean(apply((object$Assignments == MTD) , 1, sum)/ apply(object$Assignments > 0, 1, sum))
  aboveMTD <- 100*mean(apply((object$Assignments > MTD) , 1, sum)/ apply(object$Assignments > 0, 1, sum))
  if(MTD > 0){
    Delta <- mean(pmax(object$PI[MTD] - object$PI[object$MTD], 0) + 2*pmax(object$PI[object$MTD] - object$PI[MTD], 0))
  }else{
    Delta <- mean(2*pmax(object$PI[object$MTD[object$MTD > 0]] - object$target, 0))
  }
  summ <- list(meanMTD = meanMTD, Delta = Delta, meanTox = meanTox, numAssigned = numAssigned, duration = duration,
               belowMTD = belowMTD, atMTD = atMTD, aboveMTD = aboveMTD, doseAssign = object$MTD, numDoses = ncol(object$Patients))
  class(summ) <- "summary.isotite"
  return(summ)
}

#' @export
print.summary.isotite <- function(x, ...){
  cat("Probability of Correct Dose:", x$meanMTD)
  cat("\nAverage Number of Toxicites per trial:", x$meanTox)
  cat("\nDuration:", round(x$duration, digits=2))
  cat("\nPercent Treated below MTD:", round(x$belowMTD, digits=1))
  cat("\nPercent Treated at MTD:", round(x$atMTD, digits=1))
  cat("\nPercent Treated above MTD:", round(x$aboveMTD, digits=1), "\n")
}

#' @export
print.isotite <- function(x, ...){
  nsim <- length(x$MTD)
  if(nsim > 1){
    doses <- ncol(x$Patients)
    Selected <- Ntox <- rep(0, ncol(x$Patients))
    for(i in 1:doses){
      Selected[i] <- mean(x$MTD == i)
      Ntox[i] <- sum(x$Toxicities[x$Assignments == i])/nsim

    }
    Nexpt <- apply(x$Patients, 2, mean)
    oc <- rbind( "Pr(DLT)" = x$PI, Selected, Nexpt, Ntox)
    colnames(oc) <- 1:doses
    print(round(oc, digits=2))
  }
}

