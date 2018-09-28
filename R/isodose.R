
#' Dose assignment for TITE-IR designs
#'
#' Calculate the next dose assignment for a TITE-IR design.
#'
#' @param followup A vector of followup times
#' @param DLT A vector of DLT results. \code{FALSE} or 0 is interpreted as no observed DLT and \code{TRUE} or 1 is interpreted as observed DLT.
#' @param assignment a vector of dose assignments. Doses should be labeled in consecutive integers from 1 to number of dose levels.
#' @param obswin The observation window with respect to which the MTD is defined.
#' @param doses An integer providing the number of doses.
#' @param target Target DLT rate
#' @param safety The safety factor to prevent overly aggressive escalation
#'
#' @return an integer specifying the recommended dose level
#' @seealso \code{\link{isotitesim}} for simulations
#' @examples
#' isotitedose(followup = c(6, 5, 4, 3, 2, 1), DLT = c(0, 0, 0, 0, 0, 0),
#'  assignment = c(1, 1, 1, 2, 2, 2), obswin = 6, doses = 6)
#' @export
isotitedose <- function(followup, DLT, assignment, obswin, doses, target = 1/3, safety = 0.05){
  n <- length(followup)
  currdose <- assignment[n]
  patients <- tmpDLT <- currassign <- rep(0, doses)
  for(i in 1:n){
    if((followup[i] >= obswin) | (DLT[i] == 1)){
      patients[assignment[i]] <- patients[assignment[i]] + 1
      tmpDLT[assignment[i]] <- tmpDLT[assignment[i]] + DLT[i]
      currassign[assignment[i]] <- sum(assignment[i] == assignment)
    }else if((followup[i] < obswin) & (DLT[i] == 0)){
      tmpDLT[assignment[i]] <- tmpDLT[assignment[i]] + (safety+target)*(obswin - followup[i]) / obswin
      currassign[assignment[i]] <- sum(assignment[i] == assignment)
    }
  }

  probvec <- rep(0,doses)
  # Placeholder probability if prior dose hasn't been tried yet
  for(j in 2:doses){
    if(currassign[j-1] == 0) probvec[j] <- 1
  }
  # get raw pseudo-probabilities
  probvec[currassign > 0] <- tmpDLT[currassign > 0]/currassign[currassign > 0]

  tmpProbs <- pava(probvec, w=currassign) #Calculate tentative estimates with isotonic regression


  # Keep coherence as per Cheung, Chappell
  coherent <- TRUE
  curRes <- DLT[assignment==currdose]
  restrict = TRUE
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
    if(currdose < doses){
      # If next dose is closer to target, and we have 3 or more assigned to current dose, then coherence is satisfied, then escalate
      if(((target - tmpProbs[currdose]) >= (tmpProbs[currdose+1] - target)) & (sum(assignment == currdose) > 2) & (coherent)) currdose <- currdose + 1
    }
  }else if(tmpProbs[currdose] >= target){
    if(currdose > 1){
      # If previous dose is closer to target, and we have at least 3 assigned to current dose, then de-escalate
      if(((target - tmpProbs[currdose - 1]) < (tmpProbs[currdose] - target)) & (sum(assignment == currdose) > 2)) currdose <- currdose - 1
    }
  }

  return(currdose)
}


