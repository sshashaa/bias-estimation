###############################################
###############################################
########## Adaptive Sampling (A) ##############
############# Oct 18 2021 #####################
############# Kimia Vahdat ####################
###############################################
###############################################

packages <- c("caret", "dplyr", "GA", "doParallel","extraDistr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

############# Libraries #####################
library(caret)
library(GA)
library(dplyr)
library(doParallel)
library(extraDistr)
library(data.table)
##### End loading packages##################

##### Registering cores (Clusters) ####
clusters<-20
cl<-makeCluster(clusters)
registerDoParallel(cl)
##### End ########################


##### A function to select best subset of features using different method ############
best_subset<- function(data_corr,k){
  require('caret')
  require('GA')
  ## Setting the parameters
  N<-10 # number of Macro-replication
  P<-1000 # number of post-processing
  B<-15 # maximum number of microreps
  q=1000
  maxiter_param <- 5000
  OOB_err<-0
  
  random_seed_train=seeds[,k]
  ## Start the time
  ptm<-Sys.time()
  
  ## Depending on what method we want, best_subset function runs one of the following if procedures below

    N<-10 # number of Macro-replication
    P<-1000 # number of post-processing
    B<-15 # maximum number of microreps
    q=1000
    ## creating different train sets for Our GA
    # train<-matrix(0,q,B)
    # for (j in 1:B) {
    #   set.seed(random_seed_train[j])
    #   train[,j]<-sample(1: nrow(data_corr), q,replace = T)
    # }
    ### Defining a new GA function with Ranking and Selection ####
    RSga<- function(type = c("binary", "real-valued", "permutation"), 
                    fitness, ...,
                    lower, upper, nBits,
                    population = gaControl(type)$population,
                    selection = gaControl(type)$selection,
                    crossover = gaControl(type)$crossover, 
                    mutation = gaControl(type)$mutation,
                    popSize = 50, 
                    pcrossover = 0.8, 
                    pmutation = 0.1, 
                    elitism = base::max(1, round(popSize*0.05)), 
                    updatePop = FALSE,
                    postFitness = NULL,
                    maxiter = 100,
                    run = maxiter,
                    maxFitness = Inf,
                    names = NULL,
                    suggestions = NULL,
                    optim = FALSE,
                    optimArgs = list(method = "L-BFGS-B", 
                                     poptim = 0.05,
                                     pressel = 0.5,
                                     control = list(fnscale = -1, maxit = 100)),
                    keepBest = FALSE,
                    parallel = FALSE,
                    monitor = if(interactive()) gaMonitor else FALSE,
                    seed = NULL) 
    {
      
      call <- match.call()
      
      type <- match.arg(type, choices = eval(formals(ga)$type))
      
      if(!is.function(population)) population <- get(population)
      if(!is.function(selection))  selection  <- get(selection)
      if(!is.function(crossover))  crossover  <- get(crossover)
      if(!is.function(mutation))   mutation   <- get(mutation)
      
      if(missing(fitness))
      { stop("A fitness function must be provided") }
      if(!is.function(fitness)) 
      { stop("A fitness function must be provided") }
      if(popSize < 10) 
      { warning("The population size is less than 10.") }
      if(maxiter < 1) 
      { stop("The maximum number of iterations must be at least 1.") }
      if(elitism > popSize) 
      { stop("The elitism cannot be larger that population size.") }
      elitism <- as.integer(elitism)
      if(pcrossover < 0 | pcrossover > 1)
      { stop("Probability of crossover must be between 0 and 1.") }
      if(is.numeric(pmutation))
      { 
        if(pmutation < 0 | pmutation > 1)
        { stop("If numeric probability of mutation must be between 0 and 1.") }
        else if(!is.function(population))
        { stop("pmutation must be a numeric value in (0,1) or a function.") }
      }
      
      # check for min and max arguments instead of lower and upper
      callArgs <- list(...)
      if(any("min" %in% names(callArgs)))
      {
        lower <- callArgs$min
        callArgs$min <- NULL
        warning("'min' arg is deprecated. Use 'lower' instead.")
      }
      if(any("max" %in% names(callArgs)))
      {
        upper <- callArgs$max
        callArgs$max <- NULL
        warning("'max' arg is deprecated. Use 'upper' instead.")
      }
      
      if(missing(lower) & missing(upper) & missing(nBits))
      { stop("A lower and upper range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!") }
      
      # check GA search type 
      switch(type, 
             "binary"      = { nBits <- as.vector(nBits)[1]
             lower <- upper <- NA
             nvars <- nBits 
             if(is.null(names))
               names <- paste0("x", 1:nvars)
             },
             "real-valued" = { lnames <- names(lower)
             unames <- names(upper)
             lower <- as.vector(lower)
             upper <- as.vector(upper)
             nBits <- NA
             if(length(lower) != length(upper))
               stop("lower and upper must be vector of the same length!")
             nvars <- length(upper)
             if(is.null(names) & !is.null(lnames))
               names <- lnames
             if(is.null(names) & !is.null(unames))
               names <- unames
             if(is.null(names))
               names <- paste0("x", 1:nvars)
             },
             "permutation" = { lower <- as.vector(lower)[1]
             upper <- as.vector(upper)[1]
             nBits <- NA
             nvars <- length(seq.int(lower,upper))
             if(is.null(names))
               names <- paste0("x", 1:nvars)
             }
      )
      
      # check suggestions
      if(is.null(suggestions))
      { suggestions <- matrix(nrow = 0, ncol = nvars) }
      else
      { if(is.vector(suggestions)) 
      { if(nvars > 1) suggestions <- matrix(suggestions, nrow = 1)
      else          suggestions <- matrix(suggestions, ncol = 1) }
        else
        { suggestions <- as.matrix(suggestions) }
        if(nvars != ncol(suggestions))
          stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
      }
      
      # check monitor arg
      if(is.logical(monitor))
      { if(monitor) monitor <- gaMonitor }
      if(is.null(monitor)) monitor <- FALSE
      
      # if optim merge provided and default args for optim()
      if(optim)
      { # merge default and provided parameters
        optimArgs.default <- eval(formals(ga)$optimArgs)
        optimArgs.default$control[names(optimArgs$control)] <- optimArgs$control
        optimArgs$control <- NULL
        optimArgs.default[names(optimArgs)] <- optimArgs
        optimArgs <- optimArgs.default; rm(optimArgs.default)
        if(any(optimArgs$method == c("L-BFGS-B", "Brent")))
        { optimArgs$lower <- lower
        optimArgs$upper <- upper }
        else
        { optimArgs$lower <- -Inf
        optimArgs$upper <- Inf }
        optimArgs$poptim <- min(max(0, optimArgs$poptim), 1)
        optimArgs$pressel <- min(max(0, optimArgs$pressel), 1)
        optimArgs$control$maxit <- as.integer(optimArgs$control$maxit)
        # ensure that optim maximise the fitness
        if(is.null(optimArgs$control$fnscale))
          optimArgs$control$fnscale <- -1
        if(optimArgs$control$fnscale > 0)
          optimArgs$control$fnscale <- -1*optimArgs$control$fnscale
      }
      
      # Start parallel computing (if needed)
      if(is.logical(parallel))
      { if(parallel) 
      { parallel <- startParallel(parallel)
      stopCluster <- TRUE }
        else
        { parallel <- stopCluster <- FALSE } 
      }
      else
      { stopCluster <- if(inherits(parallel, "cluster")) FALSE else TRUE
      parallel <- startParallel(parallel) 
      }
      on.exit(if(parallel & stopCluster)
        stopParallel(attr(parallel, "cluster")) )
      # define operator to use depending on parallel being TRUE or FALSE
      `%DO%` <- if(parallel && requireNamespace("doRNG", quietly = TRUE)) 
        doRNG::`%dorng%`
      else if(parallel) `%dopar%` else `%do%`
      # set seed for reproducibility  
      if(!is.null(seed)) set.seed(seed)
      i. <- NULL # dummy to trick R CMD check 
      # we need to keep track of performance of each iteration
      fitnessSummary <- matrix(as.double(NA), nrow = maxiter, ncol = 6)
      colnames(fitnessSummary) <- names(gaSummary(rnorm(10)))
      bestSol <- if(keepBest) vector(mode = "list", length = maxiter)
      else         list()
      # Fitness <- rep(NA, popSize)
      Fitness <- array(NA, dim=c(2+nrow(data_corr),popSize))
      # Fitness <- vector("list", length = popSize)
      # Fitness <- lapply(Fitness, function(x) NA)
      # Fitness <- list(fitness= rep(NA, popSize), fitvar= rep(NA, popSize), IF= array(NA, dim= c(nrow(data_corr),popSize)))
      object <- new("ga", 
                    call = call, 
                    type = type,
                    lower = lower, 
                    upper = upper, 
                    nBits = nBits, 
                    names = if(is.null(names)) character() else names,
                    popSize = popSize,
                    iter = 0, 
                    run = 1, 
                    maxiter = maxiter,
                    suggestions = suggestions,
                    population = matrix(), 
                    elitism = elitism, 
                    pcrossover = pcrossover, 
                    pmutation = if(is.numeric(pmutation)) pmutation else NA,
                    optim = optim,
                    fitness = Fitness[1,],
                    #fitvar = Fitness[2,],
                    fitder = Fitness[-c(1),],
                    summary = fitnessSummary,
                    bestSol = bestSol)
      
      if(maxiter == 0)
        return(object)
      
      # generate beginning population
      Pop <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
      ng <- min(nrow(suggestions), popSize)
      if(ng > 0) # use suggestion if provided
      { Pop[1:ng,] <- suggestions }
      # fill the rest with a random population
      if(popSize > ng)
      { Pop[(ng+1):popSize,] <- population(object)[1:(popSize-ng),] }
      object@population <- Pop
      
      # start iterations
      for(iter in seq_len(maxiter))
      {
        object@iter <- iter
        
        # evalute fitness function (when needed) 
        if(!parallel)
        { for(i in seq_len(popSize))
          if(is.na(Fitness[i]))
          { fit <- do.call(fitness, c(list(Pop[i,]), callArgs)) 
          if(updatePop)
            Pop[i,] <- attributes(fit)[[1]]
          Fitness[i] <- fit
          }
        }
        else
        { Fitness <- foreach(i. = seq_len(popSize), .combine = "cbind") %DO%
          { if(is.na(Fitness[1,i.])) 
            do.call(fitness, c(list(Pop[i.,]), callArgs)) 
            else                   
              Fitness[,i.] 
          }
        }
        #print(Fitness[,1:10])
        # update object
        object@population <- Pop
        object@fitness <- Fitness[1,]
        #object@fitvar <- Fitness[2,]
        object@fitder <- Fitness[-c(1),]
        # Local search optimisation
        if(optim & (type == "real-valued"))
        {
          if(optimArgs$poptim > runif(1))
          { # perform local search from random selected solution
            # with prob proportional to fitness
            i <- sample(1:popSize, size = 1, 
                        prob = optimProbsel(Fitness, q = optimArgs$pressel))
            # run local search
            opt <- try(suppressWarnings(
              do.call(stats::optim, 
                      c(list(fn = fitness,
                             par = Pop[i,],
                             method = optimArgs$method,
                             lower = optimArgs$lower,
                             upper = optimArgs$upper,
                             control = optimArgs$control), 
                        callArgs))
            ), silent = TRUE)
            if(is.function(monitor))
            { if(!inherits(opt, "try-error"))
              cat("\b | Local search =", 
                  format(opt$value, digits = getOption("digits")))
              else cat("\b |", opt[1])
              cat("\n")
            }
            if(!inherits(opt, "try-error"))
            { Pop[i,] <- opt$par
            Fitness[i] <- opt$value }
            # update object
            object@population <- Pop
            object@fitness <- Fitness
            # update iterations summary
            fitnessSummary[iter,] <- gaSummary(object@fitness)
            object@summary <- fitnessSummary
          }
        }
        
        if(keepBest) 
        { 
          object@bestSol[[iter]] <- unique(Pop[Fitness[1,] == max(Fitness[1,], na.rm = TRUE),, drop=FALSE]) 
        }
        
        # apply a user's defined function to update the GA object
        if(is.function(postFitness))
        { 
          object <- do.call(postFitness, c(object, callArgs))
          Fitness[1,] <- object@fitness
          #Fitness[2,] <- object@fitvar
          Fitness[-c(1),] <- object@fitder
          Pop <- object@population
        }
        
        # update iterations summary
        fitnessSummary[iter,] <- gaSummary(object@fitness)
        object@summary <- fitnessSummary
        
        if(is.function(monitor)) 
        { monitor(object) }
        
        # check stopping criteria
        if(iter > 1)
          object@run <- garun(fitnessSummary[seq(iter),1])
        if(object@run >= run) break  
        if(max(Fitness[1,], na.rm = TRUE) >= maxFitness) break
        if(object@iter == maxiter) break  
        
        ord <- order(Fitness[1,], decreasing = TRUE)
        PopSorted <- Pop[ord,,drop=FALSE]
        FitnessSorted <- Fitness[,ord]
        
        # selection
        if(is.function(selection))
        { 
          sel <- selection(object)
          # sel <- do.call(selection, c(object, callArgs))
          Pop <- sel$population
          Fitness[1,] <- sel$fitness
          #Fitness[2,] <- sel$fitvar
          Fitness[-c(1),] <- sel$fitder
        }
        else
        { sel <- sample(1:popSize, size = popSize, replace = TRUE)
        Pop <- object@population[sel,]
        Fitness[1,] <- object@fitness[sel]
        #Fitness[2,] <- object@fitvar[sel]
        Fitness[-c(1),] <- object@fitder[,sel]
        }
        object@population <- Pop
        object@fitness <- Fitness[1,]
        #object@fitvar <- Fitness[2,]
        object@fitder <- Fitness[-c(1),] 
        
        # crossover
        if(is.function(crossover) & pcrossover > 0)
        { nmating <- floor(popSize/2)
        mating <- matrix(sample(1:(2*nmating), size = (2*nmating)), ncol = 2)
        for(i in seq_len(nmating))
        { if(pcrossover > runif(1))
        { parents <- mating[i,]
        Crossover <- crossover(object, parents)
        Pop[parents,] <- Crossover$children
        Fitness[1,parents] <- Crossover$fitness
        Fitness[2,parents] <- Crossover$fitness
        Fitness[-c(1,2),parents] <- Crossover$fitness
        # print(Crossover$fitness)
        }
        }             
        object@population <- Pop
        object@fitness <- Fitness[1,]
        #object@fitvar <- Fitness[2,]
        object@fitder <- Fitness[-c(1),] 
        }
        
        # mutation
        pm <- if(is.function(pmutation)) pmutation(object) else pmutation
        if(is.function(mutation) & pm > 0)
        { for(i in seq_len(popSize)) 
        { if(pm > runif(1)) 
        { Mutation <- mutation(object, i)
        Pop[i,] <- Mutation
        Fitness[,i] <- rep(NA, nrow(Fitness))
        }
        }
          object@population <- Pop
          object@fitness <- Fitness[1,]
         # object@fitvar <- Fitness[2,]
          object@fitder <- Fitness[-c(1),] 
        }
        
        # elitism
        if(elitism > 0) # (elitism > 0 & iter > 1) 
        { ord <- order(object@fitness, na.last = TRUE)
        u <- which(!duplicated(PopSorted, margin = 1))
        Pop[ord[1:elitism],] <- PopSorted[u[1:elitism],]
        Fitness[,ord[1:elitism]] <- FitnessSorted[,u[1:elitism]]
        object@population <- Pop
        object@fitness <- Fitness[1,]
        #object@fitvar <- Fitness[2,]
        object@fitder <- Fitness[-c(1),] 
        } 
        
      }
      
      # if optim is required perform a local search from the best 
      # solution at the end of GA iterations
      if(optim & (type == "real-valued"))
      { 
        optimArgs$control$maxit <- rev(optimArgs$control$maxit)[1]
        i <- which.max(object@fitness)
        # if not provided suggest approx parscale
        # if(is.null(optimArgs$control$parscale))
        #   optimArgs$control$parscale <- 10^round(log10(abs(object@population[i,])+1))
        # run local search
        opt <- try(suppressWarnings(
          do.call(stats::optim, 
                  c(list(fn = fitness,
                         par = object@population[i,],
                         method = optimArgs$method,
                         lower = optimArgs$lower,
                         upper = optimArgs$upper,
                         control = optimArgs$control), 
                    callArgs))
        ), silent = TRUE)
        if(is.function(monitor))
        { if(!inherits(opt, "try-error"))
          cat("\b | Final local search =",
              format(opt$value, digits = getOption("digits")))
          else cat("\b |", opt[1])
        }
        if(!inherits(opt, "try-error"))
        { object@population[i,] <- opt$par
        object@fitness[i] <- opt$value }
      }
      
      # if(is.function(monitor)) 
      #   { cat("\n"); flush.console() }
      
      # in case of premature convergence remove NAs from summary 
      # fitness evalutations
      object@summary <- na.exclude(object@summary)
      attr(object@summary, "na.action") <- NULL
      
      # get solution(s)
      object@fitnessValue <- max(object@fitness, na.rm = TRUE)
      valueAt <- which(object@fitness == object@fitnessValue)
      solution <- object@population[valueAt,,drop=FALSE]
      if(nrow(solution) > 1)
      { # find unique solutions to precision given by default tolerance
        eps <- gaControl("eps")
        solution <- unique(round(solution/eps)*eps, margin = 1)
      }
      colnames(solution) <- parNames(object)
      object@solution <- solution
      if(keepBest)
        object@bestSol <- object@bestSol[!sapply(object@bestSol, is.null)]  
      
      # return an object of class 'ga'
      return(object)
    }
    
    setClassUnion("numericOrNA", members = c("numeric", "logical"))
    
    setClass(Class = "ga", 
             representation(call = "language",
                            type = "character",
                            lower = "numericOrNA", 
                            upper = "numericOrNA", 
                            nBits = "numericOrNA", 
                            names = "character",
                            popSize = "numeric",
                            iter = "numeric", 
                            run = "numeric", 
                            maxiter = "numeric",
                            suggestions = "matrix",
                            population = "matrix",
                            elitism = "numeric", 
                            pcrossover = "numeric", 
                            pmutation = "numericOrNA",
                            optim = "logical",
                            fitness = "numericOrNA",
                            fitvar = "numericOrNA",
                            fitder = "matrix",
                            summary = "matrix",
                            bestSol = "list",
                            fitnessValue = "numeric",
                            solution = "matrix"
             ),
             package = "GA" 
    ) 
    
    setMethod("print", "ga", function(x, ...) str(x))
    
    setMethod("show", "ga",
              function(object)
              { cat("An object of class \"ga\"\n")
                cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
                cat("Available slots:\n")
                print(slotNames(object))
              }) 
    
    summary.ga <- function(object, ...)
    {
      nvars <- ncol(object@population)
      varnames <- parNames(object)
      domain <- NULL
      if(object@type == "real-valued")
      { domain <- rbind(object@lower, object@upper)
      rownames(domain) <- c("lower", "upper")
      if(ncol(domain) == nvars) 
        colnames(domain) <- varnames
      }
      suggestions <- NULL
      if(nrow(object@suggestions) > 0) 
      { suggestions <- object@suggestions
      dimnames(suggestions) <- list(1:nrow(suggestions), varnames) 
      }
      
      out <- list(type = object@type,
                  popSize = object@popSize,
                  maxiter = object@maxiter,
                  elitism = object@elitism,
                  pcrossover = object@pcrossover,
                  pmutation = object@pmutation,
                  domain = domain,
                  suggestions = suggestions,
                  iter = object@iter,
                  fitness = object@fitnessValue,
                  solution = object@solution)
      class(out) <- "summary.ga"
      return(out)
    }
    
    setMethod("summary", "ga", summary.ga)
    
    print.summary.ga <- function(x, digits = getOption("digits"), ...)
    {
      dotargs <- list(...)
      if(is.null(dotargs$head)) dotargs$head <- 10
      if(is.null(dotargs$tail)) dotargs$tail <- 2
      if(is.null(dotargs$chead)) dotargs$chead <- 10
      if(is.null(dotargs$ctail)) dotargs$ctail <- 2
      
      cat(cli::rule(left = crayon::bold("Genetic Algorithm"), 
                    width = min(getOption("width"),40)), "\n\n")
      # cat("+-----------------------------------+\n")
      # cat("|         Genetic Algorithm         |\n")
      # cat("+-----------------------------------+\n\n")
      
      cat("GA settings: \n")
      cat(paste("Type                  = ", x$type, "\n"))
      cat(paste("Population size       = ", x$popSize, "\n"))
      cat(paste("Number of generations = ", x$maxiter, "\n"))
      cat(paste("Elitism               = ", x$elitism, "\n"))
      cat(paste("Crossover probability = ", format(x$pcrossover, digits = digits), "\n"))
      cat(paste("Mutation probability  = ", format(x$pmutation, digits = digits), "\n"))
      #
      if(x$type == "real-valued")
      { cat(paste("Search domain = \n"))
        do.call(".printShortMatrix", 
                c(list(x$domain, digits = digits), 
                  dotargs[c("head", "tail", "chead", "ctail")]))
      }
      #
      if(!is.null(x$suggestions))
      { cat(paste("Suggestions =", "\n"))
        do.call(".printShortMatrix", 
                c(list(x$suggestions, digits = digits), 
                  dotargs[c("head", "tail", "chead", "ctail")]))
      }
      #
      cat("\nGA results: \n")
      cat(paste("Iterations             =", format(x$iter, digits = digits), "\n"))
      cat(paste("Fitness function value =", format(x$fitness, digits = digits), "\n"))
      if(nrow(x$solution) > 1) 
      { cat(paste("Solutions = \n")) }
      else
      { cat(paste("Solution = \n")) }
      do.call(".printShortMatrix", 
              c(list(x$solution, digits = digits), 
                dotargs[c("head", "tail", "chead", "ctail")]))
      #
      invisible()
    }
    
    
    plot.ga <- function(x, y, ylim, cex.points = 0.7,
                        col = c("green3", "dodgerblue3", adjustcolor("green3", alpha.f = 0.1)),
                        pch = c(16, 1), lty = c(1,2), legend = TRUE,
                        grid = graphics:::grid, ...)
    {
      object <- x  # Argh.  Really want to use 'object' anyway
      is.final <- !(any(is.na(object@summary[,1])))
      iters <- if(is.final) 1:object@iter else 1:object@maxiter
      summary <- object@summary
      if(missing(ylim)) 
      { ylim <- c(max(apply(summary[,c(2,4)], 2, 
                            function(x) min(range(x, na.rm = TRUE, finite = TRUE)))),
                  max(range(summary[,1], na.rm = TRUE, finite = TRUE))) 
      }
      
      plot(iters, summary[,1], type = "n", ylim = ylim, 
           xlab = "Generation", ylab = "Fitness value", ...)
      if(is.final & is.function(grid)) 
      { grid(equilogs=FALSE) }
      points(iters, summary[,1], type = ifelse(is.final, "o", "p"),
             pch = pch[1], lty = lty[1], col = col[1], cex = cex.points)
      points(iters, summary[,2], type = ifelse(is.final, "o", "p"),
             pch = pch[2], lty = lty[2], col = col[2], cex = cex.points)
      if(is.final)
      { polygon(c(iters, rev(iters)), 
                c(summary[,4], rev(summary[,1])), 
                border = FALSE, col = col[3]) }
      else
      { title(paste("Iteration", object@iter), font.main = 1) }
      if(is.final & legend)
      { inc <- !is.na(col)
      legend("bottomright", 
             legend = c("Best", "Mean", "Median")[inc], 
             col = col[inc], pch = c(pch,NA)[inc], 
             lty = c(lty,1)[inc], lwd = c(1,1,10)[inc], 
             pt.cex = c(rep(cex.points,2), 2)[inc], 
             inset = 0.02) }
      
      out <- data.frame(iter = iters, summary)
      invisible(out)
    }
    
    setMethod("plot", "ga", plot.ga)
    
    setGeneric(name = "parNames", 
               def = function(object, ...) { standardGeneric("parNames") }
    )
    
    setMethod("parNames", "ga",
              function(object, ...)
              { 
                names <- object@names
                nvars <- ncol(object@population)
                if(length(names) == 0)
                { names <- paste("x", 1:nvars, sep = "") }
                return(names)
              })
    
    gaSummary <- function(x, ...)
    {
      # compute summary for each step
      x <- na.exclude(as.vector(x))
      q <- fivenum(x)
      c(max = q[5], mean = mean(x), q3 = q[4], median = q[3], q1 = q[2], min = q[1])
    }
    
    #### defining a selection function with elimination process ####
    ga_selectionRS <- function(object){
     if(object@iter>1){
      r <- 2/(object@popSize * (object@popSize - 1))
      q <- 2/object@popSize
      rank <- (object@popSize+1) - rank(object@fitness, ties.method = "min")
      prob <- 1 + q - (rank-1)*r
      
      ## comparing with the best
      b_ind= which.max(object@fitness)
      ## elimination
      # z = combn(1:object@popSize,2)
      # w <- matrix(NA, nrow=object@popSize, ncol=object@popSize)
      w<- array(0, dim=c(object@popSize))
      for( l in 1:object@popSize){
        if(l==b_ind){
          next
        }
      #  w[z[1,l],z[2,l]]= sqrt(var(object@fitder[,z[1,l]]-object@fitder[,z[2,l]])/nrow(object@fitder))
      #  w[z[2,l],z[1,l]]=w[z[1,l],z[2,l]]    
        w[l]=sqrt(var(object@fitder[,b_ind]-object@fitder[,l])/nrow(object@fitder))
      }
      
     # w <- matrix(data=c(rep(1,object@popSize),object@fitvar^2), nrow = object@popSize, ncol=2, byrow = FALSE)%*%matrix(data=c(object@fitvar^2,rep(1,object@popSize)), nrow = 2, ncol=object@popSize, byrow = TRUE)
     #  w <- sqrt(w)
      # fitMat= matrix(data= rep(c(object@fitness),object@popSize), nrow =object@popSize, ncol = object@popSize, byrow = TRUE )
      # diffMat = fitMat-t(fitMat)
      # diffMat[lower.tri(diffMat, diag = TRUE)]<- 0 
      # w=w*qt(0.9,object@popSize-1)#/sqrt(10)
      # fitMat= matrix(data= rep(c(object@fitness),object@popSize), nrow =object@popSize, ncol = object@popSize, byrow = TRUE )
      diffMat = object@fitness[b_ind]-object@fitness
      w=w*qt(0.9,object@popSize-1)#/sqrt(10)
     # w2 <- t(as.matrix(object@fitder))%*%as.matrix(object@fitder)/nrow(object@fitder)
     #  w <- w+w2
      #w[lower.tri(w, diag = TRUE)] <- 0
      # shouldn't it be >0?
     # inferriors= which((diffMat-w>=0) & (diffMat+w>0),arr.ind = TRUE)[,1]
     #  if(length(unique(inferriors))>0.5*object@popSize){
     #    prob[sample(inferriors,0.5*length(inferriors))]<-0
     #  }else{
     #    prob[inferriors]=0}
      #print(unique(inferriors))
      #print(sum(prob>0))
      #print(c(min(object@fitness, na.rm = T),max(object@fitness, na.rm = T)))
      inferriors= which((diffMat-w>0) & (diffMat+w>0))
       if(length(unique(inferriors))>0.5*object@popSize){
         prob[sample(inferriors,0.5*length(inferriors))]<-0
       }else{
         prob[inferriors]=0}
      print(length(inferriors))
      if(max(object@fitness, na.rm = T)>0){
        print("WARNING:")
        print(object@population[which.max(object@fitness),,drop=FALSE])
      }
      prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
      sel <- sample(1:object@popSize, size = object@popSize, 
                    prob = prob, replace = TRUE)
      out <- list(population = object@population[sel,,drop=FALSE],
                  fitness = object@fitness[sel],
                  # fitvar = object@fitvar[sel],
                  fitder = object@fitder[,sel])
      return(out)}else{
        r <- 2/(object@popSize * (object@popSize - 1))
        q <- 2/object@popSize
        rank <- (object@popSize+1) - rank(object@fitness, ties.method = "min")
        prob <- 1 + q - (rank-1)*r

        prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
        sel <- sample(1:object@popSize, size = object@popSize,
                      prob = prob, replace = TRUE)
        out <- list(population = object@population[sel,,drop=FALSE],
                    fitness = object@fitness[sel],
                    # fitvar = object@fitvar[sel],
                    fitder = object@fitder[,sel])
         return(out)
       }
    }
    # ga_postfit <- function(object, ...){
    #   if(object@iter==object@maxiter){
    #     
    #   }else{
    #     return(object)
    #   }
    # }
    ### End Dynamic GA function ####
    #N=matrix(NA,nrow = maxiter_param,ncol = popSize_param)
    ### Defining fitness function for Our GA
    fitness_lm_full <- function(string) {
      #OOB_err<-0
      #beta=10
      # vec<-array(0,dim = samples_param)
      inc <- which(string == 1)
      cc <- colnames(data_corr)[c(inc,ncol(data_corr))]
      # function(data_corr,B1,R,B2,g){
      B1=10
      R=10
      B2=10
      data_tmp= data_corr[,cc]
      ## temporary
      rr=0.5
      rr2=0.5
      q=floor(min(nrow(data_tmp),max(ncol(data_tmp)/((0.632*rr2)^3), floor(nrow(data_tmp)*rr)))) 
      ##
      
      # initializing variables
      tst1<-matrix(0,q,B1)
      tst2<-array(0,dim=c(ceiling(q*rr2),B1,R))
      tst3<-array(0,dim=c(ceiling(rr2^2*q),B1,R,B2))
      tst4<-array(0,dim=c(ceiling(rr2^3*q),B1,R,B2))
      Y1 <- matrix(0,nrow(data_tmp),B1)
      Y2 <- array(0,dim=c(nrow(data_tmp),B1,R))
      Y3 <- array(0,dim=c(nrow(data_tmp),B1,R,B2))
      Y4 <- array(0,dim=c(nrow(data_tmp),B1,R,B2))
      phi1 <- array(0,dim=c(B1))
      phi2 <- array(0,dim=c(B1,R))
      phi3 <- array(0,dim=c(B1,R,B2))
      phi4 <- array(0,dim=c(B1,R,B2))
      
      phi1_ub <- array(0,dim=c(B1))
      phi2_ub <- array(0,dim=c(B1,R))
      S1 <- array(0,dim = c(nrow(data_tmp),B1))
      S2 <- array(0,dim=c(nrow(data_tmp),nrow(data_tmp),B1))
      S3 <- array(0,dim=c(nrow(data_tmp),nrow(data_tmp),B1))
      IF1 <- array(0,dim = c(nrow(data_tmp),B1))
      IF2 <- array(0,dim=c(nrow(data_tmp),nrow(data_tmp),B1))
      S22 <- array(0,dim=c(B1))
      for (j in 1:B1) {
        set.seed(random_seed_train[j])
        tst1[,j]<-sample(1:nrow(data_tmp),q,replace = T)
        for(l in 1:nrow(data_tmp)){
          Y1[l,j]=sum(tst1[,j]==l)
        }
        S1[,j]=q*nrow(data_tmp)*(Y1[,j]/q-1/nrow(data_tmp))
        S2[,,j]=(q*nrow(data_tmp)^2)*(Y1[,j]/q-1/nrow(data_tmp))*t((Y1[,j]/q-1/nrow(data_tmp)))/5
        S3[,,j]=q*nrow(data_tmp)*(-1.2+0.8*q/nrow(data_tmp)+nrow(data_tmp)/5)*(Y1[,j]/q-1/nrow(data_tmp))^2*t((Y1[,j]/q-1/nrow(data_tmp)))/5
        mod1<- lm(y~., data=data_tmp[-tst1[,j],])
        phi1[j]<- sum(Y1[tst1[,j],j]*(predict(mod1,data_tmp[tst1[,j],])-data_tmp[tst1[,j],"y"])^2)/sum(Y1[tst1[,j],j])
        for (o in 1:R){
          set.seed(random_seed_train[B1+o])
          tst2[,j,o]=sample(tst1[,j],ceiling(q*rr2),replace=T)
          for(l in 1:nrow(data_tmp)){
            Y2[l,j,o]=sum(tst2[,j,o]==l)
          }
          tr2=tst1[which(!(tst1[,j] %in% tst2[,j,o])),j]
          mod2<- lm(y~., data=data_tmp[tr2,])
          phi2[j,o] <- sum(Y2[tst2[,j,o],j,o]*(predict(mod2,data_tmp[tst2[,j,o],])-data_tmp[tst2[,j,o],"y"])^2)/sum(Y2[tst2[,j,o],j,o])
          for (u in 1:B2){
            set.seed(random_seed_train[B1+R+u])
            tst3[,j,o,u]=sample(tst2[,j,o],ceiling(rr2^2*q),replace = T)
            tst4[,j,o,u]=sample(tst3[,j,o,u],ceiling(rr2^3*q),replace = T)
            for(l in 1:nrow(data_tmp)){
              Y3[l,j,o,u]=sum(tst3[,j,o,u]==l)
              Y4[l,j,o,u]=sum(tst4[,j,o,u]==l)
            }
            tr3<- tst2[which(!(tst2[,j,o] %in% tst3[,j,o,u])),j,o]
            mod3<- lm(y~., data=data_tmp[tr3,])
            phi3[j,o,u]<- sum(Y3[tst3[,j,o,u],j,o,u]*(predict(mod3,data_tmp[tst3[,j,o,u],])-data_tmp[tst3[,j,o,u],"y"])^2)/sum(Y3[tst3[,j,o,u],j,o,u])
            tr4<- tst2[which(!(tst2[,j,o] %in% tst4[,j,o,u])),j,o]       
            mod4<- lm(y~., data=data_tmp[tr4,])
            phi4[j,o,u]<- sum(Y4[tst4[,j,o,u],j,o,u]*(predict(mod4,data_tmp[tst4[,j,o,u],])-data_tmp[tst4[,j,o,u],"y"])^2)/sum(Y4[tst4[,j,o,u],j,o,u])
            
          }
          ## bias estimation
          bias_tmp=phi3[j,o,]-phi4[j,o,]
          if(sd(bias_tmp)<(2*mean(bias_tmp)*sqrt(B2))){
            phi2_ub[j,o] = phi2[j,o]+mean(bias_tmp)
          }else{
            phi2_ub[j,o]= phi2[j,o]}
          #phi2_ub[j,which(phi2_ub[j,]<0)]=phi2[j,which(phi2_ub[j,]<0)]
        }
        # CV
        cst=R*mean((phi2_ub[j,]-phi2[j,])*(phi2[j,]-phi1[j]))/sum((phi2[j,]-phi1[j])^2)
        phi2_ub[j,]=phi2_ub[j,]+cst*(phi2[j,]-phi1[j])
        IF1[,j]= mean(phi2_ub[j,])*S1[,j]
        S22[j]=sum(S2[,,j]*S2[,,j]-S3[,,j])/((q*nrow(data_tmp))^2)
      }
      
      
      ## computing the bias
      phi2_ub_avg=rowSums(phi2_ub)/R
      # phi2_avg=rowSums(phi2)/R
      cov_mtx=(phi2_ub_avg-mean(phi2_ub_avg))%*%t(S22)
      diag(cov_mtx)=NA
      bias=0.5*sum(cov_mtx,na.rm = T)/(B1*(B1-1))
      ## Computing the variances
      # phi2_red=matrix(0,nrow=R,ncol = B1)
      # phi2_ub_red=matrix(0,nrow=R,ncol = B1)
      # for (id in 1:R) {
      #   # phi2_red[id,]=phi2[,id]-phi2_avg
      #   phi2_ub_red[id,]=phi2_ub[,id]-phi2_ub_avg
      # }
      #var_phi2_Sim=sum((phi2_red)^2)/(R*(R-1)*B1)
      #var_phi2_tot=sum((phi2_avg-mean(phi2_avg))^2)/(B1-1)
      #var_phi2_ub_Sim=sum((phi2_ub_red)^2)/(R*(R-1)*B1)
      var_phi2_ub_tot=sum((phi2_ub_avg-mean(phi2_ub_avg))^2)/(B1-1)
      
      return(-mean(phi2_ub_avg))
      
    }
    fitness_Debiased <- function(string) {
      inc <- which(string == 1)
      cc <- colnames(data_corr)[c(inc,ncol(data_corr))]
      # function(data_corr,B1,R,B2,g){
      B1=10
      R=10
      B2=10
      data_tmp= data_corr[,cc]
      ## temporary
      rr=0.5
      rr2=0.5
      q=floor(min(nrow(data_tmp),max(ncol(data_tmp)/((0.632*rr2)^3), floor(nrow(data_tmp)*rr)))) 
      ##
      
      # initializing variables
      tst1<-matrix(0,q,B1)
      tst2<-array(0,dim=c(ceiling(q*rr2),B1,R))
      tst3<-array(0,dim=c(ceiling(rr2^2*q),B1,R,B2))
      tst4<-array(0,dim=c(ceiling(rr2^3*q),B1,R,B2))
      Y1 <- matrix(0,nrow(data_tmp),B1)
      Y2 <- array(0,dim=c(nrow(data_tmp),B1,R))
      Y3 <- array(0,dim=c(nrow(data_tmp),B1,R,B2))
      Y4 <- array(0,dim=c(nrow(data_tmp),B1,R,B2))
      phi1 <- array(0,dim=c(B1))
      phi2 <- array(0,dim=c(B1,R))
      phi3 <- array(0,dim=c(B1,R,B2))
      phi4 <- array(0,dim=c(B1,R,B2))
      
      phi1_ub <- array(0,dim=c(B1))
      phi2_ub <- array(0,dim=c(B1,R))
      S1 <- array(0,dim = c(nrow(data_tmp),B1))
      S2 <- array(0,dim=c(nrow(data_tmp),nrow(data_tmp),B1))
      S3 <- array(0,dim=c(nrow(data_tmp),nrow(data_tmp),B1))
      IF1 <- array(0,dim = c(nrow(data_tmp),B1))
      IF2 <- array(0,dim=c(nrow(data_tmp),nrow(data_tmp),B1))
      S22 <- array(0,dim=c(B1))
      for (j in 1:B1) {
        set.seed(random_seed_train[j])
        tst1[,j]<-sample(1:nrow(data_tmp),q,replace = T)
        for(l in 1:nrow(data_tmp)){
          Y1[l,j]=sum(tst1[,j]==l)
        }
        S1[,j]=q*nrow(data_tmp)*(Y1[,j]/q-1/nrow(data_tmp))
        S2[,,j]=(q*nrow(data_tmp)^2)*(Y1[,j]/q-1/nrow(data_tmp))*t((Y1[,j]/q-1/nrow(data_tmp)))/5
        S3[,,j]=q*nrow(data_tmp)*(-1.2+0.8*q/nrow(data_tmp)+nrow(data_tmp)/5)*(Y1[,j]/q-1/nrow(data_tmp))^2*t((Y1[,j]/q-1/nrow(data_tmp)))/5
        mod1<- lm(y~., data=data_tmp[-tst1[,j],])
        phi1[j]<- sum(Y1[tst1[,j],j]*(predict(mod1,data_tmp[tst1[,j],])-data_tmp[tst1[,j],"y"])^2)/sum(Y1[tst1[,j],j])
        for (o in 1:R){
          set.seed(random_seed_train[B1+o])
          tst2[,j,o]=sample(tst1[,j],ceiling(q*rr2),replace=T)
          for(l in 1:nrow(data_tmp)){
            Y2[l,j,o]=sum(tst2[,j,o]==l)
          }
          tr2=tst1[which(!(tst1[,j] %in% tst2[,j,o])),j]
          mod2<- lm(y~., data=data_tmp[tr2,])
          phi2[j,o] <- sum(Y2[tst2[,j,o],j,o]*(predict(mod2,data_tmp[tst2[,j,o],])-data_tmp[tst2[,j,o],"y"])^2)/sum(Y2[tst2[,j,o],j,o])
          for (u in 1:B2){
            set.seed(random_seed_train[B1+R+u])
            tst3[,j,o,u]=sample(tst2[,j,o],ceiling(rr2^2*q),replace = T)
            tst4[,j,o,u]=sample(tst3[,j,o,u],ceiling(rr2^3*q),replace = T)
            for(l in 1:nrow(data_tmp)){
              Y3[l,j,o,u]=sum(tst3[,j,o,u]==l)
              Y4[l,j,o,u]=sum(tst4[,j,o,u]==l)
            }
            tr3<- tst2[which(!(tst2[,j,o] %in% tst3[,j,o,u])),j,o]
            mod3<- lm(y~., data=data_tmp[tr3,])
            phi3[j,o,u]<- sum(Y3[tst3[,j,o,u],j,o,u]*(predict(mod3,data_tmp[tst3[,j,o,u],])-data_tmp[tst3[,j,o,u],"y"])^2)/sum(Y3[tst3[,j,o,u],j,o,u])
            tr4<- tst2[which(!(tst2[,j,o] %in% tst4[,j,o,u])),j,o]       
            mod4<- lm(y~., data=data_tmp[tr4,])
            phi4[j,o,u]<- sum(Y4[tst4[,j,o,u],j,o,u]*(predict(mod4,data_tmp[tst4[,j,o,u],])-data_tmp[tst4[,j,o,u],"y"])^2)/sum(Y4[tst4[,j,o,u],j,o,u])
            
          }
          ## bias estimation
          bias_tmp=phi3[j,o,]-phi4[j,o,]
          if(sd(bias_tmp)<(2*mean(bias_tmp)*sqrt(B2))){
            phi2_ub[j,o] = phi2[j,o]+mean(bias_tmp)
          }else{
            phi2_ub[j,o]= phi2[j,o]}
          #phi2_ub[j,which(phi2_ub[j,]<0)]=phi2[j,which(phi2_ub[j,]<0)]
        }
        # CV
        cst=R*mean((phi2_ub[j,]-phi2[j,])*(phi2[j,]-phi1[j]))/sum((phi2[j,]-phi1[j])^2)
        phi2_ub[j,]=phi2_ub[j,]+cst*(phi2[j,]-phi1[j])
        IF1[,j]= mean(phi2_ub[j,])*S1[,j]*S1[,j]/(q*nrow(data_tmp))
        S22[j]=sum(S2[,,j]*S2[,,j]-S3[,,j])/((q*nrow(data_tmp))^2)
      }
      
     # print(mean(phi2),mean(phi2_ub))
      ## computing the bias
      phi2_ub_avg=rowSums(phi2_ub)/R
      # phi2_avg=rowSums(phi2)/R
      cov_mtx=(phi2_ub_avg-mean(phi2_ub_avg))%*%t(S22)
      diag(cov_mtx)=NA
      bias=0.5*sum(cov_mtx,na.rm = T)/(B1*(B1-1))
      ## Computing the variances
      phi2_red=matrix(0,nrow=R,ncol = B1)
      phi2_ub_red=matrix(0,nrow=R,ncol = B1)
      for (id in 1:R) {
        # phi2_red[id,]=phi2[,id]-phi2_avg
        phi2_ub_red[id,]=phi2_ub[,id]-phi2_ub_avg
      }
      #var_phi2_Sim=sum((phi2_red)^2)/(R*(R-1)*B1)
      #var_phi2_tot=sum((phi2_avg-mean(phi2_avg))^2)/(B1-1)
      #var_phi2_ub_Sim=sum((phi2_ub_red)^2)/(R*(R-1)*B1)
      #var_phi2_ub_tot=sum((phi2_ub_avg-mean(phi2_ub_avg))^2)/(B1-1)
      eps_est=phi1-mean(phi1)
      # return(c(-mean(phi2_ub_avg),sqrt(var_phi2_ub_tot-var_phi2_ub_Sim),rowMeans(IF1)/((q*nrow(data_tmp))^2)))
      return(c(-mean(phi2_ub_avg),eps_est+colMeans(IF1)))
    }
    ### Calling the GA function with defined fitness function
    garf <- RSga("binary", fitness = fitness_Debiased , nBits = ncol(data_corr)-1,
                  names = colnames(data_corr)[1:(ncol(data_corr)-1)],
                  selection = ga_selectionRS,
                  maxiter = 10, optim = T, popSize=118,
                  pcrossover = 0.8, pmutation = 0.35,parallel = T, run=10
    )
    # garf <- ga("binary", fitness = fitness_lm_full , nBits = ncol(data_corr)-1,
    #              names = colnames(data_corr)[1:(ncol(data_corr)-1)], 
    #              #selection = ga_selectionRS,
    #              maxiter = 10, optim = T, popSize=50,
    #              pcrossover = 0.5, pmutation = 0.35,parallel = T, run=150
    # )
    plot(garf)
    
    ### Saving the results in selected.cols variable
    #selected.cols <- na.omit(colnames(data_corr)[garf@solution[1,] == 1])
    sorted.pop= garf@population[order(garf@fitness,decreasing = TRUE),]
    # selecting variables that appear in more than 70% of the population
    #selected.cols <- na.omit(colnames(data_corr)[which(colSums(sorted.pop)> (nrow(sorted.pop)*0.75))])
    selected.cols <- na.omit(colnames(data_corr)[which(sorted.pop[1,]==1)])
    ### Saving the result in a csv file
    # write.table(garf@iter, "OPTE_Iter_IdeaH15_cv8.csv", sep = ",",row.names = k ,append = T)
    # write.csv(selected.cols,file="OurGACols.csv")

  
  ### check to make sure none of the results contain variable "y"
  if(any(selected.cols=="y")){
    idx<-which(selected.cols=="y")
    selected.cols<-selected.cols[-idx]
  }
  
  ### Building the model using the rffunc defined above
  model <- lm(y~., data = data_corr[,c(selected.cols,"y")]) 
  
  ## Calculating in sample MSE
  is.MSE<-mean((data_corr$y-predict(model,data_corr[,c(selected.cols,"y")],OOB=TRUE))**2)
  ## Calculating in sample MAE
  # is.Abs<-mean(abs(data_corr$y-predict(model,data_corr[,c(selected.cols,"y")],OOB=TRUE)))
  ### End time
  ptm6<-Sys.time()-ptm
  ## returning the values
  return(c(selected.cols,is.MSE,ptm6))
}
##### End Function ####################

##### A function to calculate Out of Sample error ############  
OOS_err<- function(test_data,selected_subset,k,N){
  N<-10 # number of Macro-replication
  P<-100 # number of post-processing
  B<-30 # maximum number of microreps
  ratio=0.02
  abs_err1=0
  mse_1=0
  r_adj=0
  ## Calculating absolute and squared error with test data set
  for (j in 1:P) {
    set.seed(seeds[1+B+j,k])
    b=sample(1:nrow(test_data),floor(ratio*nrow(test_data)), replace=F)
    test_temp=test_data[b,c(selected_subset,"y")]
    b_test=sample(1:nrow(test_temp),floor(0.5*nrow(test_temp)), replace = F)
    model<- lm(y~.,data = test_temp[-b_test,])
    test.y=test_temp[b_test,c("y")]
    pred=predict(model, test_temp[b_test,],OOB=TRUE)
    
    SS.total      <- sum((test.y - mean(test.y))^2)
    SS.residual   <- sum((test.y - pred)^2)
    #SS.regression <- sum((pred - mean(test.y))^2)
    r_adj=r_adj+(1-SS.residual/SS.total)
    #write.csv(b_test,paste("Postreps_P=",j,"_N=",k,"_n=",n,"_2.csv",sep=""))
    abs_err1=abs_err1+mean(abs(test.y-pred))
    mse_1=mse_1+mean((test.y-pred)^2)
  }
  #abs_err<-abs_err1/P
  MSE<-mse_1/P
  return(c(MSE,r_adj/P))
}
##### End Function ###################################


##### Main ############
N<-10 # number of Macro-replication
P<-100 # number of post-processing
B<-30 # maximum number of microreps

setwd("/home/kvahdat/Data/")
full=fread("./SimDataL.csv")
full=as.data.frame(full)
# full=full[1:1000,c(1:50,ncol(full))]
# colnames(full)[ncol(full)]="y"
## Initializing variables
trueVars=c(1,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,1,rep(0,200))
## Initializing variables
sampleresult<-array(0,dim = c(8,N),
                    dimnames =list(c("OOS MSE","OOS MAE","R^2 Adj","# vars","IS MAE","IS MSE","User time","FalsePos"),c(1:N)) )
colresults<- array(NA,dim = c(dim(full)[2],N),dimnames =list(c(1:dim(full)[2]),c(1:N)))
colresults.bin<- array(0,dim = c(dim(full)[2],N),dimnames =list(c(1:dim(full)[2]),c(1:N)))

set.seed(1)
seeds=matrix(NA, nrow = 1+P+B,ncol = N)

## generating the seeds
for(l in 1:N){
  seeds[,l]=rdunif(n=1+P+B,min=0, max= 10e7)
}

setwd("/home/kvahdat/Data/ch4")
############################################
############ Paralleling Macros ############
############################################
k=1

#wrapper <- function(k,seeds){
  N<-10 # number of Macro-replication
  P<-100 # number of post-processing
  B<-30 # maximum number of microreps
  trueVars=c(1,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,1,rep(0,200))
  ## Initializing variables
  sampleresult<-array(0,dim = c(6,N),
                      dimnames =list(c("OOS MSE","R^2 Adj","# vars","IS MSE","User time","FalsePos"),c(1:N)) )
  colresults<- array(NA,dim = c(dim(full)[2],N),dimnames =list(c(1:dim(full)[2]),c(1:N)))
  colresults.bin<- array(0,dim = c(dim(full)[2],N),dimnames =list(c(1:dim(full)[2]),c(1:N)))
  
  set.seed(seeds[1,k])
  ## Sampling from the combined dataset
  r<-sample(1:nrow(full),floor(nrow(full)*0.5)) # we used 50% for test and 50% for training
  data_corr<-full[r,]
  test_data<-full[-r,]
  OOS<-array(0,dim =6,
              dimnames = list(c("OOS MSE","R^2 Adj","# vars","IS MSE","User time","FalsePos")))
  i=1
  result<-best_subset(data_corr,k) # calling the first function (best subset)
  OOS[4:5]<-round(as.numeric(result[(length(result)-1):length(result)]),digits = 3)
  selectedcols<-result[1:(length(result)-2)]
  colresults[1:length(selectedcols),k]<-selectedcols
  temp2<-which(colnames(data_corr) %in% selectedcols) #indeces
  selectbin<-array(0,dim = dim(data_corr)[2]-1) ### saving the binary coding of the selected variables 
  for (p in 1:length(temp2)) {
    selectbin[temp2[p]]<-1
  }
  colresults.bin[1:length(selectbin),k]<-selectbin
  OOS[6]<-sum(colresults.bin[1:20,k]==trueVars[1:20] & trueVars[1:20]==1)/12
  OOS[1:2]<-round(OOS_err(test_data,selectedcols,k,N),digits = 3)
  OOS[3]<-length(selectedcols)
  write.table(colresults[,k], "SRLMColresults(RSGA_pop118).csv", sep = "," ,append = T)
  write.table(colresults.bin[,k], "SRLMColBin(RSGA_pop118).csv", sep = "," ,append = T)
  write.table(OOS, "SRLMResults(RSGA_pop118).csv", sep = "," ,append = T)
#  return(OOS)
#}

# resss= OOS
#resss=rbind(resss,OOS)
#resss
# final<- foreach(k=1:N, .combine=cbind, .packages = c("caret", "dplyr", "GA", "doParallel","extraDistr")) %dopar% {
#   tempMatrix = wrapper(k,seeds) #calling a function
#   tempMatrix 
# }
# write.csv(final,"./SRResults(B_Res).csv")
# ## End Loop#############
# write.csv(colresults,file = paste("./new/SRLMColresults(B_Res).csv",sep=""))
# write.csv(colresults.bin,file = paste("./new/SRLMColBin(B_Res).csv",sep=""))
# save.image(file = "./SRLMResults(B_Res).RData") 
# 
# 
#wrapper(1,seeds)
