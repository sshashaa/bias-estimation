###############################################
###############################################
########## RSGA Adaptive Sampling #############
############# Feb 13 2023 #####################
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
best_subset<- function(data_corr,Exp_Name,k){
  require('caret')
  require('GA')
  ## Setting the parameters
  B1=10
  R0=5
  delta=2
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
  train<-matrix(0,nrow = floor(0.5*nrow(data_corr)), ncol = B1)
  for (j in 1:B1) {
    set.seed(random_seed_train[j])
    train[,j]<-sample(1: nrow(data_corr), floor(0.5*nrow(data_corr)),replace = T)
  }
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
    Fitness <- array(NA, dim=c(B1*2+1,popSize))
    sampSize <- array(R0, dim=c(B1,popSize))
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
                  fitvar = Fitness[-c(1:B1+1),],
                  fitdet = Fitness[c(2:(B1+1)),],
                  sampSize = sampSize,
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
      object@fitvar <- Fitness[-c(1:(B1+1)),]
      object@fitdet <- Fitness[c(2:(B1+1)),]
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
        Fitness[1,]<- object@fitness
        Fitness[-c(1:(B1+1)),] <- object@fitvar
        Fitness[c(2:(B1+1)),] <- object@fitdet
        # Fitness[1,] <- object@fitness
        # Fitness[2,] <- object@fitvar
        # Fitness[-c(1),] <- object@fitdet
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
      sampSizeSorted <- sampSize[,ord]
      
      # selection
      if(is.function(selection))
      { 
        sel <- selection(object)
        # sel <- do.call(selection, c(object, callArgs))
        Pop <- sel$population
        Fitness[1,]<- sel$fitness
        Fitness[-c(1:(B1+1)),] <- sel$fitvar
        Fitness[c(2:(B1+1)),] <- sel$fitdet
        sampSize <- sel$sampSize
        #Fitness[1,] <- sel$fitness
        #Fitness[2,] <- sel$fitvar
        #Fitness[-c(1),] <- sel$fitdet
      }
      else
      { sel <- sample(1:popSize, size = popSize, replace = TRUE)
      Pop <- object@population[sel,]
      Fitness[1,]<- object@fitness[sel]
      Fitness[-c(1:(B1+1)),] <- object@fitvar[,sel]
      Fitness[c(2:(B1+1)),] <- object@fitdet[,sel]
      sampSize <- object@sampSize[,sel]
      #Fitness[1,] <- object@fitness[sel]
      #Fitness[2,] <- object@fitvar[sel]
      #Fitness[-c(1),] <- object@fitdet[,sel]
      }
      object@population <- Pop
      object@fitness <- Fitness[1,]
      object@fitvar <- Fitness[-c(1:(B1+1)),]
      object@fitdet <- Fitness[c(2:(B1+1)),]
      
      object@sampSize <- sampSize
      # object@fitness <- Fitness[1,]
      # #object@fitvar <- Fitness[2,]
      # object@fitdet <- Fitness[-c(1),] 
      
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
      Fitness[-c(1:(B1+1)),parents] <- Crossover$fitness
      Fitness[c(2:(B1+1)),parents] <- Crossover$fitness
      sampSize[,parents]<- R0
      # print(Crossover$fitness)
      }
      }             
      object@population <- Pop
      object@fitness <- Fitness[1,]
      object@fitvar <- Fitness[-c(1:(B1+1)),]
      object@fitdet <- Fitness[c(2:(B1+1)),]
      object@sampSize <- sampSize
      #object@fitness <- Fitness[1,]
      #object@fitvar <- Fitness[2,]
      #object@fitdet <- Fitness[-c(1),] 
      }
      
      # mutation
      pm <- if(is.function(pmutation)) pmutation(object) else pmutation
      if(is.function(mutation) & pm > 0)
      { for(i in seq_len(popSize)) 
      { if(pm > runif(1)) 
      { Mutation <- mutation(object, i)
      Pop[i,] <- Mutation
      Fitness[,i] <- rep(NA, nrow(Fitness))
      sampSize[,i]<-R0
      }
      }
        object@population <- Pop
        object@fitness <- Fitness[1,]
        object@fitvar <- Fitness[-c(1:(B1+1)),]
        object@fitdet <- Fitness[c(2:(B1+1)),]
        object@sampSize <- sampSize
        #object@fitness <- Fitness[1,]
        # object@fitvar <- Fitness[2,]
        #object@fitdet <- Fitness[-c(1),] 
      }
      
      # elitism
      if(elitism > 0) # (elitism > 0 & iter > 1) 
      { ord <- order(object@fitness, na.last = TRUE)
      u <- which(!duplicated(PopSorted, margin = 1))
      Pop[ord[1:elitism],] <- PopSorted[u[1:elitism],]
      Fitness[,ord[1:elitism]] <- FitnessSorted[,u[1:elitism]]
      sampSize[,ord[1:elitism]] <- sampSizeSorted[,u[1:elitism]]
      object@population <- Pop
      object@fitness <- Fitness[1,]
      object@fitvar <- Fitness[-c(1:(B1+1)),]
      object@fitdet <- Fitness[c(2:(B1+1)),]
      object@sampSize <- sampSize
      #object@fitness <- Fitness[1,]
      #object@fitvar <- Fitness[2,]
      #object@fitdet <- Fitness[-c(1),] 
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
                          fitness = "numericOrNA", # needs to be a single value
                          fitvar = "matrix",  # can be a vector
                          fitdet = "matrix",       # can be a vector
                          sampSize= "matrix", # tracks the sample size used 
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
    r <- 2/(object@popSize * (object@popSize - 1))
    q <- 2/object@popSize
    
    rank <- (object@popSize+1) - rank(object@fitness, ties.method = "min")
    prob <- 1 + q - (rank-1)*r
    
    prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
    sel <- sample(1:object@popSize, size = object@popSize,
                  prob = prob, replace = TRUE)
    out <- list(population = object@population[sel,,drop=FALSE],
                fitness = object@fitness[sel],
                fitvar = object@fitvar[,sel],
                fitdet = object@fitdet[,sel],
                sampSize = object@sampSize[,sel])
    return(out)
  }
  ga_postfit <- function(object, ...){
    # this condition checks the optimality gap
    
    # find the critical solutions
    sorted_fit= (object@popSize+1) -rank(object@fitness, ties.method = "min")
    best_id= sorted_fit[1]
    # print(best_id)
    if(object@iter>1){
      best_sol=object@summary[object@iter-1,1]
    }else{
    best_sol=object@fitness[best_id]}
    # print((abs((best_sol-mean(object@fitness[-best_id]))/mean(object@fitness[-best_id]))<0.5) & (sqrt(mean(object@fitvar[,best_id]))/abs(object@fitness[best_id])>2))
  #  if((abs((best_sol-mean(object@fitness[-best_id]))/mean(object@fitness[-best_id]))<0.5) & (sqrt(abs(mean(object@fitvar[,best_id])))/abs(best_sol)>2)){#sqrt(var(object@fitness))/abs(mean(object@fitness))<2.5){
      #print(abs((best_sol-mean(object@fitness[-best_id]))/mean(object@fitness[-best_id])))
      # selecting the worst case in all other cases.
      opt_gap = abs((best_sol-mean(object@fitness))/mean(object@fitness))
      est_err = sqrt(abs(mean(object@fitvar[,best_id]))/mean(object@sampSize[,best_id]))
      # print(paste0("The opt gap is: ", round(opt_gap,2), " the est err is: ", round(est_err,2)))
      while((max(colSums(object@sampSize))<60) & (est_err>(abs(best_sol)*opt_gap/2))){
        print(paste0("The opt gap is: ", round(opt_gap,2), " the est err is: ", round(est_err,2)))
        sorted_fit= (object@popSize+1) -rank(object@fitness, ties.method = "min")
        best_id= sorted_fit[1]
        rl= apply(object@fitdet,2,which.min)

        # need a list of sample sizes used for each input model and each variance.
        # Refer to Gao et. al. (2017) for more info
        A1_tst= sum(object@sampSize[,best_id]^2/object@fitvar[,best_id])
        A2_tst=0
        for(l in 1:object@popSize){
      
          if(l==best_id)
            A2_tst=A2_tst+0
          else
            A2_tst=A2_tst+((object@sampSize[rl[l],l])^2)/object@fitvar[rl[l],l]
         
        }
        R_b1=array(NA, dim=c(B1))
        R_b2=array(NA, dim=c(object@popSize))
        for(h in 1:B1){
          des_tmp_id=which.min(object@fitdet[h,-best_id])
          R_b1[h]=(object@fitdet[h,best_id]-object@fitdet[h,-best_id][des_tmp_id])^2/(object@fitvar[h,best_id]/object@sampSize[h,best_id]+object@fitvar[h,-best_id][des_tmp_id]/object@sampSize[h,-best_id][des_tmp_id])
        }
        for(l in 1:object@popSize){
          R_b2[l]=(object@fitdet[rl[l],best_id]-object@fitdet[rl[l],l])^2/(object@fitvar[rl[l],best_id]/object@sampSize[rl[l],best_id]+object@fitvar[rl[l],l]/object@sampSize[rl[l],l])
          if(l==best_id){
            R_b2[l]=10e10
          }
        }
        if(A1_tst<A2_tst){
          # assign extra budget to minimum value of R for scenario of the best design
          inp_ind= which.min(R_b1)
          des_ind= best_id
        }else{
          
          des_ind=c(1:object@popSize)[which.min(R_b2)]
          inp_ind= rl[des_ind]
        }
        # Call the fitness again for those solutions
        new_out=fitness_Debiased_adap(object@population[des_ind,],delta,inp_ind)
        # update the object
        object@fitdet[inp_ind,des_ind]=(new_out[1]*delta+object@fitdet[inp_ind,des_ind]*object@sampSize[inp_ind,des_ind])/(object@sampSize[inp_ind,des_ind]+delta)
        object@fitness[des_ind]=mean(object@fitdet[,des_ind])
        object@fitvar[inp_ind,des_ind]= (object@fitvar[inp_ind,des_ind]*object@sampSize[inp_ind,des_ind]+delta*new_out[2])/(object@sampSize[inp_ind,des_ind]+delta)^2
        object@sampSize[inp_ind,des_ind]=object@sampSize[inp_ind,des_ind]+delta
        opt_gap = abs((best_sol-mean(object@fitness))/mean(object@fitness))
        est_err = sqrt(abs(mean(object@fitvar[,des_ind]))/mean(object@sampSize[,des_ind]))}
      return(object)
    # }else{
    #   return(object)
    # }
  }
  ### End Dynamic GA function ####
  
  fitness_Debiased_fixed <- function(string) {
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
    #var_phi2_ub_Sim=sum((phi2_ub_red)^2)/(R*(R-1))
    #var_phi2_ub_tot=sum((phi2_ub_avg-mean(phi2_ub_avg))^2)/(B1-1)
    #eps_est=phi1-mean(phi1)
    var_phi2_Sim=rowSums((phi2_ub_red)^2)/(R*(R-1))
    # return(c(-mean(phi2_ub_avg),sqrt(var_phi2_ub_tot-var_phi2_ub_Sim),rowMeans(IF1)/((q*nrow(data_tmp))^2)))
    return(c(-mean(phi2_ub_avg),-phi2_ub_avg,var_phi2_Sim))
  }
  
  fitness_Debiased_adap <- function(string,delta=1,b1_ind) {
    inc <- which(string == 1)
    cc <- colnames(data_corr)[c(inc,ncol(data_corr))]
    # function(data_corr,B1,R,B2,g){
    B1=10
    R=delta
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
    for(j in c(b1_ind)){
      #set.seed(random_seed_train[j])
      #tst1[,j]<-sample(1:nrow(data_tmp),q,replace = T)
      tst1[,j]=train[,j]
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
    phi2_ub_avg=sum(phi2_ub[j,])/R
    # phi2_avg=rowSums(phi2)/R
    #cov_mtx=(phi2_ub_avg-mean(phi2_ub_avg))%*%t(S22)
    # diag(cov_mtx)=NA
    # bias=0.5*sum(cov_mtx,na.rm = T)/(B1*(B1-1))
    ## Computing the variances
    phi2_red=matrix(0,nrow=R,ncol = B1)
    phi2_ub_red=matrix(0,nrow=R,ncol = B1)
    for (id in 1:R) {
      # phi2_red[id,]=phi2[,id]-phi2_avg
      phi2_ub_red[id,]=phi2_ub[,id]-phi2_ub_avg
    }
    var_phi2_Sim=sum((phi2_ub_red)^2)/(R*(R-1))
    #var_phi2_tot=sum((phi2_avg-mean(phi2_avg))^2)/(B1-1)
    #var_phi2_ub_Sim=sum((phi2_ub_red)^2)/(R*(R-1)*B1)
    #var_phi2_ub_tot=sum((phi2_ub_avg-mean(phi2_ub_avg))^2)/(B1-1)
    #eps_est=phi1-mean(phi1)
    # return(c(-mean(phi2_ub_avg),sqrt(var_phi2_ub_tot-var_phi2_ub_Sim),rowMeans(IF1)/((q*nrow(data_tmp))^2)))
    return(c(-mean(phi2_ub_avg),var_phi2_Sim))
  }
  ### Calling the GA function with defined fitness function
  garf <- RSga("binary", fitness = fitness_Debiased_fixed , nBits = ncol(data_corr)-1,
               names = colnames(data_corr)[1:(ncol(data_corr)-1)],
               selection = ga_selectionRS, # postFitness = ga_postfit,
               maxiter = 100, optim = T, popSize=floor(ncol(data_corr)-1),
               pcrossover = 0.8, pmutation = 0.35,parallel = T, run=20
  )
  # garf <- ga("binary", fitness = fitness_lm_full , nBits = ncol(data_corr)-1,
  #              names = colnames(data_corr)[1:(ncol(data_corr)-1)], 
  #              #selection = ga_selectionRS,
  #              maxiter = 10, optim = T, popSize=50,
  #              pcrossover = 0.5, pmutation = 0.35,parallel = T, run=150
  # )
  # plot(garf)
  
  ### Saving the results in selected.cols variable
  #selected.cols <- na.omit(colnames(data_corr)[garf@solution[1,] == 1])
  #sorted.pop= garf@population[order(garf@fitness,decreasing = TRUE),]
  # selecting variables that appear in more than 70% of the population
  #selected.cols <- na.omit(colnames(data_corr)[which(colSums(sorted.pop)> (nrow(sorted.pop)*0.75))])
  selected.cols <- na.omit(colnames(data_corr)[-ncol(data_corr)][which(garf@solution[1,]==1)])
  saveRDS(garf,paste0("/home/kvahdat/Data/ch4/GA_",Exp_Name,"_",k,".rds"))
  ### End time
  ptm6<-Sys.time()-ptm
  ## returning the values
  results1= c(as.numeric(ptm6),selected.cols)
  return(results1)
}
##### End Function ####################


##### Main ############
N<-10 # number of Macro-replication
P<-100 # number of post-processing
B<-30 # maximum number of microreps

setwd("/home/kvahdat/Data/")
full=fread("./SimDataL.csv")
full=as.data.frame(full)
## Initializing variables
trueVars=c(1,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,1,rep(0,200))

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
#k=1

#wrapper <- function(k,seeds){
N<-10 # number of Macro-replication
P<-100 # number of post-processing
B<-30 # maximum number of microreps
trueVars=c(1,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,0,1,1,rep(0,200))
trueVars_ind=which(trueVars==1)
c1<- c(1:25,110:120)#sample(c(1:200)[-trueVars_ind],20))
## Initializing variables
sampleresult<-array(0,dim = c(3,N),
                    dimnames =list(c("User time","TruePos","Num"),c(1:N)) )
colresults<- array(NA,dim = c(ncol(full),N),dimnames =list(c(1:ncol(full)),c(1:N)))
colresults.bin<- array(0,dim = c(ncol(full),N),dimnames =list(c(1:ncol(full)),c(1:N)))

# Exp_Name="EstErrx2GTBstSolxOptGap"
# Done on Mar 12 
Exp_Name="NoOCBA_Debiased"
for(k in c(1:N)){
  print(paste0("Iteration ",k))
  set.seed(seeds[1,k])
  ## Sampling from the combined dataset
  r<-sample(1:nrow(full),floor(nrow(full)*0.2)) 

  data_corr<-full[r,c(c1,ncol(full))]
  write.csv(data_corr,paste0("InpDat",Exp_Name,"_",k,".csv"))
  # test_data<-full[-r,]
  results1<- best_subset(data_corr,Exp_Name,k)
  sampleresult[1,k]<- results1[1]
  selectedcols<- results1[-1]
  colresults[1:length(selectedcols),k]<-selectedcols
  sampleresult[2,k]<-sum(colnames(full)[trueVars_ind] %in% selectedcols)/12
  sampleresult[3,k]<-length(selectedcols)
  }

write.csv(colresults, paste0("RSGA_Col_",Exp_Name,".csv") )
write.csv(sampleresult, paste0("RSGA_Res_",Exp_Name,".csv"))
#  return(OOS)
# 
# for(k in c(1:N)){
# garf=readRDS(paste0("/home/kvahdat/Data/ch4/GA_",Exp_Name,"_",k,".rds"))
# selectedcols=na.omit(colnames(data_corr)[-ncol(data_corr)][garf@solution==1])
# colresults[1:length(selectedcols),k]<-selectedcols
# sampleresult[2,k]<-sum(colnames(full)[trueVars_ind] %in% selectedcols)/12
# sampleresult[3,k]<-length(selectedcols)
# }
# print(sampleresult)
# 
# rowMeans(sampleresult[,-c(8,9,10)])
