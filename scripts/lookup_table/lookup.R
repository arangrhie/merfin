## GenomeScope2 for Merfin
##
## This is a modified version of Genomescope useful to compute a lookup table of read multiplicities and associated probabilities 

## Install packages if necessary

if (!require("minpack.lm")) install.packages('minpack.lm', repos = "http://cran.us.r-project.org")
library(minpack.lm)

## 0 if NA
###############################################################################

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

## Given mean +/- stderr, report min and max value within 2 SE
###############################################################################

min_max <- function(table){
  ##return (c( abs(table[1]) - 2*abs(table[2]) , abs(table[1])+ 2*abs(table[2])))
  return (c(table[1] - 2*table[2], table[1]+ 2*table[2]))
}

## Use nls to fit 4 peak model
###############################################################################

nls_4peak<-function(x, y, k, estKmercov, estLength, max_iterations){
  model4 = NULL
  best_deviance = Inf
  d_min = 0
  d_initial = 0.10
  d_max = 1
  r_min = 0.00001
  r_initial=0.001
  num_r = 1
  r_max = 1
  kmercov_min = 0
  kmercov_initial = estKmercov
  kmercov_max = Inf
  bias_min = 0
  bias_initial = 0.5
  bias_max = Inf
  length_min = 0
  length_max = Inf
  
  if (VERBOSE) { cat("trying nls_4peak standard algorithm\n") }
  
  try(model4 <- nlsLM(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                           (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length + 
                           (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length + 
                           (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length), 
                    start = list(d=d_initial, r=r_initial, kmercov=estKmercov, bias = bias_initial, length=estLength),
                	lower   = c(c(d_min), rep(r_min, num_r), c(kmercov_min, bias_min, length_min)),
                    upper   = c(c(d_max), rep(r_max, num_r), c(kmercov_max, bias_max, length_max)),
                    control = list(minFactor=1e-12, maxiter=max_iterations, factor=0.1), trace=FALSE), silent = TRUE)
  
  return(model4)
}

## Use nls to fit 2 peak model (haploid case)
###############################################################################

nls_2peak<-function(x, y, k, estKmercov, estLength, max_iterations){
  model2 = NULL
  best_deviance = Inf
  d_min = 0
  d_initial = 0.10
  d_max = 1
  r_min = 0.00001
  r_initial=0.001
  num_r = 1
  r_max = 1
  kmercov_min = 0
  kmercov_initial = estKmercov
  kmercov_max = Inf
  bias_min = 0
  bias_initial = 0.5
  bias_max = Inf
  length_min = 0
  length_max = Inf
    
  if (VERBOSE) { cat("trying nls_2peak standard algorithm\n") }
  
  try(model2 <- nlsLM(y ~ (((1-d)*((1-r)^k))  * dnbinom(x, size = kmercov * 2 / bias, mu = kmercov * 2) * length + 
                        (d*(1-r)^(2*k))    * dnbinom(x, size = kmercov * 4 / bias, mu = kmercov * 4) * length), 

                    start = list(d=d_initial, r=r_initial, kmercov=estKmercov, bias = bias_initial, length=estLength),
                	lower   = c(c(d_min), rep(r_min, num_r), c(kmercov_min, bias_min, length_min)),
                    upper   = c(c(d_max), rep(r_max, num_r), c(kmercov_max, bias_max, length_max)),
                    control = list(minFactor=1e-12, maxiter=max_iterations, factor=0.1), trace=FALSE), silent = TRUE)
  
  return(model2)
}

## score model by number and percent of residual errors after excluding sequencing errors
#########################################################################################

score_model<-function(kmer_hist_orig, nls, round, foldername){
  x = kmer_hist_orig[[1]]
  y = kmer_hist_orig[[2]]
  
  pred=predict(nls, newdata=data.frame(x))
  model_sum=summary(nls)
  kcovfloor = floor(min_max(model_sum$coefficients['kmercov',])[[1]])
  
  ## Compute error rate, by counting kmers unexplained by model through first peak
  ## truncate errors as soon as it goes to zero, dont allow it to go back up
  error_xcutoff = kcovfloor
  error_xcutoff_ind = which(x==error_xcutoff)
  
  error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]
  
  first_zero = -1
  
  for (i in 1:error_xcutoff_ind)
  {
    if (first_zero == -1)
    {
      if (error_kmers[i] < 1.0)
      {
        first_zero = i
        if (VERBOSE) { cat(paste("Truncating errors at", i, "\n")) }
      }
    }
    else
    {
      error_kmers[i] = 0
    }
  }
  
  if (first_zero == -1)
  {
    first_zero = error_xcutoff_ind
  }
  
  ## The fit is residual sum of square error, excluding sequencing errors
  model_fit_all    = c(sum(as.numeric(y[first_zero:length(y)]     - pred[first_zero:length(y)])     ** 2), first_zero, x[length(y)])
  model_fit_full   = c(sum(as.numeric(y[first_zero:(5*kcovfloor)] - pred[first_zero:(5*kcovfloor)]) ** 2), first_zero, (5*kcovfloor))
  model_fit_unique = c(sum(as.numeric(y[first_zero:(3*kcovfloor)] - pred[first_zero:(3*kcovfloor)]) ** 2), first_zero, (3*kcovfloor))
  
  ## The score is the percentage of unexplained kmers, excluding sequencing errors
  model_fit_allscore    = c(1-sum(abs(as.numeric(y[first_zero:length(y)]     - pred[first_zero:length(y)])))     / sum(as.numeric(y[first_zero:length(y)])),     first_zero, x[length(y)])
  model_fit_fullscore   = c(1-sum(abs(as.numeric(y[first_zero:(5*kcovfloor)] - pred[first_zero:(5*kcovfloor)]))) / sum(as.numeric(y[first_zero:(5*kcovfloor)])), first_zero, (5*kcovfloor))
  model_fit_uniquescore = c(1-sum(abs(as.numeric(y[first_zero:(3*kcovfloor)] - pred[first_zero:(3*kcovfloor)]))) / sum(as.numeric(y[first_zero:(3*kcovfloor)])), first_zero, (3*kcovfloor))
  
  fit = data.frame(all  = model_fit_all,      allscore  = model_fit_allscore,
                   full = model_fit_full,     fullscore = model_fit_fullscore, 
                   unique = model_fit_unique, uniquescore = model_fit_uniquescore)
  
  return (fit)
}

## Pick between the two model forms, resolves ambiguity between which is the homozygous and which is the heterozygous peak
###############################################################################

eval_model<-function(kmer_hist_orig, nls1, nls2, round, foldername){
  nls1score = -1
  nls2score = -1
  
  ## Evaluate the score the nls1
  if (!is.null(nls1))
  {
    nls1score = score_model(kmer_hist_orig, nls1, round+0.1, foldername)
    
    if(VERBOSE){ cat(paste("nls1score$all:\t", nls1score$all[[1]], "\n"))}
    
    if (VERBOSE)
    {
      mdir = paste(foldername, "/round", round, ".1", sep="")
      dir.create(mdir, showWarnings=FALSE)
      report_results(kmer_prof_orig,kmer_prof_orig, k, (list(nls1, nls1score)) , mdir, ploidy)
    }
  }
  else
  {
    if (VERBOSE) { cat("nls1score failed to converge\n") }
  }
  
  
  ## Evaluate the score of nls2
  if (!is.null(nls2))
  {
    nls2score = score_model(kmer_hist_orig, nls2, round+0.2, foldername)
    
    if(VERBOSE){ cat(paste("nls2score$all:\t", nls2score$all[[1]], "\n"))}
    
    if (VERBOSE)
    {
      mdir = paste(foldername, "/round", round, ".2", sep="")
      dir.create(mdir, showWarnings=FALSE)
      report_results(kmer_prof_orig, kmer_prof_orig, k, (list(nls2, nls2score)) , mdir, ploidy)
    }
  }
  else
  {
    if (VERBOSE) { cat("nls2score failed to converge\n") }
  }
  
  ## Return the better of the scores
  if (!is.null(nls1))
  {
    if (!is.null(nls2))
    {
      pdiff = abs(nls1score$all[[1]] - nls2score$all[[1]]) / max(nls1score$all[[1]], nls2score$all[[1]])
      
      if (pdiff < SCORE_CLOSE)
      {
        het1 = summary(nls1)$coefficients['r',][[1]]
        het2 = summary(nls2)$coefficients['r',][[1]]
        
        if (het2 * SCORE_HET_FOLD_DIFFERENCE < het1)
        {
          if (VERBOSE) { cat(paste("returning nls1, similar score, higher het\n")) }
          return (list(nls1, nls1score))
        }
        else if (het1 * SCORE_HET_FOLD_DIFFERENCE < het2)
        {
          if (VERBOSE) { cat(paste("returning nls2, similar score, higher het\n")) }
          return (list(nls2, nls2score))
        }
      }
      
      if (nls1score$all[[1]] < nls2score$all[[1]])
      {
        if (VERBOSE) { cat(paste("returning nls1, better score\n")) }
        return (list(nls1, nls1score))
      }
      else
      {
        if (VERBOSE) { cat(paste("returning nls2, better score\n")) }
        return (list(nls2, nls2score))
      }
    }
    else
    {
      if (VERBOSE) { cat(paste("returning nls1, nls2 fail\n")) }
      return (list(nls1, nls1score))
    }
  }
  
  if (VERBOSE) { cat(paste("returning nls2 by default\n")) }
  return (list(nls2, nls2score))
}

## Wrapper function to try fitting 4 peak model with 2 forms
###############################################################################

estimate_Genome_4peak2<-function(kmer_hist_orig, x, y, k, round, foldername, ploidy){

  ## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
  numofKmers   = sum(as.numeric(x)*as.numeric(y))
  estKmercov1  = x[which(y==max(y))][1]
  estLength1   = numofKmers/estKmercov1
  
  if (VERBOSE) { cat(paste("trying with kmercov: ", estKmercov1, "\n")) }

  if (!(ploidy==1)){

    nls1    = nls_4peak(x, y, k, estKmercov1, estLength1, MAX_ITERATIONS)

  } else {

    nls1    = nls_2peak(x, y, k, estKmercov1, estLength1, MAX_ITERATIONS)

  }

  if (VERBOSE) { print(summary(nls1)) }
  
  ## Second we half the max kmercoverage (typically the heterozygous peak)
  estKmercov2  = estKmercov1 / 2 ##2.5
  estLength2   = numofKmers/estKmercov2
  
  if (VERBOSE) { cat(paste("trying with kmercov: ", estKmercov2, "\n")) }
  nls2 = nls_4peak(x, y, k, estKmercov2, estLength2, MAX_ITERATIONS)
  if (VERBOSE) { print(summary(nls2)) }
  
  return(eval_model(kmer_hist_orig, nls1, nls2, round, foldername))
}

## Format numbers
###############################################################################
bp_format<-function(num) {
  paste(formatC(round(num),format="f",big.mark=",", digits=0), "bp",sep=" ")
}

percentage_format<-function(num) {
  paste(signif(num,6)*100,"%",sep="")
}
X_format<-function(num) {
  paste(signif(num,4),"X",sep="")  
}

## Report results and make plots
###############################################################################

report_results<-function(kmer_hist,kmer_hist_orig, k, container, foldername, ploidy)
{
  x=kmer_hist_orig[[1]]
  y=kmer_hist_orig[[2]]
  
  #automatically zoom into the relevant regions of the plot, ignore first 15 positions
  xmax=length(x)
  start=which(y == min(y[1:TYPICAL_ERROR]))
  zoomx=x[start:(xmax-1)]
  zoomy=y[start:(xmax-1)]
  
  ## allow for a little space above max value past the noise
  y_limit = max(zoomy[start:length(zoomy)])*1.1
  
  x_limit = which(y == max(y[start:length(zoomx)])) 
  
  if (min(zoomy) > zoomy[1]){
    x_limit=max(which(zoomy<zoomy[1])[2],600)
  }
  
  if (!is.null(container[[1]]))
  {
    model_sum=summary(container[[1]])
    kcov = min_max(model_sum$coefficients['kmercov',])[1]
    x_limit = max(kcov*5.1, x_limit)
  }
  
  ## Uncomment this to enforce a specific number
  # x_limit=150
  
  ## Features to report
  het=c(-1,-1)
  total_len=c(-1,-1)
  repeat_len=c(-1,-1)
  unique_len=c(-1,-1)
  dups=c(-1,-1)
  error_rate=c(-1,-1)
  model_status="fail"
  
  model_fit_unique      = c(0,0,0)
  model_fit_full        = c(0,0,0)
  model_fit_all         = c(0,0,0)
  model_fit_allscore    = c(0,0,0)
  model_fit_fullscore   = c(0,0,0)
  model_fit_uniquescore = c(0,0,0)
  
  plot_size=7
  font_size=1.2
  resolution=300
  
  ## Plot the distribution, and hopefully with the model fit
  # png(paste(foldername, "/plot.png", sep=""),width=plot_size,height=plot_size, res=resolution)
  pdf(paste(foldername, "/plot.pdf", sep=""),width=plot_size,height=plot_size)
  plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n", xlab="Coverage", ylab="Frequency", ylim=c(0,y_limit), xlim=c(0,x_limit),cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
  points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
  ## if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
  ##    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
  ##}
  box(col="black")
  
  ## Make a second plot in log space over entire range
  # png(paste(foldername, "/plot.log.png", sep=""),width=plot_size,height=plot_size,res=resolution)
  pdf(paste(foldername, "/plot.log.pdf", sep=""),width=plot_size,height=plot_size)
  plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n", xlab="Coverage", ylab="Frequency", log="xy",cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
  rect(1e-10, 1e-10, max(kmer_hist_orig[[1]])*10 , max(kmer_hist_orig[[2]])*10, col=COLOR_BGCOLOR)
  points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
  }
  box(col="black")
  
  if(!is.null(container[[1]])) 
  {
    x=kmer_hist[[1]]
    y=kmer_hist[[2]]
    
    ## The model converged!
    pred=predict(container[[1]], newdata=data.frame(x))
    
    ## Compute the genome characteristics
    model_sum=summary(container[[1]])
    
    ## save the model to a file
    capture.output(model_sum, file=paste(foldername,"/model.txt", sep=""))
    
    ## Identify key values
    het  = min_max(model_sum$coefficients['r',])
    dups = min_max(model_sum$coefficients['bias',])
    kcov = min_max(model_sum$coefficients['kmercov',])
    mlen = min_max(model_sum$coefficients['length',])
    md   = min_max(model_sum$coefficients['d',])
    
    amlen = (mlen[1] + mlen[2]) / 2
    ahet  = (het[1]  + het[2])  / 2
    amd   = (md[1]   + md[2])   / 2
    akcov = (kcov[1] + kcov[2]) / 2
    adups = (dups[1] + dups[2]) / 2
    
    ## Compute error rate, by counting kmers unexplained by model through first peak
    ## truncate errors as soon as it goes to zero, dont allow it to go back up
    error_xcutoff = floor(kcov[1])
    error_xcutoff_ind = which(x==error_xcutoff*2)
    
    error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]
    
    first_zero = -1
    
    for (i in 1:error_xcutoff_ind)
    {
      if (first_zero == -1)
      {
        if (error_kmers[i] < 1.0)
        {
          first_zero = i
          if (VERBOSE) { cat(paste("Truncating errors at", i, "\n")) }
        }
      }
      else
      {
        error_kmers[i] = 0
      }
    }
    
    if (first_zero == -1)
    {
      first_zero = error_xcutoff_ind
    }
    
    ## Rather than "0", set to be some very small number so log-log plot looks okay
    error_kmers = pmax(error_kmers, 1e-10)
    
    total_error_kmers = sum(as.numeric(error_kmers * x[1:error_xcutoff_ind]))
    total_kmers = sum(as.numeric(x*y))
    
    f1 <- function(x){
      i=seq(1,k) 
      h=(1-x)^(k-i)*x^i*choose(k,i)
      sum(h)*total_kmers-total_error_kmers
    }
    
    error_rate_root = try(uniroot(f1, c(0,1))$root)
    
    if (class(error_rate_root) == "try-error")
    {
      error_rate  = c(total_error_kmers/total_kmers/k, total_error_kmers/total_kmers/k)
    }
    else
    {
      error_rate  = c(error_rate_root, error_rate_root)
    }
    
    total_len = (total_kmers-total_error_kmers)/(2*kcov)
    
    ## find kmers that fit the 2 peak model (no repeats)

    if (!(ploidy==1)){

      unique_hist <- (2 * (1 - amd) * (1 - (1 - ahet)^k))                         * dnbinom(x, size = akcov     / adups, mu = akcov)     * amlen +
        ((amd * (1 - (1 - ahet)^k)^2) + (1 - amd) * ((1 - ahet)^k))  * dnbinom(x, size = akcov * 2 / adups, mu = akcov * 2) * amlen 

    } else {

      unique_hist <- ((amd * (1 - (1 - ahet)^k)^2) + (1 - amd) * ((1 - ahet)^k))  * dnbinom(x, size = akcov * 2 / adups, mu = akcov * 2) * amlen   
    
    }


    one_hist <- ((2*(1-amd)*(1-(1-ahet)^k)) + (2*amd*(1-(1-ahet)^k)^2) + (2*amd*((1-ahet)^k)*(1-(1-ahet)^k))) * dnbinom(x, size = akcov   / adups, mu = akcov)   
    two_hist <- (((1-amd)*((1-ahet)^k)) + (amd*(1-(1-ahet)^k)^2))                                             * dnbinom(x, size = akcov*2 / adups, mu = akcov * 2)
    thr_hist <- (2*amd*((1-ahet)^k)*(1-(1-ahet)^k))                                                           * dnbinom(x, size = akcov*3 / adups, mu = akcov * 3)
    fou_hist <- (amd*(1-ahet)^(2*k))                                                                          * dnbinom(x, size = akcov*4 / adups, mu = akcov * 4)
    
    total_kmers = sum(as.numeric(x)*as.numeric(y))
    unique_kmers = sum(as.numeric(x)*unique_hist)
    total_error_kmers = sum(as.numeric(error_kmers * x[1:error_xcutoff_ind]))
    repeat_kmers = total_kmers - unique_kmers - total_error_kmers
    repeat_len=repeat_kmers/(2*kcov)
        
    repeat_len=repeat_kmers/(2*kcov)
    unique_len=unique_kmers/(2*kcov)
    
    score = container[[2]]
    
    model_fit_allscore    = score$allscore
    model_fit_fullscore   = score$fullscore
    model_fit_uniquescore = score$uniquescore
    
    model_fit_all    = score$all
    model_fit_full   = score$full
    model_fit_unique = score$unique
    
    residual = y - pred
    
    ## Finish Log plot
    title(paste("\nlen:",  prettyNum(total_len[1], big.mark=","), 
                "bp", 
                " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                "% ", 
                " het:",  format(100*ahet, digits=3), 
                "%", 
                " kcov:", format(akcov, digits=3), 
                " err:",   format(100*error_rate[1], digits=3),
                "% ", 
                " dup:",  format(adups, digits=3),
                "% ",  
                " k:",   format(k, digits=3), 
                sep=""), 
          cex.main=.85)
    
    ## Mark the modes of the peaks
    akcov
    abline(v=akcov * c(1,2,3,4), col=COLOR_KMERPEAK, lty=2)
    
    ## Draw just the unique portion of the model
    #lines(x, unique_hist, col=COLOR_2PEAK, lty=1, lwd=3)
    #lines(x, pred, col=COLOR_4PEAK, lwd=3)
    lines(x[1:error_xcutoff_ind], error_kmers, lwd=1, col=COLOR_ERRORS)
    
    ## Add legend
    if(length(kmer_hist[,1])==length(kmer_hist_orig[,1])){
      legend(exp(.65 * log(max(x))), 1.0 * max(y),
             
             legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
             lty=c("solid", "solid", "solid", "solid", "dashed"),
             lwd=c(3,3,3,3,3),
             col=c(COLOR_HIST, COLOR_4PEAK, COLOR_2PEAK, COLOR_ERRORS, COLOR_KMERPEAK),
             bg="white")
    }
    else
    {
      legend("topright",
             ##legend(exp(.65 * log(max(x))), 1.0 * max(y),
             legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks","cov-threshold"),
             lty=c("solid", "solid", "solid", "solid", "dashed", "dashed"),
             lwd=c(3,3,3,3,2,3),
             col=c(COLOR_HIST, COLOR_4PEAK, COLOR_2PEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
             bg="white")
    }
    
    dev.set(dev.next())
    
    ## Finish Linear Plot
    title(paste("\nlen:",  prettyNum(total_len[1], big.mark=","), 
                "bp", 
                " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                "% ", 
                " het:",  format(100*ahet, digits=3), 
                "%", 
                " kcov:", format(akcov, digits=3), 
                " err:",   format(100*error_rate[1], digits=3),
                "% ", 
                " dup:",  format(adups, digits=3),
                "% ",  
                " k:",   format(k, digits=3), 
                sep=""), 
          cex.main=.85)
    
    ## Mark the modes of the peaks
    abline(v=akcov * c(1,2,3,4), col=COLOR_KMERPEAK, lty=1)
    
    ## Draw just the unique portion of the model
    lines(x, unique_hist, col=COLOR_2PEAK, lty=1, lwd=3)
    
    #lines(x, pred, col=COLOR_4PEAK, lwd=3)
    lines(x[1:error_xcutoff_ind], error_kmers, lwd=1, col=COLOR_ERRORS)
    
    ## Add legend
    legend(.65 * x_limit, 1.0 * y_limit,
           legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
           lty=c("solid", "solid", "solid", "solid", "dashed"),
           lwd=c(3,3,3,3,2),
           col=c(COLOR_HIST, COLOR_4PEAK, COLOR_2PEAK, COLOR_ERRORS, COLOR_KMERPEAK),
           bg="white")
    
    model_status="done"
    
    cat(paste("Model converged het:", format(ahet, digits=3), 
              " kcov:", format(akcov, digits=3), 
              " err:", format(error_rate[1], digits=3), 
              " model fit:", format(adups, digits=3), 
              " len:", round(total_len[1]), "\n", sep="")) 
  }
  else
  {
    title("\nFailed to converge")
    dev.set(dev.next())
    title("\nFailed to converge")
    cat("Failed to converge")
  }
  
  dev.off()
  dev.off()
  
  ## Write key values to summary file
  summaryFile <- paste(foldername,"/summary.txt",sep="")
  
  format_column_1 = "%-30s"
  format_column_2 = "%-18s"
  format_column_3 = "%-18s"
  
  cat(paste("GenomeScope version 2.0", sep=""),                                                                                                                                                               file=summaryFile, sep="\n") 
  cat(paste("k = ", k,sep=""),                                                                                                                                                                                file=summaryFile, sep="\n", append=TRUE) 
  cat(paste("\n",sprintf(format_column_1,"property"),         sprintf(format_column_2,"min"), sprintf(format_column_3,"max"), sep=""),                                                                        file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Heterozygosity"),        sprintf(format_column_2,percentage_format(het[1])), sprintf(format_column_3,percentage_format(het[2])), sep=""),                                file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Genome Haploid Length"), sprintf(format_column_2,bp_format(total_len[2])), sprintf(format_column_3,bp_format(total_len[1])), sep=""),                                    file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Genome Repeat Length"),  sprintf(format_column_2,bp_format(repeat_len[2])), sprintf(format_column_3,bp_format(repeat_len[1])), sep=""),                                  file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Genome Unique Length"),  sprintf(format_column_2,bp_format(unique_len[2])), sprintf(format_column_3,bp_format(unique_len[1])), sep=""),                                  file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Model Fit "),            sprintf(format_column_2,percentage_format(model_fit_allscore[1])), sprintf(format_column_3,percentage_format(model_fit_fullscore[1])), sep=""), file=summaryFile, sep="\n", append=TRUE)
  cat(paste(sprintf(format_column_1,"Read Error Rate"),       sprintf(format_column_2,percentage_format(error_rate[1])), sprintf(format_column_3,percentage_format(error_rate[2])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
  if (VERBOSE) {
    cat(paste("\nPercent Kmers Modeled (All Kmers) = ",  percentage_format(model_fit_allscore[1]),    " [", model_fit_allscore[2],    ", ", model_fit_allscore[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Percent Kmers Modeled (Full Model) = ",   percentage_format(model_fit_fullscore[1]),   " [", model_fit_fullscore[2],   ", ", model_fit_fullscore[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Percent Kmers Modeled (Unique Kmers) = ", percentage_format(model_fit_uniquescore[1]), " [", model_fit_uniquescore[2], ", ", model_fit_uniquescore[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    
    cat(paste("\nModel RSSE (All Kmers) = ",  model_fit_all[1],    " [", model_fit_all[2],    ", ", model_fit_all[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model RSSE (Full Model) = ",   model_fit_full[1],   " [", model_fit_full[2],   ", ", model_fit_full[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
    cat(paste("Model RSSE (Unique Model) = ", model_fit_unique[1], " [", model_fit_unique[2], ", ", model_fit_unique[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)	
  }
  
  ## Finalize the progress
  progressFilename=paste(foldername,"/progress.txt",sep="")
  cat(model_status, file=progressFilename, sep="\n", append=TRUE)
 
  ## Merfin probabilities

  colors <- c("black","red","green","purple","blue")

  if (!(ploidy==1)){

  	fitted_hist=data.frame(cbind(one_hist,two_hist,thr_hist,fou_hist))
  	
  } else {
  
    fitted_hist=data.frame(cbind(two_hist,fou_hist))
  
  }

  ## Fitted histogram  

  png(paste(foldername, "/fitted_hist.png", sep=""), height = plot_size*300, width = plot_size*300, res=300)
  
  layout(matrix(c(1,2), nrow=2, byrow = TRUE),heights=lcm(c(11,5.5)))
  par(mar=c(0,5,1,1))
  
  plot(kmer_hist_orig, type="n", ylab="Frequency", 
       xlab=NULL, ylim=c(0,y_limit), xlim=c(0,akcov*5),
        cex.main=font_size, cex.sub=font_size, cex.axis=font_size, cex.lab=font_size,
       xaxt='n', yaxt='n')
  myTicks = axTicks(4)
  myTicks[1]<-0
  axis(2, at=myTicks, labels=myTicks, cex.axis=font_size,
       cex.lab=font_size)
  
  peaks<-akcov * 1:6

  points(kmer_hist, type="h", col=COLOR_HIST, lwd=2)
  lines(error_kmers, col="black", lwd=1.5)
  
  abline(v=peaks, col="black", lty=2, lwd=0.3)
  
  for (i in 1:ncol(fitted_hist)) {
    lines(fitted_hist[,i]*amlen, type="l", lwd=1.5, col=colors[i+1])
  }

  lines(rowSums(fitted_hist*amlen), col="darkgray", lwd=2, lty=2)
  
  legend_names=c("0-copy", "1-copy", "2-copy", "3-copy", "4-copy")
  legend_lty=c("solid", "solid", "solid", "solid", "solid")
  legend_lwd=c(3,3,3,3,3)
  
  legend("topright",
		 ##legend(exp(.65 * log(max(x))), 1.0 * max(y),
		 legend=c("Observed",legend_names[1:(ncol(fitted_hist)+1)],"Full model"),
		 lty=c("solid",legend_lty[1:(ncol(fitted_hist)+1)], "dashed"),
		 lwd=c(2, legend_lwd[1:(ncol(fitted_hist)+1)], 3),
		 col=c(COLOR_HIST,colors[1:(ncol(fitted_hist)+1)], "darkgray"),
		 bg="white")

  ## Generate lookup_table
  
  lookup_table <- NULL
  plot_table <- NULL
  
  fitted_hist[(akcov*5):length(one_hist),1:ncol(fitted_hist)] <- 0
  
  fitted_hist <- na.zero(fitted_hist)
  
  for (i in 1:(akcov*5-1)){
    
    totalP<-sum(na.zero(error_kmers[i]), rowSums(fitted_hist[i,]*amlen))
 
    prob<-c(na.zero(error_kmers[i]/totalP),as.numeric(fitted_hist[i,]*amlen)/totalP)
    
    plot_table<-rbind(plot_table,prob)
      
    max.p<-max(prob)  
    readK<-which.max(prob)-1
    
    lookup_table<-rbind(lookup_table,c(readK,max.p))
    
  }
  
  lookup_table<-data.frame(lookup_table)
  rownames(plot_table) <- make.names(plot_table[,1], unique = TRUE)
  
  write.table(lookup_table, paste(foldername, "/lookup_table.txt", sep=""), sep=",", row.names = FALSE, col.names = FALSE)

  ## Plot lookup values
  
  par(mar=c(5,5,0,1))
  
  plot_table<-data.frame(plot_table)
  plot(plot_table$X1,type="n", xlab="Coverage", ylab="Probability",xlim=c(0,akcov*5), ylim=c(0,1), cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size, yaxt='n')

  abline(v=peaks, col="black", lty=2, lwd=0.3)
    
  for (i in 1:(ncol(fitted_hist)+1)) {
    lines(plot_table[,i], type="l", lwd=1.5, col=colors[i])
  }
  axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), cex.axis=font_size,
       cex.lab=font_size)
  
  dev.off()
}

## Main program starts here
###############################################################################

args<-commandArgs(TRUE)

if(length(args) < 3) {
	cat("USAGE: lookup.R histogram_file k-mer_length output_dir [ploidy] [verbose]\n")
	cat("ploidy=1|2\n")
	cat("verbose=0|1\n")
	quit()
} else{

	histfile   <- args[[1]]
	k          <- as.numeric(args[[2]])
	foldername <- args[[3]]

	if (length(args) > 3){
		ploidy <- args[[4]]
	} else {
		ploidy <- 2
	}
}

if (ploidy==1){n=2}else{n=4}

if(length(args) > 4) {
	verbose <- as.numeric(args[[5]])
} else {
	verbose <- as.numeric(0)
}

## Colors for plots
COLOR_BGCOLOR  = "light grey"
COLOR_HIST     = "#56B4E9"
COLOR_4PEAK    = "black"
COLOR_2PEAK    = "#F0E442"
COLOR_ERRORS   = "#D55E00"
COLOR_KMERPEAK = "black"
COLOR_RESIDUAL = "black"
COLOR_COVTHRES = "red"

## Number of rounds before giving up
NUM_ROUNDS=4

## Coverage steps to trim off between rounds
START_SHIFT=5

## Typical cutoff for sequencing error
TYPICAL_ERROR = 15

## Max rounds on NLS
MAX_ITERATIONS=200

## Overrule if two scores are within this percent (0.05 = 5%) but larger difference in het
SCORE_CLOSE = 0.20

## Overrule heterozygosity if there is a large difference in het rate
SCORE_HET_FOLD_DIFFERENCE = 10

maxCovGenomeLen = -1
  
VERBOSE = verbose

dir.create(foldername, showWarnings=FALSE)

kmer_prof <- read.csv(file=histfile,sep="\t", header=FALSE) 

minkmerx = 1;
if (kmer_prof[1,1] == 0) {
	if (VERBOSE) { cat("Histogram starts with zero, reseting minkmerx\n");  }
	minkmerx = 2;
}

kmer_prof <- kmer_prof[c(minkmerx:(length(kmer_prof[,2])-1)),] #get rid of the last position
kmer_prof_orig <- kmer_prof

## try to find the local minimum between errors and the first (heterozygous) peak
start <- which(kmer_prof[, 2] == min(kmer_prof[1:TYPICAL_ERROR, 2]))
maxCovIndex <- length(kmer_prof[, 1])

## terminate after NUM_ROUND iterations, store best result so far in container
round <- 0
best_container <- list(NULL, 0)

while (round < NUM_ROUNDS)
{
  cat(paste("round", round, "trimming to", start, "trying", n, "peak model... "),
      sep = "")
  if (VERBOSE) {
    cat(paste(
      "round",
      round,
      "trimming to",
      start,
      "trying", n, "peak model... \n"
    ))
  }
  
  ## Reset the input trimming off low frequency error kmers
  kmer_prof = kmer_prof_orig[1:maxCovIndex, ]
  x <- kmer_prof[start:maxCovIndex, 1]
  y <- kmer_prof[start:maxCovIndex, 2]
  
  model_4peaks <-
    estimate_Genome_4peak2(kmer_prof, x, y, k, round, foldername, ploidy)
  
  if (!is.null(model_4peaks[[1]])) {
    cat(paste("converged. score: ", model_4peaks[[2]]$all[[1]]), sep = "\n")
    
    if (VERBOSE)
    {
      mdir = paste(foldername, "/round", round, sep = "")
      dir.create(mdir, showWarnings = FALSE)
      report_results(kmer_prof, kmer_prof_orig, k, model_4peaks, mdir, ploidy)
    }
  } else {
    cat(paste("unconverged"))
  }
  
  #check if this result is better than previous
  if (!is.null(model_4peaks[[1]]))
  {
    if (is.null(best_container[[1]]))
    {
      if (VERBOSE) {
        cat(paste("no previous best, updating best"))
      }
      best_container = model_4peaks
    }
    else
    {
      pdiff = abs(model_4peaks[[2]]$all[[1]] - best_container[[2]]$all[[1]]) / max(model_4peaks[[2]]$all[[1]], best_container[[2]]$all[[1]])
      
      if (pdiff < SCORE_CLOSE)
      {
        hetm = summary(model_4peaks[[1]])$coefficients['r', ][[1]]
        hetb = summary(best_container[[1]])$coefficients['r', ][[1]]
        
        if (hetb * SCORE_HET_FOLD_DIFFERENCE < hetm)
        {
          if (VERBOSE) {
            cat(
              paste(
                "model has significantly higher heterozygosity but similar score, overruling"
              )
            )
          }
          best_container = model_4peaks
        }
        else if (hetm * SCORE_HET_FOLD_DIFFERENCE < hetb)
        {
          if (VERBOSE) {
            cat(
              paste(
                "previous best has significantly higher heterozygosity and similar score, keeping"
              )
            )
          }
        }
        else if (model_4peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])
        {
          if (VERBOSE) {
            cat(
              paste(
                "score is marginally better but het rate is not extremely different, upating"
              )
            )
          }
          best_container = model_4peaks
        }
      }
      else if (model_4peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])
      {
        if (VERBOSE) {
          cat(paste("score is significantly better, upating"))
        }
        best_container = model_4peaks
      }
    }
  }
  
  ## Ignore a larger number of kmers as errors
  start <- start + START_SHIFT
  round <- round + 1
}
## Report the results, note using the original full profile
report_results(kmer_prof, kmer_prof_orig, k, best_container, foldername, ploidy)
