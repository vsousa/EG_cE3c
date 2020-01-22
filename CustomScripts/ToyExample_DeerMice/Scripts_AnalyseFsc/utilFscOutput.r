#' lookForTheBestLhoodSim
#' Function to plot the distribution among runs and produce the file
#' with the BestLhood among the different runs
#' @param fileName name of the file with the lhood and parameter estimates of all runs
#' @param prefix 
#' @param population population tag. dirname is usually population-model
#' @param model model tag. dirname is usually population-model
#' @return maxlhood
lookForTheBestLhoodSim <- function(fileName, prefix, population, model) {
  
  pop <- population
  # read the bestlhood files
  bestlhoods<-read.table(fileName, sep="\t", header=T, stringsAsFactors = F)
  # get the header of the file
  header.tbl<-colnames(bestlhoods)
  # get number of parameters
  col.nb<-1:length(header.tbl)
  
  # check only the runs where the MaxEstLhood was obtained without error
  evaluate <- bestlhoods$MaxEstLhood<0
  
  # get the maxlhood
  maxlhood <- max(bestlhoods$MaxEstLhood[evaluate], na.rm=T)
  
  maxlhood
}

#' log10toln
#' convert log_10 to log_e
#' AIC related functions - Laurent
#' Author: Laurent Excoffier Computation of AIC and relative lhoods
#' INPUT
#' @param  l10 numeric value with the likelihood in log10 scale
log10toln<-function(l10) {
  logn=l10/log10(exp(1))
}

#' AIC
#' compute the AIC given the likelihood in log10 scale and number of parameters
#' AIC related functions - Laurent
#' Author: Laurent Excoffier Computation of AIC and relative lhoods
#' INPUT
#' @param log10L likelihood value in log10 units
#' @param k number of parameters of the model
AIC <- function(log10L, k) {
  lognL=log10toln(log10L)
  2*k-2*lognL
}

#' relLhood
#' compute the relative likelihood of models based on AIC
#' @param AICs vector with the AIC of several models
#' @return vector with relative likelihood of each model
relLhood<-function(AICs) {
  minAIC=min(AICs)
  ws=exp(0.5*(minAIC-AICs))
  ws
}

# COMPUTELHOOD
# computes the loglikelihood given an obs.sfs and exp.sfs
# INPUT
#   obs.sfs : numeric vector with the observed SFS
#   exp.sfs : numeric vector with the expected SFS
# RETURN
#   log likelihood computed as SUM m_i log10(p_i), 
#   where m_i is the ith entry of obs.sfs and p_i is the ith entry of exp.sfs
# NOTE: for entries where m_i > 0 and p_i=0, we replace p_i by a small value (penalty)
computelhood <- function(obs.sfs, exp.sfs) {
  
  lhood <- 0
  
  # remove the first and last entries
  obs.sfs <- obs.sfs[-c(1,length(obs.sfs))]
  exp.sfs <- exp.sfs[-c(1,length(exp.sfs))]
  
  # Get the valid entries, i.e. entries where obs.SFS > 0
  eval <- which(obs.sfs > 0)
  
  # Calculate expected SFS with the penaltie for entries where obs.SFS > 0 and exp.SFS == 0
  if(sum(exp.sfs[eval]==0) > 0) {
    # Settings (penalty for exo SFS entries with zero)
    penalty <- 1e-10
    minpenalty <- 1e-8
    
    penalty <- min(exp.sfs[exp.sfs>0])/100
    if(penalty > minpenalty) {
      penalty <- minpenalty
    } 
    
    # Get the entries which are zero at the obs SFS to have the penalty
    tmp.exp <- exp.sfs[eval]  # note that the length of tmp.exp is length(eval)
    tmp.exp[tmp.exp==0] <- penalty 
    exp.sfs[eval] <- tmp.exp
  }
  
  # check that the sum of exp.sfs is 1
  exp.sfs <- exp.sfs/sum(exp.sfs)
  
  # compute the likelihood
  if(sum(exp.sfs[eval]==0) > 0) {
    print("ERROR: still entries with exp.sfs==0!!!!")
  } else {
    lhood <- sum(obs.sfs[eval]*log10(exp.sfs[eval]))    
  }
  
  lhood
}



# GET_SFSCOORDINATES
# get the coordinates for the x axis of the plots, where 0 is coded as ancestral, 1 is het derived and 2 is homozygous derived
# INPUT:
#   popsinmodel : number of populations in the model
#   pop.sizes   : sample sizes in number of gene copies for each population
# OUTPUT:
#   matrix.sfs : vector of strings with size prod(pop.sizes+1) with the coordinates for each entry of the multidimensional SFS
get_sfscoordinates <- function(popsinmodel, pop.sizes) {
  if (popsinmodel == 7) {
    
    matrix.sfs <-c(rep(0, prod(pop.sizes+1)))
    k<-1
    
    for (i in 0:pop.sizes[1]) {
      
      for (j in 0:pop.sizes[2]) {
        
        for (r in 0:pop.sizes[3]) {
          
          for (x in 0:pop.sizes[4]) {
            
            for (y in 0:pop.sizes[5]) {
              
              for (z in 0:pop.sizes[6]) {
                
                for (v in 0:pop.sizes[7]) {
                  matrix.sfs[k]<-c(paste(i,j,r,x,y,z,v, sep=","))
                  k<-k+1
                }
              }
            }
          }  
        }    
      }
    }
    
  } 
  
  
  if (popsinmodel == 6) {
    
    matrix.sfs <-c(rep(0, prod(pop.sizes+1)))
    k<-1
    
    for (i in 0:pop.sizes[1]) {
      
      for (j in 0:pop.sizes[2]) {
        
        for (r in 0:pop.sizes[3]) {
          
          for (x in 0:pop.sizes[4]) {
            
            for (y in 0:pop.sizes[5]) {
              
              for (z in 0:pop.sizes[6]) {
                
                matrix.sfs[k]<-c(paste(i,j,r,x,y,z, sep=","))
                k<-k+1
              }
            }
          }  
        }    
      }
    }
    
  } 
  if (popsinmodel == 5) {
    
    matrix.sfs <-c(rep(0, prod(pop.sizes+1)))
    k<-1
    
    for (i in 0:pop.sizes[1]) {
      
      for (j in 0:pop.sizes[2]) {
        
        for (r in 0:pop.sizes[3]) {
          
          for (x in 0:pop.sizes[4]) {
            
            for (y in 0:pop.sizes[5]) {
              
              matrix.sfs[k]<-c(paste(i,j,r,x,y, sep=","))
              k<-k+1              
            }
          }  
        }    
      }
    }
  } 
  if(popsinmodel == 4) {
    
    matrix.sfs <-c(rep(0, prod(pop.sizes+1)))
    k<-1
    
    for (i in 0:pop.sizes[1]) {
      
      for (j in 0:pop.sizes[2]) {
        
        for (r in 0:pop.sizes[3]) {
          
          for (x in 0:pop.sizes[4]) {
            
            
            matrix.sfs[k]<-c(paste(i,j,r,x, sep=","))
            k<-k+1
          }
        }
      }  
    }        
  }
  matrix.sfs
}

#This function creates a color scale for use with the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "horiz" argument
#defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
image.scale <- function(z, zlim, col = rainbow(12), breaks, horiz=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){ylim<-c(0,1); xlim<-range(breaks)}
  if(!horiz){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

#' error.bar
#' Function to add error bars - copied from http://monkeysuncle.stanford.edu/?p=485
#' @param x vector with the values of x-axis
#' @param upper vector with y-axis for upper values for corresponding x values
#' @param lower vector with y-axis for lower values for corresponding x values
#' @param length length of the error bar
error.bar <- function(x, upper, lower=upper, length=0.1,...){
  if(length(x) !=length(lower) | length(lower) != length(upper)) {
    stop("vectors must be same length")
  }    
  arrows(x,upper, x, lower, angle=90, code=3, length=length, ...)
}

#' GETINDSLABEL
#' get the individual label as a string with pop index and ind for a pair of selected populations
#' INPUT:
#' @param poplabels : vector of strings of size numpop with the population names
#' @param popindex  : vector of size numpop with pop index in fastsimcoal2 output format, e.g. c(1,3,5) means that pops 2 and 4 were unsampled, and hence we do not have lineages for those
#' @param popindexperlin: vector of size number od lineages with the population label, e.g. c(1, 1, 2, 2) if first two lineages belong to pop 1 and last two lineages belong to pop 2
#' @param indindex     :  vector of size numlineages with index of individuals (tips of the tree), usually 1:number of lineages
#' @param selectedpops :  vector of strings of size two, with the two selected populations e.g. c("YRB", "SAR")
#' RETURN
#' @return  vector of strings with the labels of individuals from the selectedpops, e.g. c("6.5", "9.7") mean lineage 6 from pop5 and lineage 9 from pop 7
getindslabel <- function(poplabels, popindex, popindexperlin, indindex, selectedpops) {
  selectedpopsindex <- c(which(poplabels==selectedpops[1]),which(poplabels==selectedpops[2]))
  
  indslabel <- paste(indindex[popindexperlin==popindex[selectedpopsindex[1]]], 
                     popindexperlin[popindexperlin==popindex[selectedpopsindex[1]]], sep=".")
  for(i in 2:length(selectedpopsindex)) {
    indslabel <- c(indslabel, paste(indindex[popindexperlin==popindex[selectedpopsindex[i]]], 
                                    popindexperlin[popindexperlin==popindex[selectedpopsindex[i]]], sep="."))
  }  
  
  indslabel
}






#--- Function to remove separators within a string
removeTrailingSep=function(string, sep) {
  temp=strsplit(string, split=sep)
  temp2=temp[[1]][nchar(temp[[1]])>0]
  cleanStr=temp2[1]
  if (length(temp2)>1) {
    for (i in 2:length(temp2)) {
      cleanStr=paste(cleanStr, temp2[i], sep=sep)
    }
  }
  cleanStr
}

#--- Replace Keep by -9999
replaceKeep=function(string) {
  if (grepl("keep", string)) {
    return(gsub("keep", "-9999", string))
  }
  return(string)
}

#--- Reading numbers on separate lines -----
getNumbers=function(start, parFile, numSamples) {  
  for (i in 1:numSamples) {
    curnum=as.numeric(unlist(strsplit(parFile[start+i], split=separator))[1])
    if (i==1) {
      num=curnum 
    } else {
      num=c(num, curnum)
    }
  }
  num
}

readSampleSizesTimesAndInbreedingLevel=function(start, parFile, numSamples) {
  for (i in 1:numSamples) {
    curLine=unlist(strsplit(parFile[start+i], split=separator))
    curSampSize=as.numeric(curLine[1])
    curSampTime=0
    curInbreeding=0
    if (length(curLine)>1) curSampTime=as.numeric(curLine[2])
    if (length(curLine)>2) curInbreeding=as.numeric(curLine[3])
    if (i==1) {
      sampSize=curSampSize
      sampTime=curSampTime
      inbreeding=curInbreeding
    } else {
      sampSize=c(sampSize,curSampSize)
      sampTime=c(sampTime,curSampTime)
      inbreeding=c(inbreeding,curInbreeding)
    }
  }
  list(ss=sampSize, st=sampTime, inb=inbreeding)
}

#--- Read migration matrix
readMigMat=function(start, parFile, numSamples) {
  for (i in 1:numSamples) {
    tmpseparator=separator
    if(is.na(as.numeric(unlist(strsplit(parFile[start+i], split=separator)))[1])) {
      tmpseparator=separator2
    }
    curmigs=as.numeric(unlist(strsplit(parFile[start+i], split=tmpseparator)))
    if (i==1) {
      migs=curmigs 
    } else {
      migs=rbind(migs, curmigs)
    }
  }
  rownames(migs)=1:numSamples
  migs 
}


#' GETALLFILESSERVER
#' function to get the ALL file with the parameters for each run from server 
#' INPUT
#' @param settings list with pathtofolder, pop, model info
getALLFilesServer <- function(settings_input) {
  # get the relevant info from settings
  pop <- settings_input$poptag
  model <- settings_input$modeltag
  pathtofolder <- settings_input$pathtofolder
  
  # define the file name with bestlhood from all runs
  file.name<-paste(pop, "-", model, "_ALL.param", sep="")
  
  # Get the file with ALL_params
  system(paste("scp -r ", pathtofolder, "/", file.name, " .", sep=""))
  # Check that file was copied
  if(!file.exists(file.name)) {
    stop(paste("File ", pathtofolder, "/", file.name, " does not exist!",sep=""))
  }
  # move the file to the folder with results
  file.copy(file.name, paste(pop, "-", model ,sep=""))
  file.remove(file.name)
}

#' GETBESTPARAM
#' function to read the file with each run and get the best run parameters and likelihood
#' INPUT
#' @param  settings list with pathtofolder, pop, model, nosingletons
#' OUTPUT:
#' @return  list with the header, bestrun params, run and run2 (indicate the best run)
getbestparam <- function(settings_input) {
  # get relevant info from settings
  pop <- settings_input$poptag
  model <- settings_input$modeltag
  pathtofolder <- settings_input$pathtofolder
  
  file.name<-paste(pop, "-",model, "_ALL.param", sep="") 
  
  # read file with results of each run
  params <- read.table(file.name, header=T)
  # number of parameters is the number of columns - 4 (1st column, Rescaling factor, MaxEstLhood, MaxObsLhood)
  nparam <- ncol(params)-4
  
  #Read the models and generate the boxplot .pdf files
  header <- names(params)
  # clean the header of "X." and "." due to using "$"
  param_i <- which(grepl("X.", header, fixed=TRUE))
  header[param_i] <- unlist(strsplit(header[param_i], "X."))[seq(from=2,to=length(param_i)*2,by=2)]
  
  # Save the number of the run for each bootstrap replicate
  run <- as.numeric(rownames(params))
  run2 <- as.numeric(unlist(strsplit(as.character(params[,1]), split="run"))[seq(from=2,to=nrow(params)*2,by=2)])
  
  # Find the run with the highest likelihood
  maxlhood_i <- which(params$MaxEstLhood==max(params$MaxEstLhood))
  
  # Save the best likelihoods in a variable  
  bestlik <- as.numeric(params[maxlhood_i,2:length(params)])  
  names(bestlik) <- header[-1]
  
  # Output results
  list(header=header[-1], bestlik=bestlik, run=run[maxlhood_i], run2=run2[maxlhood_i])
}


#' GETRUNFILESSERVER
#' function to get the files corresponding to the best run from server (maxL.par, bestlhoods, etc.)
#' INPUT
#'  @param settings_input list with pathtofolder, pop, model, nosingletons
#'  @param results_input  list with results
getRUNFilesServer <- function(settings_input, results_input) {
  
  pop <- settings_input$poptag
  model <- settings_input$modeltag
  pathtofolder <- settings_input$pathtofolder
  file.name<-paste(pop, "-",model, "_ALL.param", sep="")
  
  run <- results_input$run
  run2 <- results_input$run2
  
  # Copy the folder with the best run
  if(sum(run-run2)==0) {
    system(paste("scp -r ", pathtofolder, "/run", run2, " ./", pop, "-", model,sep=""))    
  } else {
    system(paste("scp -r ", pathtofolder, "/run", run2, " .", pop, "-", model,sep=""))    
    warning(paste("Runs missing from ",file.name))
  }
  
  # Check that best run folder was copied
  if(!file.exists( paste( pop, "-", model, "./run",run2, sep="")) ) {
    stop(paste("File with best run: ", pop, "-", model, "./run",run2, " not found!"))  
  }
  
  # # Create folder with maxPar files
  # folder <- paste("./", pop, "-",model,sep="")
  # dir.create(paste(folder, "/maxParfiles_", pop, "-",model, sep=""))
  # for(i in 1:length(run2)) {
  #   # get the maxParFiles for the best run of each replicate
  #   filenameMaxParRun <- paste(folder, "/run",run2[i], "/", pop, "-",model, "/", pop, "-",model, "_maxL.par", sep="")
  #   filenameMaxPar <- paste(folder, "/maxParfiles_", pop, "-",model,"/", pop,"-", model, "_maxL.par", sep="")
  #   if(file.exists(filenameMaxParRun)) {
  #     file.copy(filenameMaxParRun, filenameMaxPar)      
  #   }
  #   # get the Rescaling factor RSF files for the best run of each replicate
  #   filenameRsfRun <- paste(folder, "/run",run2[i], "/", pop, "-", model, "/", pop, "-", model, "_maxL.rsf", sep="")
  #   filenameRsf <- paste(folder, "/maxParfiles_", pop, "-", model,"/", pop, "-", model, "_maxL.rsf", sep="")
  #   if(file.exists(filenameRsfRun)) {
  #     file.copy(filenameRsfRun, filenameRsf)      
  #   }
  # }
}




#' GETBARPLOTWORSTFITSFS
#' function to look at the worst n fitted entries of the joint SFS
#' @param obs.SFS multidimensional array with the obs.SFS
#' @param rel.exp.sfs multidimensional array with the exp.SFS in probabilities
#' @param model strring with the tag for the model, usually poptag-modeltag
#' @param nworst number of worst entries to plot, default=30
#' @return plot with the worst entries fit, and a matrix with the index of worst entries
getBarplotWorstFitSFS <- function(obs.SFS, rel.exp.sfs, model, nworst=30) {
  
  # get number of populations and sample size
  ss <- dim(obs.SFS)
  npop <- length(ss)
  
  # Read and compute observed SFS, discarding entry for fixed ancestral and fixed derived
  index_mono <- matrix(rep.int(1, times=npop), nrow=1)
  obs.SFS[index_mono]<-0  
  
  # get the expected SFS in round(counts)
  exp.sfs <- rel.exp.sfs*sum(obs.SFS)
  
  # compute and plot the relative difference between the SFS
  diff.sfs<-((obs.SFS-exp.sfs)/obs.SFS)
  
  # compute the likelihood based on the expected sfs, but do this for each entry  
  print(paste("Sum of entries of expected SFS==0 is ", sum(rel.exp.sfs==0)))
  numzeroentries_expsfs <- sum(rel.exp.sfs==0)
  
  # to simplify evaluations, convert the multidimensional array to a vector
  obs.SFS_v <- aperm(obs.SFS, c(npop:1))
  obs.SFS_v <- as.vector(obs.SFS_v)
  rel.exp.sfs_v <- aperm(rel.exp.sfs, c(npop:1))
  rel.exp.sfs_v <- as.vector(rel.exp.sfs_v)
  
  eval <- obs.SFS_v > 0 & rel.exp.sfs_v > 0
  exp.log.lik <- log10(rel.exp.sfs_v[eval])*obs.SFS_v[eval]  
  
  rel.obs.sfs <- obs.SFS_v/sum(obs.SFS_v)
  numzeroentries_obssfs <- sum(rel.obs.sfs==0)
  print(paste("Sum of entries of observed SFS==0 is ", numzeroentries_obssfs))
  obs.log.lik <- log10(rel.obs.sfs[eval])*obs.SFS[eval]
  
  # compute the difference in likelihood between expected and obs for each entry
  diff.log.lik <- exp.log.lik-obs.log.lik  
  
  # Get the 30 entries with worst SFS fit in terms of effect on likelihood
  threshold <- as.numeric(quantile(abs(diff.log.lik), (length(diff.log.lik)-nworst)/length(diff.log.lik)))
  
  # get outlier entries  
  outl <- which(diff.log.lik > threshold | diff.log.lik < threshold*(-1))  
  outlierentries <- which(eval)[outl]
  
  # get a matrix with the codes of the SFS entries
  matrix.sfs <- as.matrix(expand.grid(sapply(ss[npop:1], function(x) {1:x}))[,npop:1])
  
  # Read the files with entries with bias
  par(mfrow=c(1,1), mar=c(10,5,2,1), las=2)  
  # get the entries with bias
  entries.SFS <- matrix.sfs[outlierentries,,drop=FALSE]
  
  # plot the diff likelihood for those entries, based on relative SFS and on the observed counts
  rel_expSFS <- exp.sfs[entries.SFS]/sum(exp.sfs)
  rel_obsSFS <- obs.SFS[entries.SFS]/sum(obs.SFS)
  exp.log.lik <- log10(rel_expSFS)*obs.SFS[entries.SFS]
  obs.log.lik <- log10(rel_obsSFS)*obs.SFS[entries.SFS]
  diff.log.lik <- exp.log.lik-obs.log.lik  
  labels_x <- apply(entries.SFS, 1, function(x) {paste(x-1, collapse=",")})
  barplot(diff.log.lik, names.arg=labels_x, ylab="diff. lhood. exp - obs")  
  
  # plot the barplots comparing counts for the expected and observed SFS (note that the expected SFS was multiplied by the number of SNPs)
  res <- matrix(c(exp.sfs[entries.SFS], obs.SFS[entries.SFS]), ncol=nrow(entries.SFS), byrow=T)        
  barx <- barplot(res, beside=T, legend.text=c("exp","obs"), names.arg=labels_x, ylim=c(0,max(res)*1.5), angle=90, main=model)
  
  # plot relative bias in the entries
  relres <- apply(res, 2, function(col){col[1]/col[2]})  
  barx <- barplot(relres, names.arg=labels_x, ylim=c(0,max(relres)*1.15), angle=90, main="Relative fit", ylab="Relative fit=exp/obs")
  abline(h=1)
  
  # return the matrix with the worst entries
  entries.SFS
}


# TABLEPARAMS
# returns a table with the parameters and point estimates
# INPUT
#   results_input : list with header_rescaled, etc.
#   paramstring  : string matching the suffix of a given type of parameters, e.g. "N_" for Ne parameters
# RETURNS
#   table with parameters
tableparams <- function(results_input,paramstring, ndigits, mult) {
  
  header_rescaled <- results_input$header_rescaled 
  header_rescaled <- as.character(sapply(header_rescaled, function(x) {strsplit(x, "[.]")[[1]][2]}))
  
  pIndexes<-as.vector(which(sapply(header_rescaled, function(x) {substr(x,1,nchar(paramstring))})==paramstring))
  
  params <- data.frame(header_rescaled[pIndexes], format(results_input$bestlikConverted[pIndexes], digits = ndigits, scientific=F))
  names(params)<-c("param","point estimate")
  kable(params)  
}

#' computelhoodAIC
#' compute the likelihood and AIC for a SFS with unlinked sites
#'  @param settings_input list with pathtofolder, pop, model, nosingletons
#'  @param results_input  list with results
computelhoodAIC <- function(settings_input, results_input) {
  obsfile <- settings_input$obsfilename_unlinkedSNPs
  popmodel <- paste(settings_input$poptag, "-", settings_input$modeltag, sep="")
  pop.names <- settings_input$pop.names
  minentry <- settings_input$minentry
  
  # read observed SFS and convert to vector
  tmp <- scan(obsfile, skip=1, nlines = 1)
  npop <- tmp[1]  # number of pops
  ss <- tmp[2:(npop+1)]+1 # sample size
  obssfs <- scan(obsfile, skip=2) # read multisfs
  dim(obssfs) <- ss[c(npop:1)]
  obssfs <- aperm(obssfs, c(npop:1))
  
  # read expected SFS and convert to vector
  expsfsfile <- list.files(path=paste("./", popmodel, "./run", results_input$run2, "/", popmodel,sep=""), pattern="*.txt")
  expsfsfile <- paste("./", popmodel, "./run", results_input$run2, "/", popmodel,"/",expsfsfile,sep="")
  if(file.exists(expsfsfile)) {
    expsfs <- scan(expsfsfile, skip=2)
    dim(expsfs) <- ss[c(npop:1)]
    expsfs <- aperm(expsfs, c(npop:1))
  } else {
    stop(paste("Expected SFS file", expsfsfile," not found"))
  }
  
  # check that the obs and exp.sfs have the same size
  if(length(obssfs)!=length(expsfs)) {
    stop("Error: observed and expected SFS have different lengths.")
  }
  
  # compute the likelihood for the dataset with unlinked SNPs
  lhood <- computelhood(obssfs, expsfs, minentry)
  
  # compute the AIC
  nparam <- length(results_input$bestlik)-3
  aic <- AIC(lhood, k=nparam) 
  
  # output a list with the likelihood and AIC
  list(AIC=aic, loglhood=lhood)
}

#' GETFITOBSEXP
#' get the fit of expected SFS to observed SFS
#' @param settings_input list with info about settings and path of files
#' @param results_input list with info about results
#' @return matrix with expected SFS
getfitobsexp <- function(settings_input, results_input) {
  obsfile <- settings_input$obsfilename
  popmodel <- paste(settings_input$poptag, "-", settings_input$modeltag, sep="")
  pop.names <- settings_input$pop.names
  minentry <- settings_input$minentry
  if(settings_input$multiSFS) {
    # read observed SFS
    tmp <- scan(obsfile, skip=1, nlines = 1)
    npop <- tmp[1]  # number of pops
    ss <- tmp[2:(npop+1)]+1 # sample size
    obssfs <- scan(obsfile, skip=2) # read multisfs
    dim(obssfs) <- ss[c(npop:1)]
    obssfs <- aperm(obssfs, c(npop:1))
    
    # read expected SFS
    expsfsfile <- list.files(path=paste("./", popmodel, "./run", results_input$run2, "/", popmodel,sep=""), pattern="*.txt")
    expsfsfile <- paste("./", popmodel, "./run", results_input$run2, "/", popmodel,"/",expsfsfile,sep="")
    if(file.exists(expsfsfile)) {
      expsfs <- scan(expsfsfile, skip=2)
      dim(expsfs) <- ss[c(npop:1)]
      expsfs <- aperm(expsfs, c(npop:1))
    } else {
      stop(paste("Expected SFS file", expsfsfile," not found"))
    }
    
    # check that the obs and exp.sfs have the same size
    if(length(obssfs)!=length(expsfs)) {
      stop("Error: observed and expected SFS have different lengths.")
    }
    
    # plot the entries of the SFS with worst fit
    entries.sfs <- getBarplotWorstFitSFS(obssfs, expsfs, popmodel, nworst=30)
    
    # plot the fit to the marginal 2D SFS
    pairobs <- getmarginal_pairwise_2dsfs(obssfs)
    pairexp <- getmarginal_pairwise_2dsfs(expsfs)
    
    # go through all possible pairwise combinations
    tmp <- sapply(c(1:ncol(pairobs$pairwise)), function(i){
      # plot the 2D marginal SFS fit
      plot2dSFS(pairobs$sfs2d[[i]], pairexp$sfs2d[[i]], 
                xtag=pop.names[pairobs$pairwise[2,i]], 
                ytag=pop.names[pairobs$pairwise[1,i]],
                minentry=minentry);
      # Plot the relative 2D-SFS fit
      plot_relDiff2dSFS(pairobs$sfs2d[[i]], pairexp$sfs2d[[i]], 
                        xtag=pop.names[pairobs$pairwise[2,i]], 
                        ytag=pop.names[pairobs$pairwise[1,i]], 
                        minentry=minentry)
    })
    
    # plot the fit to the marginal 1D SFS
    obs1d <- getmarginal_1dsfs(obssfs)
    exp1d <- getmarginal_1dsfs(expsfs)
    
    par(mfrow=c(1,1))
    tmp <- sapply(c(1:npop), function(i) {
      plot_fitSFS_1d(obs1d[[i]], exp1d[[i]], 
                     c(0,max(obs1d[[i]],exp1d[[i]])*1.1), 
                     pop.names[i], log10scale=TRUE)
    })
  } else if(length(pop.names)>1) {
    # pairwise 2D sfs files
    stop("Function to plot the fit to observed SFS only works for multi-dimensional SFS (--multiSFS option)")
    # TO DO!!!
  } else {
    # pairwise 2D sfs files
    stop("Function to plot the fit to observed SFS only works for multi-dimensional SFS (--multiSFS option)")
    # TO DO!!!
  }
}

# PLOT_FITSFS_1D
# plot the 1d sfs fit to the data
# INPUT
#   obs.sfs : observed SFS
#   exp.sfs : expected SFS
#   ylims   : limits for the y axis
#   pop.name : pop.name
#   log10scale : boolean, if TRUE plot is given in log10 scale
plot_fitSFS_1d <- function(obs.sfs, exp.sfs, ylims, pop.name, log10scale) {
  mymat <- matrix(c(as.numeric(obs.sfs), as.numeric(exp.sfs)*sum(obs.sfs)), ncol=2)
  if(log10scale){
    barplot(t(log10(mymat)), beside = TRUE, names.arg = 0:(length(obs.sfs)-1), legend.text =c("obs","exp"), xlab = pop.name, ylab="log10(count SNPs)", main="Marginal SFS log10 scale", ylim=log10(c(min(mymat[mymat>0])*0.95,max(mymat)*1.05)), xpd = FALSE)
    # lines(0:(length(obs.sfs)-1), log10(exp.sfs*sum(obs.sfs)), type="s", col=3, lty=1, lwd=2)
    # legend("topright", c("obs","exp"), col=c(1,3), lty=c(1,1), lwd=2)
  } else {
    barplot(t(mymat), beside = TRUE, names.arg = 0:(length(obs.sfs)-1), legend.text =c("obs","exp") , xlab = pop.name, ylab="count SNPs", main="Fit SFS")
    # plot(0:(length(obs.sfs)-1), obs.sfs, type="s", ylim=ylims, lwd=2, xlab=pop.name, ylab="SNP counts")
    # lines(0:(length(obs.sfs)-1), exp.sfs*sum(obs.sfs), type="s", col=3, lty=1, lwd=2)
    # legend("topright", c("obs","exp"), col=c(1,3), lty=c(1,1), lwd=2)    
    # plot(0:(length(obs.sfs)-1), log10(obs.sfs), type="s", ylim=log10(ylims), lwd=2, xlab=pop.name, ylab="log10(SNPcounts)")
  }
}

# PLOT 2D SFS
# Plot 2D SFS for observed and expected SFS
# INPUT
#    obsSFS : matrix with observed SFS (counts)
#    expSFS : matrix with expected SFS (probabilities)
#    xtag : string with the label of x-axis
#    ytag : string with the label of y-axis
#    minentry : number with the minimum entry in the SFS (all entries with less than this are pooled together)
# OUTPUT
#    plot with the observed and expected SFS
plot2dSFS <- function(obsSFS,expSFS,xtag,ytag, minentry) {
  layout(matrix(1:3, nrow = 1), widths = c(0.4,0.4,0.2))
  
  # transform data.frames into matrices
  obsSFS <- as.matrix(obsSFS)
  expSFS <- as.matrix(expSFS)
  
  library(RColorBrewer)
  nclasses <- c(0,colorRampPalette(brewer.pal(8,"OrRd"))(15))
  
  # put expected SFS in same scale as obs, by multiplying
  expSFS <- expSFS/sum(expSFS)
  expSFS_counts <- round(expSFS*sum(obsSFS))
  
  # plot the 2D SFS fit
  breaksplot <- c(0, log10(minentry),seq(log10(minentry),max(log10(obsSFS),log10(expSFS_counts)),length.out = length(nclasses)-1))
  image(0:(nrow(obsSFS)-1), 0:(ncol(obsSFS)-1), log10(obsSFS), col=nclasses, breaks=breaksplot, xlab=xtag, ylab=ytag, main="Observed 2D-SFS", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  image(0:(nrow(obsSFS)-1), 0:(ncol(obsSFS)-1), log10(expSFS_counts), col=nclasses, breaks=breaksplot, xlab=xtag, ylab=ytag, main="Expected 2D-SFS", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  image.scale(log10(obsSFS), zlim=range(breaksplot), col = nclasses, breaks=breaksplot, horiz=FALSE, ylab="log10(SFS counts)", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
}

# PLOT_RELDIFF2DSFS
# Plot the fit of the 2D SFS in terms of difference between the observed and expected SFS
# INPUT
#    obsSFS : matrix with observed SFS (counts)
#    expSFS : matrix with expected SFS (probabilities)
#    xtag : string with the label of x-axis
#    ytag : string with the label of y-axis
#    minentry : number with the minimum entry in the SFS (all entries with less than this are pooled together)
# OUTPUT
#    plot with the relative observed/expected SFS
plot_relDiff2dSFS <- function(obsSFS,expSFS,xtag,ytag, minentry) {
  layout(matrix(1:2, nrow = 1), widths = c(0.7,0.3))
  
  library(RColorBrewer)
  
  obsSFS[1,1] <- 0
  expSFS <- expSFS*sum(obsSFS)
  
  eval <- obsSFS > minentry
  diffSFS <- matrix(0,nrow=nrow(obsSFS), ncol=ncol(obsSFS))
  diffSFS[eval] <- (obsSFS[eval]-expSFS[eval])/obsSFS[eval]
  eval <- obsSFS <= minentry
  meandiff <- (sum(obsSFS[eval], na.rm=T)-sum(expSFS[eval], na.rm=T))/sum(obsSFS[eval], na.rm=T)
  diffSFS[eval] <- meandiff/sum(eval)
  
  nclasses <- rev(c(colorRampPalette(brewer.pal(8,"BrBG"))(15)))
  nclasses[8] <- 0
  breaksplot <- c(seq(min(c(-1,diffSFS), na.rm=TRUE),-0.1,length.out = (length(nclasses))/2), seq(0.1,max(c(1,diffSFS), na.rm=TRUE),length.out = (length(nclasses))/2))
  image(0:(nrow(expSFS)-1), 0:(ncol(expSFS)-1), diffSFS, col=nclasses, breaks=breaksplot, xlab=xtag, ylab=ytag, main="Relative difference obs.-exp.", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  image.scale(diffSFS, zlim=range(breaksplot), col = nclasses, breaks=breaksplot, horiz=FALSE, ylab="(obs. SFS - exp. SFS)/obs. SFS", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  
}



# COMPUTELHOOD
# computes the loglikelihood given an obs.sfs and exp.sfs
# INPUT
#   obs.sfs : numeric vector with the observed SFS
#   exp.sfs : numeric vector with the expected SFS
# RETURN
#   log likelihood computed as SUM m_i log10(p_i), 
#   where m_i is the ith entry of obs.sfs and p_i is the ith entry of exp.sfs
# NOTE: for entries where m_i > 0 and p_i=0, we replace p_i by a small value (penalty)
# COMPUTELHOOD
# computes the loglikelihood given an obs.sfs and exp.sfs
# INPUT
#   obs.sfs : numeric vector with the observed SFS
#   exp.sfs : numeric vector with the expected SFS
#   min.entry : numeric value with the minimum entry for the obs SFS. All entries i where obs.sfs[i]<=min.entry are pooled together.
# RETURN
#   log likelihood computed as SUM m_i log10(p_i), 
#   where m_i is the ith entry of obs.sfs and p_i is the ith entry of exp.sfs
# NOTE: for entries where m_i > 0 and p_i=0, we replace p_i by a small value (penalty)
computelhood <- function(obs.sfs, exp.sfs, min.entry) {
  
  lhood <- 0
  
  # remove the first and last entries
  obs.sfs <- obs.sfs[-c(1,length(obs.sfs))]
  exp.sfs <- exp.sfs[-c(1,length(exp.sfs))]
  
  # Get the valid entries, i.e. entries where obs.SFS > 0 
  eval <- which(obs.sfs > 0)
  
  # Calculate expected SFS with the penaltie for entries where obs.SFS > 0 and exp.SFS == 0
  if(sum(exp.sfs[eval]==0) > 0) {
    # Settings (penalty for exp SFS entries with zero)
    penalty <- 1e-10
    minpenalty <- 1e-8
    
    penalty <- min(exp.sfs[exp.sfs>0])/100
    if(penalty > minpenalty) {
      penalty <- minpenalty
    } 
    
    # Get the entries which are zero at the obs SFS to have the penalty
    tmp.exp <- exp.sfs[eval]  # note that the length of tmp.exp is length(eval)
    tmp.exp[tmp.exp==0] <- penalty 
    exp.sfs[eval] <- tmp.exp
  }
  
  # ensure that the sum of exp.sfs is 1
  exp.sfs <- exp.sfs/sum(exp.sfs)
  
  # compute the likelihood
  if(sum(exp.sfs[eval]==0) > 0) {
    stop("ERROR: still entries with exp.sfs==0!!!!")
  } else {
    
    # Pool together all the entries <= min.entry
    indexLargerMinEntry <- which(obs.sfs[eval] > min.entry)
    indexSmallerMinEntry <- which(obs.sfs[eval] <= min.entry)
    if(length(indexLargerMinEntry)>0) {
      lhood <- sum(obs.sfs[eval[indexLargerMinEntry]]*log10(exp.sfs[eval[indexLargerMinEntry]]))  
    } else {
      stop("ERROR: No entries larger than min.entry (-C x option)")
    }
    if(length(indexSmallerMinEntry)>0) {
      lhood <- lhood+sum(obs.sfs[eval[indexSmallerMinEntry]])*log10(sum(exp.sfs[eval[indexSmallerMinEntry]]))  
    } 
  }
  
  lhood
}


#' GETMARGINAL_PAIRWISE_2DSFS
#' obtain the marginal pairwise 2D-SFS given a joint multidimensional SFS.
#' @param jointsfs multidimensional array with the joint SFS
getmarginal_pairwise_2dsfs <- function(jointsfs) {
  # sample size per pop
  samplesize <- dim(jointsfs)
  # number of pops
  npops <- length(samplesize)
  
  # get index of pairwsise combinations
  pairwise_i <- combn(npops,2)
  # get the 2D sfs for each pairwise_i combination
  # using apply to a multidimensional array, summing the marginal of 2 dimensions (given by each column in pairwise_i)
  # where the rows refer to pop1 and cols to pop2
  # hence need to transpose the resulting matrix
  marg2dsfs <- apply(pairwise_i, 2, function(col) {t(apply(jointsfs,col,sum))})
  # output a list with the 2D-SFS and corresponding pairwise comparisons
  list(sfs2d=marg2dsfs, pairwise=pairwise_i)
}

#' GETMARGINAL_1DSFS
#' obtain the marginal 1D-SFS given a joint multidimensional SFS.
#' @param jointsfs multidimensional array with the joint SFS
getmarginal_1dsfs <- function(jointsfs) {
  # sample size per pop
  samplesize <- dim(jointsfs)
  # number of pops
  npops <- length(samplesize)
  # get the 1D sfs for each population
  marg1dsfs <- sapply(1:npops, function(i) {apply(jointsfs,i,sum)})
  # output a list with the 2D-SFS and corresponding pairwise comparisons
  marg1dsfs
}




