#' Hierarchical Partitioning of R2 for Spatial Simultaneous Autoregressive Model

#' @param  mod  Fitted spatialreg objects.
#' @param  iv  optional The relative importance of predictor groups will be assessed. The input for iv should be a list, where each element contains the names of variables belonging to a specific group. These variable names must correspond to the predictor variables defined in the model (mod).

#' @param  commonality Logical; If TRUE, the result of commonality analysis is shown, the default is FALSE. 

#' @details This function conducts hierarchical partitioning to calculate the individual contributions of spatial and each predictor towards total R2 from spatialreg package for spatial simultaneous autoregressive model. 

#' @return \item{Total.R2}{The R2 for the full model.}
#' @return \item{commonality.analysis}{If commonality=TRUE, a matrix containing the value and percentage of all commonality (2^N-1 for N predictors or matrices).}
#' @return \item{Individual.R2}{A matrix containing individual effects and percentage of individual effects for spatial and each predictor}

#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}



#' @references
#' \itemize{
#' \item Lai J.,Zhu W., Cui D.,Mao L.(2023)Extension of the glmm.hp package to Zero-Inflated generalized linear mixed models and multiple regression.Journal of Plant Ecology,16(6):rtad038<DOI:10.1093/jpe/rtad038>
#' \item Lai J.,Zou Y., Zhang S.,Zhang X.,Mao L.(2022)glmm.hp: an R package for computing individual effect of predictors in generalized linear mixed models.Journal of Plant Ecology,15(6):1302-1307<DOI:10.1093/jpe/rtac096>
#' \item Lai J.,Zou Y., Zhang J.,Peres-Neto P.(2022) Generalizing hierarchical and variation partitioning in multiple regression and canonical analyses using the rdacca.hp R package.Methods in Ecology and Evolution,13(4):782-788<DOI:10.1111/2041-210X.13800>
#' \item Chevan, A. & Sutherland, M. (1991). Hierarchical partitioning. American Statistician, 45, 90-96. doi:10.1080/00031305.1991.10475776
#' \item Nimon, K., Oswald, F.L. & Roberts, J.K. (2013). Yhat: Interpreting regression effects. R package version 2.0.0.
#' \item Nimon, Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' }

#'@export
#'@examples
#'library(spatialreg)
#'library(spdep)
#'data(oldcol, package="spdep")
#'listw <- spdep::nb2listw(COL.nb, style="W")
#'ev <- eigenw(listw)
#'W <- as(listw, "CsparseMatrix")
#'trMatc <- trW(W, type="mult")
#'COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, listw=listw,
#'method="eigen", control=list(pre_eig=ev, OrdVsign=1))
#'spatialreg.hp(COL.lag.eig)
#'spatialreg.hp(COL.lag.eig,iv=list(pre1="INC",pre2="HOVAL"))
#'spatialreg.hp(COL.lag.eig,iv=list(pre1="INC",pre2="HOVAL"),commonality=T)

#'COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
#' listw, control=list(pre_eig=ev))
#'spatialreg.hp(COL.errW.eig)
#'spatialreg.hp(COL.errW.eig,iv=list(pre1="INC",pre2="HOVAL"))
#'spatialreg.hp(COL.errW.eig,iv=list(pre1="INC",pre2="HOVAL"),commonality=T)



spatialreg.hp <- function(mod,iv=NULL,commonality=FALSE) 
{
  # initial checks
  if (!inherits(mod, "Sarlm")) stop("spalm.hp only supports Sarlm objects at the moment")
  
  Formu <- strsplit(as.character(mod$call$formula)[3],"")[[1]]
  if("*"%in%Formu)stop("Please put the interaction term as a new variable (i.e. link variables by colon(:)) and  avoid the asterisk (*) and colon(:) in the original model")
  ivname <- strsplit(gsub(" ", "", as.character(mod$call$formula)[3]), "[/+]")[[1]]
  
  
  #if(type=="adjR2")outr2  <- summary(mod)$r.sq
  #if(type=="dev")outr2  <- summary(mod)$dev.expl
  outr2  <- summary(mod,  Nagelkerke=TRUE)$NK
  
  
  dat <- eval(mod$call$data)
  #phy0 <- eval(mod$call$phy)
  
  #if(!inherits(dat, "data.frame")){stop("Please change the name of data object in the original phylolm analysis then try again.")}
  #if('phylolm' %in% class(mod)) to_del <- paste(str_split(as.character(mod$call$formula)[2],'~')[[1]][1],"~","1") else to_del <- paste(as.character(mod$call$formula)[2],"~","1")
  
  
  iv.name <- c("spatial",ivname)
  
  if(is.null(iv)) 
  { 
    nvar <- length(iv.name)
    #  stop("Analysis not conducted. Insufficient number of predictors.")
    
    totalN <- 2^nvar - 1
    binarymx <- matrix(0, nvar, totalN)
    for (i in 1:totalN) {
      binarymx <- creatbin(i, binarymx)
    }
    
    
    
    to_del <- paste(paste("-", ivname, sep= ""), collapse = " ")
    # reduced formula
    modnull<- stats::update(stats::formula(mod), paste(". ~ . ", to_del, sep=""))
    
    mod_null.1 <-  stats::update(object = mod, formula. = modnull, data = dat)
    mod_null.2 <-  lm(modnull,data = dat)
    
    
    
    outputList  <- list()
    outputList[[1]] <- outr2
    
    commonM <- matrix(nrow = totalN, ncol = 3)
    commonM[1, 2]  <- summary(mod_null.1, Nagelkerke=TRUE)$NK  
    for (i in 2:totalN) 
    {
      tmp.name <- iv.name[as.logical(binarymx[, i])]
      if(!'spatial' %in% tmp.name)
      {	
        to_add <- paste("~", paste(c(tmp.name, if (!is.null(mod$offset)) as.character(attr(mod$terms, "variables")[attr(mod$terms, "offset") + 1])), collapse = " + "))
        modnew <- stats::update(object = mod_null.2, data = dat, formula = as.formula(to_add))
        
        #to_add <- paste("~",paste(tmp.name,collapse = " + "),sep=" ")
        # modnew  <- stats::update(object = mod_null, data = dat,to_add) 
        #if(type=="dev")commonM[i, 2]  <- summary(modnew)$dev.expl
        #if(type=="adjR2")commonM[i, 2]  <- summary(modnew)$r.sq
        commonM[i, 2]  <- summary(modnew)$adj.r.squared
      }
      
      if('spatial' %in% tmp.name)
      {
        tmp.name <- tmp.name[-1]
        to_add <- paste("~", paste(c(tmp.name, if (!is.null(mod$offset)) as.character(attr(mod$terms, "variables")[attr(mod$terms, "offset") + 1])), collapse = " + "))
        modnew <- stats::update(object = mod_null.1, data = dat, formula = as.formula(to_add))
        
        #to_add <- paste("~",paste(tmp.name,collapse = " + "),sep=" ")
        # modnew  <- stats::update(object = mod_null, data = dat,to_add) 
        #if(type=="dev")commonM[i, 2]  <- summary(modnew)$dev.expl
        #if(type=="adjR2")commonM[i, 2]  <- summary(modnew)$r.sq
        commonM[i, 2]  <- summary(modnew, Nagelkerke=TRUE)$NK
      }
      
    }
    
    
    
  }
  else
  {
    nvar  <-  length(iv)+1
    ilist <- names(iv)
    if(is.null(ilist))
    {names(iv) <- paste("X",1:(nvar-1),sep="")}
    else
    {whichnoname <- which(ilist=="")
    names(iv)[whichnoname] <- paste("X",whichnoname,sep="")}
    
    ilist <- c("spatial",names(iv))
    
    
    
    ivlist <- ilist
    iv.name <- ilist
    
    ivID <- matrix(nrow = nvar, ncol = 1)
    for (i in 0:nvar - 1) {
      ivID[i + 1]  <-  2^i
    }
    
    totalN  <-  2^nvar - 1
    
    binarymx <- matrix(0, nvar, totalN)
    for (i in 1:totalN) {
      binarymx <- creatbin(i, binarymx)
    }
    
    commonM <- matrix(nrow = totalN, ncol = 3)
    
    
    to_del <- paste(paste("-", ivname, sep= ""), collapse = " ")
    # reduced formula
    modnull<- stats::update(stats::formula(mod), paste(". ~ . ", to_del, sep=""))
    
    mod_null.1 <-  stats::update(object = mod, formula. = modnull, data = dat)
    mod_null.2 <-  lm(modnull,data = dat)
    
    outputList  <- list()
    outputList[[1]] <- outr2
    
    commonM[1, 2]  <- summary(mod_null.1, Nagelkerke=TRUE)$NK    
    for (i in 2:totalN) 
    {
      tmpname <- iv.name[as.logical(binarymx[, i])]
      if(!'spatial' %in% tmpname)
      {
        tmp.name <- unlist(iv[names(iv)%in%tmpname])
        
        to_add <- paste("~", paste(c(tmp.name, if (!is.null(mod$offset)) as.character(attr(mod$terms, "variables")[attr(mod$terms, "offset") + 1])), collapse = " + "))
        modnew <- stats::update(object = mod_null.2, data = dat, formula = as.formula(to_add))
        
        #to_add <- paste("~",paste(tmp.name,collapse = " + "),sep=" ")
        # modnew  <- stats::update(object = mod_null, data = dat,to_add) 
        #if(type=="dev")commonM[i, 2]  <- summary(modnew)$dev.expl
        #if(type=="adjR2")commonM[i, 2]  <- summary(modnew)$r.sq
        commonM[i, 2]  <- summary(modnew)$adj.r.squared
      }
      
      if('spatial' %in% tmpname)
      {
        tmp.name <- unlist(iv[names(iv)%in%tmpname[-1]])
        
        to_add <- paste("~", paste(c(tmp.name, if (!is.null(mod$offset)) as.character(attr(mod$terms, "variables")[attr(mod$terms, "offset") + 1])), collapse = " + "))
        modnew <- stats::update(object = mod_null.1, data = dat, formula = as.formula(to_add))
        
        #to_add <- paste("~",paste(tmp.name,collapse = " + "),sep=" ")
        # modnew  <- stats::update(object = mod_null, data = dat,to_add) 
        #if(type=="dev")commonM[i, 2]  <- summary(modnew)$dev.expl
        #if(type=="adjR2")commonM[i, 2]  <- summary(modnew)$r.sq
        commonM[i, 2]  <- summary(modnew, Nagelkerke=TRUE)$NK
      } 
    }
  }
  
  
  commonlist <- vector("list", totalN)
  
  seqID <- vector()
  for (i in 1:nvar) {
    seqID[i] = 2^(i-1)
  }
  
  
  for (i in 1:totalN) {
    bit <- binarymx[1, i]
    if (bit == 1)
      ivname <- c(0, -seqID[1])
    else ivname <- seqID[1]
    for (j in 2:nvar) {
      bit <- binarymx[j, i]
      if (bit == 1) {
        alist <- ivname
        blist <- genList(ivname, -seqID[j])
        ivname <- c(alist, blist)
      }
      else ivname <- genList(ivname, seqID[j])
    }
    ivname <- ivname * -1
    commonlist[[i]] <- ivname
  }
  
  for (i in 1:totalN) {
    r2list <- unlist(commonlist[i])
    numlist  <-  length(r2list)
    ccsum  <-  0
    for (j in 1:numlist) {
      indexs  <-  r2list[[j]]
      indexu  <-  abs(indexs)
      if (indexu != 0) {
        ccvalue  <-  commonM[indexu, 2]
        if (indexs < 0)
          ccvalue  <-  ccvalue * -1
        ccsum  <-  ccsum + ccvalue
      }
    }
    commonM[i, 3]  <-  ccsum
  }
  
  orderList <- vector("list", totalN)
  index  <-  0
  for (i in 1:nvar) {
    for (j in 1:totalN) {
      nbits  <-  sum(binarymx[, j])
      if (nbits == i) {
        index  <-  index + 1
        commonM[index, 1] <- j
      }
    }
  }
  
  outputcommonM <- matrix(nrow = totalN + 1, ncol = 2)
  totalRSquare <- sum(commonM[, 3])
  for (i in 1:totalN) {
    outputcommonM[i, 1] <- round(commonM[commonM[i,
                                                 1], 3], digits = 4)
    outputcommonM[i, 2] <- round((commonM[commonM[i,
                                                  1], 3]/totalRSquare) * 100, digits = 2)
  }
  outputcommonM[totalN + 1, 1] <- round(totalRSquare,
                                        digits = 4)
  outputcommonM[totalN + 1, 2] <- round(100, digits = 4)
  rowNames  <-  NULL
  for (i in 1:totalN) {
    ii  <-  commonM[i, 1]
    nbits  <-  sum(binarymx[, ii])
    cbits  <-  0
    if (nbits == 1)
      rowName  <-  "Unique to "
    else rowName  <-  "Common to "
    for (j in 1:nvar) {
      if (binarymx[j, ii] == 1) {
        if (nbits == 1)
          rowName  <-  paste(rowName, iv.name[j], sep = "")
        else {
          cbits  <-  cbits + 1
          if (cbits == nbits) {
            rowName  <-  paste(rowName, "and ", sep = "")
            rowName  <-  paste(rowName, iv.name[j], sep = "")
          }
          else {
            rowName  <-  paste(rowName, iv.name[j], sep = "")
            rowName  <-  paste(rowName, ", ", sep = "")
          }
        }
      }
    }
    rowNames  <-  c(rowNames, rowName)
  }
  rowNames  <-  c(rowNames, "Total")
  rowNames <- format.default(rowNames, justify = "left")
  colNames <- format.default(c("Fractions", " % Total"),
                             justify = "right")
  dimnames(outputcommonM) <- list(rowNames, colNames)
  
  VariableImportance <- matrix(nrow = nvar, ncol = 4)
  # VariableImportance <- matrix(nrow = nvar, ncol = 2)
  for (i in 1:nvar) {
    VariableImportance[i, 3] <-  round(sum(binarymx[i, ] * (commonM[,3]/apply(binarymx,2,sum))), digits = 4)
    #VariableImportance[i, 1] <-  round(sum(binarymx[i, ] * (commonM[,3]/apply(binarymx,2,sum))), digits = 4)
  }
  
  VariableImportance[,1] <- outputcommonM[1:nvar,1]
  VariableImportance[,2] <- VariableImportance[,3]-VariableImportance[,1]
  
  total=round(sum(VariableImportance[,3]),digits = 4)
  #total=round(sum(VariableImportance[,1]),digits = 4)
  VariableImportance[, 4] <- round(100*VariableImportance[, 3]/total,2)
  #VariableImportance[, 2] <- round(100*VariableImportance[, 1]/total,2)
  #dimnames(VariableImportance) <- list(iv.name, c("Individual","I.perc(%)"))
  dimnames(VariableImportance) <- list(iv.name, c("Unique","Average.share","Individual","I.perc(%)"))
  
  if(commonality)
  {outputList[[2]]<-outputcommonM
  outputList[[3]]<-VariableImportance
  names(outputList) <- c("Total.R2","commonality.analysis","Individual.R2")
  }
  
  else
  {outputList[[2]]<-VariableImportance
  names(outputList) <- c("Total.R2","Individual.R2")
  }
  
  #outputList[[k+1]]<- c(VariableImportance[,"Individual"],phytree=outr2-sum(VariableImportance[,"Individual"]))
  #outputList[[2]]<-VariableImportance
  
  #if(type=="adjR2"){names(outputList) <- c("adjusted.R2",r2type)}
  #if(type=="dev"){names(outputList) <- c("Explained.deviance",r2type)}
  #if(inherits(mod, "lm")&!inherits(mod, "glm")){names(outputList) <- c("Total.R2",r2type)}
  #outputList$variables <- iv.name
  #if(commonality){outputList$type="commonality.analysis"}
  #if(!commonality){outputList$type="hierarchical.partitioning"}
  
  class(outputList) <- "spatialreghp" # Class definition
  outputList
}



