#' Plot for a \code{\link{spatialreg.hp}} object
#'
#' @param x A \code{\link{spatialreg.hp}} object.
#' @param  plot.perc Logical;if TRUE, the bar plot (based on ggplot2 package) of the percentage to individual effects of variables and spatial towards total explained variation, the default is FALSE to show plot with original individual effects.
#' @param color Color of variables.
#' @param  commonality Logical; If TRUE, the result of commonality analysis is shown, the default is FALSE. 
 
#' @param  dig Integer; number of decimal places in Venn diagram. 
#' @param ... unused
#' @return a ggplot object
#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}
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

#'plot(spatialreg.hp(COL.lag.eig))
#'plot(spatialreg.hp(COL.lag.eig,commonality=T),commonality=TRUE)



plot.spatialreghp <- function(x, plot.perc = FALSE, commonality=FALSE,color = NULL,dig = 4,...){
  if (!inherits(x, "spatialreghp")){
    stop("x should be the output of spatialreg.hp()")
  }
if(!commonality)
  {if (plot.perc){
    tips3 = data.frame(variable = rownames(x$Individual.R2),
                       value = as.numeric(x$Individual.R2[, "I.perc(%)"]))
    gg = ggplot2::ggplot(tips3, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "% Individual effect to Rsquare (%I)") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))+ggplot2::labs(title ="Individual.R2")
  } else {
    tips2 = data.frame(variable = rownames(x$Individual.R2),
                       value = as.numeric(x$Individual.R2[, "Individual"]))
    gg = ggplot2::ggplot(tips2, ggplot2::aes(x = stats::reorder(variable,-value), y = value)) + ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Variables", y = "Individual effect") + ggplot2::theme(axis.text = element_text(size = 10)) + ggplot2::theme(axis.title = element_text(size = 13))+ggplot2::labs(title ="Individual.R2")
  }
return(gg)
}
if(commonality)
{ 
  Var.part <- as.data.frame(x$commonality.analysis)
  #Var.part <- Var.part[which(Var.part$Fractions >=cutoff), ]
  Var.part$Fractions <- round(Var.part$Fractions,dig)
  variable <- rownames(x$Individual.R2)
  nvar <- length(variable)
  Constrained <- Var.part$Fractions[2^nvar]
  if (!nvar%in% 2:4)
    stop("Venn diagram supports only 2-4 variables")
  else if (nvar == 2)
	Var <- Var.part$Fractions[1:3]
  else if (nvar == 3)
    Var <- Var.part$Fractions[c(1:4, 6, 5, 7)]
  else if (nvar == 4)
    Var <- Var.part$Fractions[c(1:5, 7, 6, 8:10, 12, 11, 14, 13, 15)]
    vegan::showvarparts(part = nvar, bg = color, Xnames = variable, labels = as.character(c(Var, 1-Constrained)))
}

}
