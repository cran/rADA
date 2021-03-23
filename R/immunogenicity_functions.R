utils::globalVariables(c(".data", "value", "Category","..density.."))

#' Define ImmunoAssay class
#'
#' This stores the data that is used for screening cut point analysis.
#'
#' @slot data Imported data as is, used for CV analysis
#' @slot melted.data Data used for most functions
#' @slot exp.name Experiment name
#' @slot stats List of statistics, results gathered from both coefficient of variation analysis as well as plot generation
#' @slot outlier.rm Has any outlier analysis been performed on this dataset?
#' @slot outlier.rm.method If outlier removal has been performed, what method was used?
#' @slot scp.table Table of cut point information
#' @slot cv.table Table derived from coefficient of variation analysis
#'
#' @rdname ImmunoAssay-class
#'
#' @exportClass ImmunoAssay

ImmunoAssay <- setClass("ImmunoAssay", slots = list(data = "data.frame",
                                                    melted.data = "data.frame", exp.name = "character", stats = "list",
                                                    cv.table = "data.frame",
                                                    outlier.rm = "character", outlier.rm.method = "character", scp.table = "data.frame"))


#' @title Import assay as ImmunoAssay object
#'
#' @description Function to import assay information into an ImmunoAssay object for analysis.
#'
#' @param assay.df Pathname to (.csv or .xlsx files accepted) or imported data.frame consisting of the following columns: 'ID','Lot', and columns identifying the Day, Operator and Replicate like so: 'D1_Op2_3' to indicate Day 1, operator 2, replicate 3.
#' @param exp.name Experiment name. If stays NULL, will automatically name experiment as 'experiment1'.
#'
#' @author Emma Gail
#'
#' @importFrom methods new
#'
#' @return An object of the class ImmunoAssay
#'
#' @export
#'
#' @examples
#'
#' assay.df <- importAssay(assay.df = lognormAssay)

importAssay <- function(assay.df, exp.name = NULL) {

  # read in assay if it is a character
  if (typeof(assay.df) == "character")
  {
    if (endsWith(assay.df, "csv")) {
      assay.df <- utils::read.csv(assay.df)
    } else if (endsWith(assay.df, "xlsx")) {
      assay.df <- openxlsx::read.xlsx(assay.df)
    } else {
      warning("Unknown file type")
      return(NULL)
    }  #end else to if(endsWith(assay.df, 'csv')){
  }  # end if(typeof(assay.df) == 'character'){

  if (is.null(exp.name)) {
    cat("No name given, called experiment1")
    exp.name <- "experiment1"
  }
  assay.df.melted <- assayMelt(assay.df, exp.name = exp.name)

  assay.obj <- methods::new('ImmunoAssay', data = assay.df, melted.data = assay.df.melted,
                           exp.name = exp.name, outlier.rm = "no")

  return(assay.obj)

}  # end importAssay <- function(assay.df){


#' @title Calculate Coefficient of Variation
#'
#' @description The function calculates the mean, standard deviation and coefficient of variation for replicates of an immunoassay.
#'
#' @param assay.obj An ImmunoAssay object imported by importAssay
#' @param cv.threshold Threshold for re-calculation of means and standard deviation based on coefficient of variation. The default threshold is 20 (i.e., 20\% CV)
#'
#' @author Emma Gail
#'
#' @return An object of the class ImmunoAssay with calculated CV stats in the \code{cv.table} slot
#'
#' @export
#'
#' @name calcCvStats
#' @aliases calcCvStats,ImmunoAssay-method
#' @docType methods
#'
#' @examples
#'
#' assay.obj <- importAssay(lognormAssay)
#' assay.obj <- calcCvStats(assay.obj)


setGeneric("calcCvStats", function(assay.obj, cv.threshold = 20) {
  standardGeneric("calcCvStats")
})

setMethod("calcCvStats", "ImmunoAssay", function(assay.obj, cv.threshold = 20) {

  if ("cv" %in% assay.obj@outlier.rm.method) {
    warning("CV analysis has already been performed on this ImmunoAssay, returning original dataset")
    return(assay.obj)
  }

  assay.df <- stats::na.omit(assay.obj@data)
  assay.df.copy <- assay.df
  cv.df <- data.frame(matrix(nrow = nrow(assay.df)))
  data.colnames <- colnames(assay.df)[!(colnames(assay.df) %in% c("ID",
                                                                  "Lot"))]
  us_count <- max(unique(stringr::str_count(data.colnames, "_")))
  if (us_count == 1) {
    colnames.df <- t(data.frame(strsplit(data.colnames, "_")))
  } else if (us_count == 2) {
    colnames.df <- t(data.frame(strsplit(sub("_", "", data.colnames),
                                         "_")))
  }
  means.df <- data.frame(matrix(nrow = nrow(assay.df), ncol = length(unique(colnames.df[,
                                                                                        1]))))
  colnames(means.df) <- unique(colnames.df[, 1])

  sd.df <- data.frame(matrix(nrow = nrow(assay.df), ncol = length(unique(colnames.df[,
                                                                                     1]))))
  colnames(sd.df) <- unique(colnames.df[, 1])

  cv.df <- data.frame(matrix(nrow = nrow(assay.df), ncol = length(unique(colnames.df[,
                                                                                     1]))))
  colnames(cv.df) <- unique(colnames.df[, 1])

  means.adj.df <- data.frame(matrix(nrow = nrow(assay.df), ncol = length(unique(colnames.df[,
                                                                                            1]))))
  colnames(means.adj.df) <- unique(colnames.df[, 1])

  sd.adj.df <- data.frame(matrix(nrow = nrow(assay.df), ncol = length(unique(colnames.df[,
                                                                                         1]))))
  colnames(sd.adj.df) <- unique(colnames.df[, 1])

  cv.adj.df <- data.frame(matrix(nrow = nrow(assay.df), ncol = length(unique(colnames.df[,
                                                                                         1]))))
  colnames(cv.adj.df) <- unique(colnames.df[, 1])

  is.adj.df <- data.frame(matrix(nrow = nrow(assay.df), ncol = length(unique(colnames.df[,
                                                                                         1]))))
  colnames(is.adj.df) <- unique(colnames.df[, 1])


  for (ucolname in unique(colnames.df[, 1])) {
    grouped.colnames <- data.colnames[colnames.df[, 1] == ucolname]
    grouped.df <- assay.df[, grouped.colnames]

    if (FALSE %in% apply(data.frame(grouped.df), 2, is.numeric)) {
      # if there is a FALSE then it is not numeric
      next
    }

    row.means <- rowMeans(data.matrix(grouped.df))
    row.sds <- matrixStats::rowSds(data.matrix(grouped.df))
    row.cv <- (row.sds/row.means) * 100

    gdf.high.cv <- grouped.df[row.cv >= cv.threshold, ]

    high.cv.row.means <- rowMeans(data.matrix(gdf.high.cv))
    high.cv.row.sds <- matrixStats::rowSds(data.matrix(gdf.high.cv))
    gdf.high.cv2 <- gdf.high.cv

    gdf.high.cv2[] <- lapply(gdf.high.cv2, function(x) ifelse((x >
                                                                 (high.cv.row.means + high.cv.row.sds)), NA, x))
    gdf.high.cv2[] <- lapply(gdf.high.cv2, function(x) ifelse((x <
                                                                 (high.cv.row.means - high.cv.row.sds)), NA, x))

    high.cv.row.means2 <- rowMeans(data.matrix(gdf.high.cv2), na.rm = TRUE)
    high.cv.row.sds2 <- matrixStats::rowSds(data.matrix(gdf.high.cv2),
                                            na.rm = TRUE)
    high.cv.new.cv <- (high.cv.row.sds2/high.cv.row.means2) * 100

    means.df[[ucolname]] <- row.means
    sd.df[[ucolname]] <- row.sds
    cv.df[[ucolname]] <- row.cv

    means.adj.df[[ucolname]] <- rep(NA, length(row.means))
    sd.adj.df[[ucolname]] <- rep(NA, length(row.means))
    cv.adj.df[[ucolname]] <- rep(NA, length(row.means))
    is.adj.df[[ucolname]] <- rep(NA, length(row.means))

    means.adj.df[as.numeric(rownames(gdf.high.cv2)), ucolname] <- high.cv.row.means2
    sd.adj.df[as.numeric(rownames(gdf.high.cv2)), paste0(ucolname)] <- high.cv.row.sds2
    cv.adj.df[as.numeric(rownames(gdf.high.cv2)), paste0(ucolname)] <- high.cv.new.cv
    is.adj.df[as.numeric(rownames(gdf.high.cv2)), paste0(ucolname)] <- rep(TRUE,
                                                                           length(as.numeric(rownames(gdf.high.cv2))))

    assay.df.copy[rownames(gdf.high.cv2), colnames(gdf.high.cv2)] <- gdf.high.cv2

  }  #end for(ucolname in unique(colnames.df[,1])){

  assay.obj@stats <- list(means = means.df, sd = sd.df, cv.df = cv.df,
                          means.adj.df = means.adj.df, sd.adj.df = sd.adj.df, cv.adj.df = cv.adj.df,
                          is.adj.df = is.adj.df)



  #assay.obj@data <- assay.df.copy
  #assay.obj@melted.data <- assayMelt(assay.df.copy, exp.name = assay.obj@exp.name)

  assay.df.melt <- assay.obj@melted.data
  assay.df.dcast <- reshape2::dcast(assay.df.melt, DayOperator + Day + Operator + ID + Lot ~ Replicate)
  assay.df.dcast$xid <- paste0(assay.df.dcast$DayOperator,'_',assay.df.dcast$ID,'_',assay.df.dcast$Lot)
  cv.df <- (assay.obj@stats$cv.df)
  cv.df$index <- paste0(assay.df$ID,'_',assay.df$Lot)
  cv.df.melted <- reshape2::melt(cv.df)
  colnames(cv.df.melted) <- c('index','DayOperator','CV')
  cv.df.melted$xid <- paste0(cv.df.melted$DayOperator,'_',cv.df.melted$index)
  cv.df.melted$Decision <- ifelse(cv.df.melted$CV < 20,'Keep','Discard')

  discard.xid <- cv.df.melted$xid[cv.df.melted$Decision == 'Discard']
  adfmelt.xid <- paste0(assay.df.melt$DayOperator,'_',assay.df.melt$ID,'_',assay.df.melt$Lot)

  assay.df.melt$xid <- adfmelt.xid
  assay.df.melt[assay.df.melt$xid %in% discard.xid,'value'] <- NA
  assay.df.data <- reshape2::dcast(assay.df.melt, formula = ID + Lot~variable,fun.aggregate = sum,value.var = "value")


  assay.obj@data <- assay.df.data
  assay.obj@melted.data <- assay.df.melt

  assay.obj@cv.table <- cv.df.melted

  assay.obj@outlier.rm <- "yes"
  assay.obj@outlier.rm.method <- c(assay.obj@outlier.rm.method, "cv")

  return(assay.obj)

})


#' @title Melt Assay Dataset
#'
#' @description This function is a wrapper for the reshape2::melt() function and splits column of Day, Operator, and Replicate information into 3 separate columns.
#'
#' @param assay.df Imported data.frame consisting of the following columns: 'ID','Lot', and columns identifying the Day, Operator and Replicate like so: 'D1_Op2_3' to indicate Day 1, operator 2, replicate 3.
#' @param exp.name Experiment name (as a string). To be used to distinguish between experiments for when melted assays are combined using rbind().
#'
#'
#' @author Emma Gail
#'
#' @return A melted data.frame with the following columns "ID","Lot","variable","value","Day",'Operator","Replicate","Category","DayOperator"
#'
#' @examples
#' assay.df.melted <- assayMelt(assay.df = lognormAssay, exp.name = 'Experiment1')
#' head(assay.df.melted)
#'
#' @export
#'

assayMelt <- function(assay.df, exp.name) {


  if (is.numeric(assay.df$ID) == TRUE) {
    assay.df$ID <- as.character(assay.df$ID)
  }


  if ("outlier.rm" %in% names(assay.df)) {
    assay.df_melted <- suppressMessages(reshape2::melt(assay.df$outlier.rm))
  } else {
    assay.df_melted <- suppressMessages(reshape2::melt(assay.df))
  }

  split_cols <- data.frame(stringr::str_split_fixed(assay.df_melted$variable,
                                                    "_", 3))
  colnames(split_cols) <- c("Day", "Operator", "Replicate")


  assay.df_melted <- cbind(assay.df_melted, split_cols)

  assay.df_melted$Category <- rep(exp.name, nrow(assay.df_melted))

  assay.df_melted$DayOperator <- paste0(assay.df_melted$Day, assay.df_melted$Operator)

  return(assay.df_melted)

}  #end function assayMelt

#' @title Evaluate the Assays with Boxplots
#'
#' @description This function produces a boxplot based on the variable chosen in order to visualize any analytical variability.
#'
#' @param assay.obj An ImmunoAssay object imported by importAssay
#' @param var Variable to examine in the plot. Either "Day" or "Operator".
#'
#' @author Emma Gail
#'
#' @return A boxplot as generated by \code{ggplot2}
#'
#'
#' @name evalBoxplot
#' @aliases evalBoxplot,ImmunoAssay-method
#' @docType methods
#'
#' @examples
#' assay.obj <- importAssay(lognormAssay, exp.name = 'Experiment1')
#' evalBoxplot(assay.obj,var='Day') #visualize days on boxplot
#'
#' @export
#'

setGeneric('evalBoxplot',function(assay.obj, var = c('Day','Operator')){
  standardGeneric('evalBoxplot')
})

#'
setMethod('evalBoxplot','ImmunoAssay',function(assay.obj, var = c('Day','Operator')){

    assay.df.melted <- assay.obj@melted.data


  if(length(var) > 1){
    warning("Only 1 variable allowed at a time")
    return(NULL)
  }


  .e <- environment()

  p <- ggplot2::ggplot(assay.df.melted, ggplot2::aes(x=.data[[var]], y=value, fill=Category), environment = .e) +
    ggplot2::stat_boxplot(geom ='errorbar',width=0.5) + ggplot2::geom_boxplot() +
    #ggplot2::facet_wrap(~assay.df.melted$Category) +
    ggplot2::labs(title = paste0("Between ",var," Differences"), x=var, y='Relative Unit', fill="Category")

  return(p)

})


#' @title Normality Evaluation
#'
#' @description This function evaluates the normality of the melted immunoassay dataset. In order to determine whether or not the distribution is normal, two tests are performed: the Shapiro Wilk test for normality and the test for skewness. See \link[stats]{shapiro.test} and \link[e1071]{skewness} for details. In order to for a "nonparametric" recommendation to be made, the SW test must be significant (below desire value) and the absolute value skewness must be above the desired value. If only one or neither of these conditions are met, then the recommendation will be "normal".
#'
#' @param assay.obj An ImmunoAssay object imported by importAssay
#' @param category If assay.df.melted consists of more than 1 dataset, choose the category here to split dataset
#' @param data.transf Should the data should be transformed before normality is evaluated
#' @param transf.method If data.transf is TRUE, which method should be used. Can choose between 'log10' and 'ln'.
#' @param excl.outliers Should outliers be excluded from this analysis? If TRUE, data points which lie beyond the extremes of the whiskers in boxplot will be excluded, see boxplot.stats for details.
#' @param hist Should a histogram be outputted? TRUE/FALSE
#' @param p.val Value to be used for cutoff for Shapiro-Wilks test. Defaults to 0.05.
#' @param skew Value to be used to determine skewness. Defaults to 1.
#' @param return.object If FALSE, only the plot is returned and the stats are returned as a list.
#'
#'
#' @aliases evalNorm,ImmunoAssay-method
#' @docType methods
#'
#' @return If \code{return.object==FALSE}, only the plot is returned and the stats are returned as a list. Otherwise, an object of the class ImmunoAssay is returned.
#'
#' @author Emma Gail
#'
#' @examples
#' assay.obj <- importAssay(lognormAssay, exp.name = 'Experiment1')
#' assay.obj <- evalNorm(assay.obj, category = 'Experiment1',
#' data.transf = TRUE, transf.method = 'log10')
#'
#'
#' @export


setGeneric('evalNorm',function(assay.obj, category = NULL, data.transf = FALSE, transf.method = c('log10','ln'), excl.outliers = FALSE, hist = TRUE, p.val = 0.05,skew = 1, return.object=TRUE){
  standardGeneric('evalNorm')
})


setMethod('evalNorm','ImmunoAssay',function(assay.obj, category = NULL, data.transf = FALSE, transf.method = c('log10','ln'), excl.outliers = FALSE, hist = TRUE, p.val = 0.05,skew = 1, return.object=TRUE){

  assay.df.melted <- assay.obj@melted.data

  if(!is.null(category)){
    assay.df.melted <- assay.df.melted[assay.df.melted$Category == category,]
  } #end if(!is.null(category)){

  if(data.transf == TRUE){

    if(transf.method == 'log10'){
      assay.df.melted$value <- log10(assay.df.melted$value)
    } else if(transf.method == 'ln'){
      assay.df.melted$value <- log(assay.df.melted$value)
    } else {
      warning('Unknown data transformation method --> using original dataset\n')
    }

  } else { #end if(data.transf == TRUE){

    transf.method <- 'None'

  }


  if(excl.outliers == TRUE){
    assay.df.melted <- excludeOutliers(assay.df.melted)
  }


  if(hist == TRUE){

    p <- ggplot2::ggplot(assay.df.melted, ggplot2::aes(x=assay.df.melted$value)) + ggplot2::geom_histogram(ggplot2::aes(x = assay.df.melted$value, y = ..density..), colour="black", fill="white") +
      ggplot2::stat_function(fun = stats::dnorm, args = list(mean = mean(assay.df.melted$value, na.rm = TRUE), sd = stats::sd(assay.df.melted$value, na.rm = TRUE))) +
      ggplot2::labs(title = 'Normality Investigation')

  } else {#end if(hist == TRUE){

    p <- 'No plot rendered'

  }

  sw.results <- stats::shapiro.test(assay.df.melted$value)
  skew.val <- e1071::skewness(assay.df.melted$value, na.rm = TRUE)

  if(is.na(sw.results$p.value) == TRUE){
    #if is not a number --> will need to return a different recommendation
    recommendation <- 'not determined'
    warning('No p-value determined from Shapiro-Wilks test, no recommendation given')

  } else if(sw.results$p.value < p.val && abs(skew.val) > skew){
    recommendation <- 'nonparametric'
  } else {#end if(sw.results$p.value < p.val && abs(skew.val) > skew){
    recommendation <- 'normal'
  }

  assay.obj@stats$transformation = transf.method
  assay.obj@stats$sw.results=sw.results
  assay.obj@stats$skewness=skew.val
  assay.obj@stats$recommendation=recommendation

  if(hist==TRUE){
    print(p)
  }

  if(return.object==TRUE){
    return(assay.obj)
  } else {
    return(list(transf.method=transf.method,
    sw.results=sw.results,
    skew.val=skew.val,
    recommendation=recommendation))
  }


})


#' Exclude Outliers from Melted Assay Dataframe
#'
#' This function excludes outliers from the assay dataframe based on grDevices::boxplot.stats(). This outlier removal method is based on Tukey's test where outliers are removed if outside the established interquartile range.
#'
#' @param assay.df.melted A data.frame produced by assayMelt()
#' @param data.transf Should the data should be transformed before normality is evaluated
#' @param transf.method If data.transf is TRUE, which method should be used. Can choose between 'log10' and 'ln'.
#'
#' @author Emma Gail
#'
#' @aliases excludeOutliers
#'
#' @return A melted data.frame in the same format produced by assayMelt() with outliers excluded
#'
#'
#' @examples
#' assay.df.melted <- assayMelt(assay.df = lognormAssay, exp.name = 'Experiment1')
#' excludeOutliers(assay.df.melted, data.transf = TRUE, transf.method = 'log10')
#'
#' @export

excludeOutliers <- function(assay.df.melted, data.transf = FALSE, transf.method = c('log10','ln')){

  assay.df.melted.mega <- assay.df.melted

  assay.df.melted.out.rm.mega <- data.frame()

  for(category in unique(assay.df.melted.mega$Category)){

    assay.df.melted <- assay.df.melted.mega[assay.df.melted.mega$Category == category,]


    if(data.transf == TRUE){
      if(transf.method == 'log10'){
        assay.df.melted$value <- log10(assay.df.melted$value)
      }
    }


    outliers <- grDevices::boxplot.stats(assay.df.melted$value)$out

    assay.df.melted.out.rm <- assay.df.melted[!assay.df.melted$value %in% grDevices::boxplot.stats(assay.df.melted$value)$out,]

    assay.df.melted.out.rm.mega <- rbind(assay.df.melted.out.rm.mega,assay.df.melted.out.rm)

  } #end for(category in unique(assay.df_melted$Category)){


  return(assay.df.melted.out.rm.mega)

} #end excludeOutliers function



#' @title Mixed model wrapper for assay dataframe
#'
#' @description This function is a wrapper for the lmer() function to produce a table of results. Here, the sample ID is included as a random intercept effect, then the results of the fixed effect model estimates are reported together with 95\% confidence intervals and t statistics.
#'
#' @param assay.df.melted A data.frame produced by assayMelt()
#' @param var Variable to look at. Either "Day" or "Operator".
#'
#' @author Emma Gail
#' @author Lidija Turkovic
#'
#' @aliases mixedModel
#'
#' @return A data.frame with the following columns: "Parameter", "Estimate", "LowerCI", "UpperCI", "Tstat"
#'
#' @examples
#' assay.df.melted <- assayMelt(assay.df = lognormAssay, exp.name = 'Experiment1')
#' mixedModel(assay.df.melted, var = 'Day')
#'
#' @export

mixedModel <- function(assay.df.melted, var = c('Operator','Day')){

  if(var == 'Operator'){
    mod <- suppressWarnings(suppressMessages(lme4::lmer(value ~ Operator + (1|ID), data = assay.df.melted)))
  } else if(var == 'Day'){
    mod <- suppressWarnings(suppressMessages(lme4::lmer(value ~ Day + (1|ID), data = assay.df.melted)))

  } else {
    warning('No valid variable chosen, returning NA')
    return(NA)

  }

  ci <- suppressMessages(stats::confint(mod))
  ci <- ci[3:nrow(ci),]
  colnames(ci) <- c('LowerCI','UpperCI')
  tab <- as.data.frame(cbind(summary(mod)$coefficients[,1], ci, summary(mod)$coefficients[,3]))
  tab <- cbind(Parameter = rownames(tab), tab)
  rownames(tab) <-  NULL
  colnames(tab) <- c("Parameter", "Estimate", "LowerCI", "UpperCI", "Tstat")

  return(tab)

} #end mixedModel function


#' @title Calculate screening cut point values for scp()
#'
#' @description This function calculates the values needed for the output of the scp() data.frame
#'
#' @param assay.values List of selected values from the assay data.frame such as selected values from assayMelt()
#' @param conf.level Decimal describing level of confidence to be used for confidence interval calculation. Defaults to 0.95
#' @param distrib Distribution selection to determine the cut point calculation. Two options: 'nonparametric' or 'normal'
#' @param transf.method Transformation method used. The inverse will be calculated as part of the output.
#'
#' @author Emma Gail
#'
#' @return A data.frame cotaining the values: "mean", "sd", "distrib", "cp", "mean.conf.int1", "mean.conf.int2", "cp.conf.int1", "cp.conf.int2"
#'
#' @examples
#'
#' assay.df.melted <- assayMelt(assay.df = lognormAssay, exp.name = 'Experiment1')
#' assay.values <- assay.df.melted[assay.df.melted$DayOperator == 'D1Op1',]$value
#' #This function assumes that the data has already been transformed.
#' scp.df <- calcScpValues(assay.values = log10(assay.values), distrib = 'normal',
#' transf.method = 'log10')
#'
#'
#'
#' @export

calcScpValues <- function(assay.values, conf.level = 0.95, distrib = c('nonparametric','normal'), transf.method = c('log10','ln')){

  assay.mean <- mean(assay.values, na.rm = TRUE)
  assay.sd <- stats::sd(assay.values, na.rm = TRUE)


  if(distrib == 'nonparametric'){
    cp <- stats::quantile(assay.values,0.95,na.rm = TRUE)
  } else if(distrib == 'normal'){
    cp <- assay.mean + (1.645*assay.sd)
  } else {
    warning('Unknown distribution selection, please choose "nonparametric" or "normal"\n')
  }


  if(!is.na(transf.method)){

    if(transf.method == 'log10'){
      cp.inv <- 10^(cp)
    } else if(transf.method == 'ln'){
      cp.inv <- exp(cp)
    } else {
      cp.inv <- cp
      #warning('Unknown transformation method')
    }
  } else {
    cp.inv <- cp
  }

  conf.int1 <- assay.mean - (1.96*(assay.sd/sqrt(length(assay.values)))) #lower confidence interval
  conf.int2 <- assay.mean + (1.96*(assay.sd/sqrt(length(assay.values)))) #upper confidence interval


  cp.conf.int1 <- cp - (1.96*(assay.sd/sqrt(length(assay.values))))
  cp.conf.int2 <- cp + (1.96*(assay.sd/sqrt(length(assay.values))))

  if(transf.method == "log10"){
    cp.conf.int1 <- 10^(cp.conf.int1)
    cp.conf.int2 <- 10^(cp.conf.int2)
  }

  runs.stats.df <- data.frame(list(mean=assay.mean,
                                   sd=assay.sd,
                                   distrib=distrib,
                                   cp=cp,
                                   cp.inv=cp.inv,
                                   mean.conf.int1=conf.int1,
                                   mean.conf.int2=conf.int2,
                                   cp.conf.int1=cp.conf.int1,
                                   cp.conf.int2=cp.conf.int2))


  return(runs.stats.df)

} #end function calcScpValues


#' @title Calculate screening cut point
#'
#' @description This function will calculate the screening cut point from the melted assay.df
#'
#' @param assay.obj An ImmunoAssay object imported by importAssay
#' @param category If assay.obj consists of more than 1 dataset, choose the category here to split dataset
#' @param data.transf Should the data should be transformed before the cut point is calculated
#' @param distrib Distribution selection to determine the cut point calculation. Two options: 'nonparametric' or 'normal'
#' @param transf.method If data.transf is TRUE, which method should be used. Can choose between 'log10' and 'ln'.
#' @param rm.out Should outliers be excluded from this analysis?
#'
#' @author Emma Gail
#'
#' @return An object of the class ImmunoAssay
#'
#' @aliases scp,ImmunoAssay-method
#' @docType methods
#'
#' @examples
#' assay.obj <- importAssay(assay.df = lognormAssay, exp.name = 'Experiment1')
#' assay.obj <- scp(assay.obj, category = 'Experiment1', distrib = 'normal',
#' data.transf = TRUE, transf.method = 'log10', rm.out = FALSE)
#'
#' @export

setGeneric("scp", def = function(assay.obj, category = NULL, distrib = c("normal",
                                                                         "nonparametric"), data.transf = FALSE, transf.method = c("log10", "ln"),
                                 rm.out = FALSE) {
  standardGeneric("scp")
})

#'
setMethod("scp", "ImmunoAssay", function(assay.obj, category = NULL, distrib = c("normal",
                                                                                 "nonparametric"), data.transf = FALSE, transf.method = c("log10", "ln"),
                                         rm.out = FALSE) {

  assay.df.melted <- assay.obj@melted.data

  if (length(transf.method) > 1) {
    transf.method <- "else"
  }

  if (is.na(transf.method)) {
    transf.method <- "else"
  }

  assay.df.melted.copy <- assay.df.melted

  if (!is.null(category)) {
    assay.df.melted <- assay.df.melted[assay.df.melted$Category ==
                                         category, ]
  }

  day_op_combined <- paste0(assay.df.melted$Day, assay.df.melted$Operator)

  if (!is.na(transf.method))
  {
    if ((length(transf.method) == 1))
    {

      if (transf.method == "log10") {
        assay.df.melted$value <- log10(assay.df.melted$value)
      } else if (transf.method == "ln") {
        assay.df.melted$value <- log(assay.df.melted$value)
      } else {
        if (data.transf == TRUE) {
          warning("Unknown transformation method")  #only if data.tranf is TRUE will this error appear (not relevant otherwise)
        }
      }
    }  #end if((length(transf.method) == 1)){
  }  #end if(!is.na(transf.method)){


  if (rm.out == TRUE) {
    assay.df.melted <- excludeOutliers(assay.df.melted = assay.df.melted.copy,
                                       data.transf = data.transf, transf.method = transf.method)
  }

  runs.stats.df <- data.frame()


  for (uname in unique(day_op_combined)) {

    assay.values <- assay.df.melted[day_op_combined == uname, ]$value

    rs.df <- calcScpValues(assay.values = assay.values, distrib = distrib,
                           transf.method = transf.method)

    op.day.df <- data.frame(list(OperatorDay = uname))

    rs.df <- cbind(op.day.df, rs.df)

    runs.stats.df <- rbind(runs.stats.df, rs.df)

  }  #end for(uname in unique(day_op_combined)){


  runs.stats.df <- stats::na.omit(runs.stats.df)

  runs.stats.df2 <- data.frame()

  for (operator in levels(assay.df.melted$Operator)) {

    assay.values <- assay.df.melted[assay.df.melted$Operator == operator,
                                    ]$value

    if (length(assay.values) > 0)
    {

      rs.df <- calcScpValues(assay.values = assay.values, distrib = distrib,
                             transf.method = transf.method)

      op.day.df <- data.frame(list(OperatorDay = paste0(operator,
                                                        "_all")))
      rs.df <- cbind(op.day.df, rs.df)

      runs.stats.df2 <- rbind(runs.stats.df2, rs.df)

    }  #end if(length(assay.values) > 0){



  }  #end for(operator in levels(runs.stats.df$Operator)){

  assay.values <- assay.df.melted$value

  op.day.df <- data.frame(list(OperatorDay = "all"))
  rs.df <- calcScpValues(assay.values = assay.values, distrib = distrib,
                         transf.method = transf.method)

  rs.df <- cbind(op.day.df, rs.df)

  runs.stats.df2 <- rbind(runs.stats.df2, rs.df)

  runs.stats.df <- rbind(runs.stats.df, runs.stats.df2)

  if (data.transf == FALSE) {
    # if no data transformation --> remove the inverse of the cut point
    # here (same value as the cut point)
    runs.stats.df$cp.inv <- NULL
  } else {
    # if there is a data transformation --> tell the user that the mean and
    # sd have been transformed
    names(runs.stats.df)[names(runs.stats.df) == "mean"] <- "mean.transf"
    names(runs.stats.df)[names(runs.stats.df) == "sd"] <- "sd.transf"
  }  #end else to if(data.transf == FALSE){


  rownames(runs.stats.df) <- 1:nrow(runs.stats.df)

  assay.obj@scp.table <- runs.stats.df

  return(assay.obj)

})


#' @title Generate forest plot of SCP values
#'
#' @description This function creates a forest plot of the calculated screening cut points. The scp function must be called prior to this function.
#'
#' @param assay.obj An ImmunoAssay object imported by importAssay
#' @param ... Additional arguments for forestplot() function
#'
#' @author Emma Gail
#'
#' @aliases scpForestPlot,ImmunoAssay-method
#' @docType methods
#'
#' @return A forestplot
#'
#' @examples
#' assay.obj <- importAssay(assay.df = lognormAssay, exp.name = 'Experiment1')
#' assay.obj <- scp(assay.obj, category = 'Experiment1', distrib = 'normal',
#' data.transf = TRUE, transf.method = 'log10', rm.out = FALSE)
#' scpForestPlot(assay.obj)
#'
#' @export

setGeneric('scpForestPlot',def = function(assay.obj, ...){
  standardGeneric('scpForestPlot')
})

#'
setMethod('scpForestPlot','ImmunoAssay',function(assay.obj, ...){
  cp.table <- assay.obj@scp.table

  if(nrow(cp.table) == 0){
    warning('Please run scp() function first.\n')
    return(assay.obj)
  }

  # Prepare the text for the forest plot
  fp.text.df <- as.character(cp.table$OperatorDay)
  # Prepare the data for the forest plot
  if('cp.inv' %in% colnames(cp.table)){
    mean.list <- cp.table$cp.inv
  } else {
    mean.list <- cp.table$cp
  }

  fp.data.df <- data.frame(mean=mean.list,upper=cp.table$cp.conf.int2,lower=cp.table$cp.conf.int1)
  # Create a list of FALSE statements for the is.summary option
  rep_false <- rep(FALSE,length(fp.text.df)-1)

  # Render the forest plot
  fp <- forestplot::forestplot(labeltext = fp.text.df, # Text for y-axis
                   mean = fp.data.df,
                   zero = NA,
                   is.summary=c(rep_false,TRUE),
                   ...)

  # add in diamond shape

  return(fp)

})


#'Simulated Lognormal Dataset
#'
#'This is a simulated dataset that using a lognormal distribution
#'
#'@docType data
#'
#'@usage data(lognormAssay)
#'
"lognormAssay"

