## ----setup, include=FALSE-----------------------------------------------------
#knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(cache=TRUE)

## -----------------------------------------------------------------------------

assay.obj <- scp(assay.obj = assay.obj,
                category = 'Experiment1',
                distrib = 'normal',
                data.transf = TRUE,
                transf.method = 'log10',
                rm.out = FALSE)

print(assay.obj@scp.table)


## -----------------------------------------------------------------------------

assay.df.melted <- assay.obj@melted.data
mod <- aov(value ~ DayOperator, data=assay.df.melted)
aov.p <- summary(mod)[[1]][[1,"Pr(>F)"]]

varhom <- car::leveneTest(value ~ DayOperator, data=assay.df.melted)
levene.p <- varhom$`Pr(>F)`[1]

data.frame(aov=aov.p, levene=levene.p)


## -----------------------------------------------------------------------------

scp.compare <- data.frame()

for(distrib in c('normal','nonparametric')){
  for(exclOut in c(T,F)){
  for(dataTransf in c(T,F)){

      if(dataTransf == TRUE){
        for(transfMeth in c('log10')){
  
          
                 assay.obj2 <- scp(assay.obj = assay.obj,
                category = 'Experiment1',
                distrib = distrib,
                transf.method = transfMeth,
                data.transf = dataTransf,
                rm.out = exclOut)
    
    scp.row <- assay.obj2@scp.table[nrow(assay.obj2@scp.table),]
    
    if('cp.inv' %in% colnames(scp.row)){
      cp.row <- scp.row[,'cp.inv']
    } else {
      cp.row <- scp.row[,'cp']
    }
    
    scp.compare <- rbind(scp.compare, data.frame(cp=round(cp.row,digits=3), distrib=distrib,excludeOutliers=as.character(exclOut),dataTransformed=as.character(dataTransf), transfMethod=transfMeth))
          
          
  
        }
    
      } else {
        
            assay.obj2 <- scp(assay.obj = assay.obj,
                category = 'Experiment1',
                distrib = distrib,
                data.transf = dataTransf,
                rm.out = exclOut)
    
    scp.row <- assay.obj2@scp.table[nrow(assay.obj2@scp.table),]
    
    if('cp.inv' %in% colnames(scp.row)){
      cp.row <- scp.row[,'cp.inv']
    } else {
      cp.row <- scp.row[,'cp']
    }
    
    scp.compare <- rbind(scp.compare, data.frame(cp=round(cp.row,digits=3), distrib=distrib,excludeOutliers=as.character(exclOut),dataTransformed=as.character(dataTransf), transfMethod='NA'))
        
      }
    
    
  }
  } 
}

scp.compare





## -----------------------------------------------------------------------------

sessionInfo()


