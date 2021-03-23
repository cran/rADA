## ----setup, include=FALSE-----------------------------------------------------
#knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(cache=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  
#  # Current version
#  install.packages('rADA')
#  
#  # Development version
#  devtools::install_github('egmg726/rADA')
#  

## ---- message=FALSE-----------------------------------------------------------

library(rADA)
library(forestplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)


## -----------------------------------------------------------------------------

data(lognormAssay)


## -----------------------------------------------------------------------------
head(lognormAssay)

## -----------------------------------------------------------------------------

assay.obj <- importAssay(lognormAssay, exp.name = 'Experiment1')


## -----------------------------------------------------------------------------

assay.obj <- calcCvStats(assay.obj)


## -----------------------------------------------------------------------------
names(assay.obj@stats)

## -----------------------------------------------------------------------------

table(assay.obj@stats$is.adj.df$D1Op2)


## -----------------------------------------------------------------------------

table(assay.obj@stats$is.adj.df$D2Op1)


## -----------------------------------------------------------------------------

head(assay.obj@melted.data, n = 7)


## -----------------------------------------------------------------------------

evalBoxplot(assay.obj, var = 'Day')


## -----------------------------------------------------------------------------

evalBoxplot(assay.obj, var = 'Operator') + theme_minimal()


## -----------------------------------------------------------------------------

evalBoxplot(assay.obj, var = 'Replicate') + ggplot2::theme_minimal() + ggplot2::scale_fill_manual(values='#00a4b2')


## -----------------------------------------------------------------------------

# Create the boxplot for the top of the plot
p1 <- ggplot(assay.obj@melted.data, aes(x=Category,y=value)) + stat_boxplot(geom ='errorbar',width=0.2, size = 1.5) +geom_boxplot(size = 1.1) + coord_flip() +  scale_y_continuous(limits = c(0, 40)) +
    theme(
        axis.title=element_text(size=12,face="bold"),
       panel.background = element_blank(), 
       panel.grid = element_blank(),
       axis.title.y = element_blank(),
       axis.title.x = element_blank(),
       axis.ticks.x=element_blank(),
       axis.ticks.y=element_blank(),
       panel.border = element_rect(colour = "black", fill=NA, size=2),
       axis.text.x=element_blank(),
       axis.text.y=element_blank())


# Create the histogram for the bottom of the plot
p2 <- ggplot(assay.obj@melted.data, aes(x=value)) + geom_histogram(aes(x = value, y = ..density..), colour="black", fill="#6c78a7", size = 2) +  scale_x_continuous(limits = c(0, 40)) +
  stat_function(fun = dnorm, args = list(mean = mean(assay.obj@melted.data$value, na.rm = TRUE), sd = sd(assay.obj@melted.data$value, na.rm = TRUE)), size = 2) +
  theme(
        axis.title=element_text(size=12,face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
       panel.background = element_blank(), 
       panel.grid = element_blank())


## -----------------------------------------------------------------------------

grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))


## -----------------------------------------------------------------------------


evalNorm(assay.obj = assay.obj, category = 'Experiment1', data.transf = FALSE, return.object=FALSE)


## -----------------------------------------------------------------------------


assay.obj <- evalNorm(assay.obj = assay.obj, category = 'Experiment1', data.transf = TRUE, transf.method = 'log10')



## -----------------------------------------------------------------------------

names(assay.obj@stats)


## -----------------------------------------------------------------------------

# Results from the Shapiro-Wilks test
print(assay.obj@stats$sw.results)

# Calculated skewness value
print(assay.obj@stats$skewness)

# Recommendation based on the previous 2 values
print(assay.obj@stats$recommendation)


## -----------------------------------------------------------------------------

assay.obj <- evalNorm(assay.obj = assay.obj, category = 'Experiment1', data.transf = FALSE, excl.outliers = TRUE)


## -----------------------------------------------------------------------------

assay.obj <- scp(assay.obj = assay.obj,
    category = 'Experiment1',
    distrib = 'nonparametric',
    data.transf = FALSE,
    rm.out = FALSE)

print(assay.obj@scp.table)


## -----------------------------------------------------------------------------

assay.obj <- scp(assay.obj = assay.obj,
                category = 'Experiment1',
                distrib = 'normal', #assay.norm.eval$recommendation,
                data.transf = TRUE,
                transf.method = 'log10',
                rm.out = FALSE)

print(assay.obj@scp.table)


## -----------------------------------------------------------------------------

scpForestPlot(assay.obj)


## -----------------------------------------------------------------------------

sessionInfo()


