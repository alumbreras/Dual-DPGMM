# # Merge all the results.csv files
###################################
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
options(warn=1)

dataset <- 'iris'
dataset <- 'clear'
dataset <- 'overlapped'
dataset <- 'confused_features'
dataset <- "agreement"

path <- file.path("out", dataset)

# one view sees 4, the other sees 5 clusters
files.name <- list.files(path=path, pattern=".csv", recursive=TRUE, full.names=TRUE)
files.df <- lapply(files.name, read.delim, header=FALSE, sep='\t')
df <- do.call(rbind, files.df)

names(df) <- c("model", "T", "negloglike", "ARI")
df <- df[c('model', 'T','negloglike')]
levels(df$model)[levels(df$model)=="DP"] <- "dual-DP"
levels(df$model)[levels(df$model)=="fixed"] <- "dual-fixed"
levels(df$model)[levels(df$model)=="norole"] <- "single"

df <- df[df$T %in% c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150, 200,300,400,500),]

# Summarize means and variances for every group (U,T)
df <- ddply(df,.(model, T), 
                summarise, 
                neglog.mean = mean(negloglike), neglog.sd = sd(negloglike))


# Plot accuracy 
p1 <- ggplot(df, aes(x= df$T, y = neglog.mean, color=model)) + geom_line() + geom_point(size=2) +
  xlab('Threads') +
  ylab('negative loglikelihood') +
  labs(color=NULL) +
  geom_errorbar(aes(ymin=neglog.mean-neglog.sd/2, ymax=neglog.mean+neglog.sd/2), width=1.2)+
  scale_x_continuous(breaks=df$T)+
  theme(text = element_text(size = 15),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        aspect.ratio = 6/9,
        legend.position=c(0.9, 0.8))

print(p1)
g1 <- grid.grab() 


# see also:
#theme_classic()
#panel.grid.major = element_blank()
#panel.grid.minor = element_blank()


################################
# ARI
################################
# Get all the experiments with 50 users
files.name <- list.files(path=path, pattern=".csv", recursive=T, full.names=T)
files.df <- lapply(files.name, read.delim, header=FALSE, sep='\t')
df <- do.call(rbind, files.df)

names(df) <- c("model", "T", "negloglike", "ARI")
df = df[c('model','T','ARI')]
levels(df$model)[levels(df$model)=="DP"] <- "dual-DP"
levels(df$model)[levels(df$model)=="fixed"] <- "dual-fixed"
levels(df$model)[levels(df$model)=="norole"] <- "single"

df <- df[df$mode != "single",]
df <- df[df$T %in% c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150, 200,300,400,500),]

# Summarize means and variances for every group (U,T)
df <- ddply(df,.(model, T), 
            summarise, 
            ARI.mean = mean(ARI), ARI.sd = sd(ARI))


# Plot accuracy 
p2 <- ggplot(df, aes(x= df$T, y = ARI.mean, color=model)) + 
     geom_line() + geom_point(size=2) +
     geom_errorbar(aes(ymin=ARI.mean-ARI.sd/2, ymax=ARI.mean+ARI.sd/2), width=1.2) +
     xlab('Threads') +
     ylab('Adjusted Rand Index') +
     labs(color=NULL) +
     scale_x_continuous(breaks=df$T) +
     scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0,1.1)) +
     theme(text = element_text(size = 15),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        aspect.ratio = 6/9,
        legend.position=c(0.9, 0.2))

print(p2)
g2 <- grid.grab() 

grid.arrange(g1, g2, widths=c(0.5,0.5))  

