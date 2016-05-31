# # Merge all the results.csv files
###################################
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
options(warn=1)

dataset <- "agreement"
dataset <- "disagreement"
dataset <- 'iris'


path <- file.path("out", dataset)

# one view sees 4, the other sees 5 clusters
files.name <- list.files(path=path, pattern=".csv", recursive=TRUE, full.names=TRUE)
files.df <- lapply(files.name, read.delim, header=FALSE, sep='\t')
df <- do.call(rbind, files.df)

names(df) <- c("dir", "model", "T", "negloglike", "ARI")
df <- df[c('model', 'T','negloglike')]
levels(df$model)[levels(df$model)=="DP"] <- "dual-DP"
levels(df$model)[levels(df$model)=="fixed"] <- "dual-fixed"
levels(df$model)[levels(df$model)=="norole"] <- "single"

df <- df[df$T %in% c(10,20,30,40,50,60,70,80,90,100),]

# Summarize means and variances for every group (U,T)
df <- ddply(df,.(model, T), 
                summarise, 
                neglog.mean = mean(negloglike), neglog.sd = sd(negloglike))

names(df)[2] <- "threads"

# Plot accuracy 
p1 <- ggplot(df, aes(x= threads, y = neglog.mean, color=model)) + 
  geom_line(aes(linetype=model)) + 
  geom_point(size=2) +
  scale_color_manual(name= element_blank(),
                     values=c('black', 'red', 'blue'),
                     labels=c("dual-DP", "dual-fixed", "single")) +
  scale_linetype_manual(name= element_blank(),
                        values=c('solid', 'dashed', 'dotted'),
                        labels=c("dual-DP", "dual-fixed", "single")) +
  guides(color=guide_legend(override.aes=list(shape=c(NA,NA, NA), linetype=c('solid', 'dashed', 'dotted')))) +
  
  xlab('threads') +
  ylab('negative loglikelihood') +
  labs(color=NULL) +
  geom_errorbar(aes(ymin=neglog.mean-neglog.sd/2, ymax=neglog.mean+neglog.sd/2), width=1.2)+
  scale_x_continuous(breaks=df$threads)+
  theme(text = element_text(size = 13),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        aspect.ratio = 6/9,
        legend.position=c(0.8, 0.8))

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

names(df) <- c("dir", "model", "T", "negloglike", "ARI")
df = df[c('model','T','ARI')]
levels(df$model)[levels(df$model)=="DP"] <- "dual-DP"
levels(df$model)[levels(df$model)=="fixed"] <- "dual-fixed"
levels(df$model)[levels(df$model)=="norole"] <- "single"

df <- df[df$mode != "single",]
df <- df[df$T %in% c(10,20,30,40,50,60,70,80,90,100),]

# Summarize means and variances for every group (U,T)
df <- ddply(df,.(model, T), 
            summarise, 
            ARI.mean = mean(ARI), ARI.sd = sd(ARI))

names(df)[2] <- "threads"

# Plot accuracy 
p2 <- ggplot(df, aes(x= threads, y = ARI.mean, color=model)) + 
     geom_line(aes(linetype=model)) + geom_point(size=2) +
    scale_color_manual(name= element_blank(),
                       values=c('black', 'red'),
                       labels=c("dual-DP", "dual-fixed")) +
    scale_linetype_manual(name= element_blank(),
                          values=c('solid', 'dashed'),
                          labels=c("dual-DP", "dual-fixed")) +
    guides(color=guide_legend(override.aes=list(shape=c(NA,NA), linetype=c('solid', 'dashed')))) +
     geom_errorbar(aes(ymin=ARI.mean-ARI.sd/2, ymax=ARI.mean+ARI.sd/2), width=1.2) +
     xlab('threads') +
     ylab('Adjusted Rand Index') +
     labs(color=NULL) +
     scale_x_continuous(breaks=df$threads) +
     scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0,1.1)) +
     theme(text = element_text(size = 13),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank(),
        aspect.ratio = 6/9,
        legend.position=c(0.8, 0.2))



g <- plot_grid(p1, p2, ncol = 2, 
               align = 'v') 

ggsave("./doc/ComputStat submission/Fig8_results_iris_bw.eps", 
       device="eps", 
       plot=g,
       #width = 20, height = 8, units = "cm")
       width = 174, height = 70, units = "mm")

# Springer instructions
# For most journals the figures should be 39 mm, 84 mm, 129 mm, or 174 mm wide and not higher than 234 mm.