library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(LaplacesDemon)

setwd('/home/stephen/Git_Repositories/VDJ_Recombination/Visualisation')
my_dt <- function(param,n){

  pp <- param[1]
  pn <- param[2]
  pe <- param[3]
  
  if (abs(pn - pp) < 0.00001){
    pn = pe-0.00001
  }
  dens <<- (n==0)*10^-9 + (n>0)*((pe^2*pn*pp^2*((n-1)*(pe*(1-pp))^(n-2)))/((pn-pp)*(pe+pp-pe*pp)^(n+2))) + 
      ((pe*pn*pp^2*(pe*(1-pn))^n)/((pn-pp)^2*(pe-pe*pn)*(pe+pn-pe*pn)^(n+1))) - 
      ((pe*pn*pp^2 * (pe*(1-pp))^n * (pe-2*pn + 3*pp + 2* pe*pn - 3 *pe*pp))/((pn - pp)^2*(pe-pe*pp)*(pe + pp - pe*pp)^(n+2)))  
  
  dens[is.na(dens)] <- 10^-9
  dens
}

objfun <- function(param,n){
  tbl <- table(n)
#  print(-sum( log(my_dt(param, as.numeric(names(tbl))))*unname(tbl) ))
  -sum( log(my_dt(param, as.numeric(names(tbl))))*unname(tbl) )
}

objfun2 <- function(param,n){
  
  #  print(-sum( log(my_dt(param, as.numeric(names(tbl))))*unname(tbl) ))
  -sum( log(my_dt(param, n$Length))*n$Freq )
}


my_optimsearch <- function(in_sample){
  rez<-optim(c(0.25,0.1,0.35),objfun,n=in_sample,method="L-BFGS-B",lower = c(0.01, 0.01, 0.01), upper=c(0.99,0.99,0.99))    
  rez
}

all_files <- list.files(pattern = "\\.csv$")
ds_name <- c("HD07_N_A79HP",
             "HD09_N_A93F3",
             "HD09_N_AB8KB",
             "HD10_N_AB8KB",
             "HD13_N_A93F3",
             "HD13_N_AB8KB")

plist = list()
counter <- 0
for (f in all_files){
  counter <- counter + 1
  dataset <- read.csv(f, header = TRUE) 
  plist[[counter]] <- dataset
}

rez_df = data.frame(sample_size = c(), RSSE = c())

for(sample_size in c(10,100,1000,2000,4000,8000,16000)){
  er <- sapply(c(1:160), function(k){
    ds_num <- sample(x = c(2:length(plist)), 1, replace = T, prob = rep(1,5)/5)
    dataset <- plist[[ds_num]]
    subsample <- sample(x = dataset$Length, sample_size, replace = T, prob = dataset$Freq_norm)
    rez <- my_optimsearch(subsample)
#    print(rez$par)    
    emp <- sapply(c(1:50), function(x) sum(subsample==x))
    emp <- emp/sample_size
    sqrt(sum((my_dt(rez$par,c(1:50)) - emp)^2))
  })
  new_df <- data.frame(sample_size = rep(sample_size,40), RSSE = er)
  rez_df <- rbind(rez_df,new_df)
  print(er)
  print(mean(er))
}

rez_df$sample_size <- as.factor(rez_df$sample_size)
rez_df$logRSSE = log(rez_df$RSSE)
ggplot(rez_df, aes(x=sample_size, y=logRSSE, fill=sample_size)) + 
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        plot.title = element_text(size=18,face="bold"),
        legend.position="none") +
  labs(x="Repertoire size", y = "log(RSSE)")

means <- sapply( c(10,100,1000,2000,4000,8000,16000), 
                 function(x) mean(rez_df[rez_df$sample_size==x,]$logRSSE))

sds <- sapply( c(10,100,1000,2000,4000,8000,16000), 
                 function(x) sd(rez_df[rez_df$sample_size==x,]$logRSSE))

appl_df <- data.frame(sample_size = log10(c(10,100,1000,2000,4000,8000,16000)),
                      mean = means,
                      ub = means + 1.64*sds)

experim <- read.csv('/home/stephen/Git_Repositories/VDJ_Recombination/Model/experimental_data.csv')
mod <- read.csv('/home/stephen/Git_Repositories/VDJ_Recombination/Model/model_data.csv')



error <- function(vector1,vector2) log( sqrt( sum( (as.numeric(vector2) - as.numeric(vector1))^2 ) ) )
error_notlog <- function(vector1,vector2)  sqrt( sum( (as.numeric(vector2) - as.numeric(vector1))^2 ) ) 

#############

ID = c('SRR5888724','SRR5888725','SRR5888726',
       'SRR5888727','SRR5888728','SRR5888729',
       'SRR5888730','SRR5888731','SRR5888732',
       'SRR5888733','SRR5888734','SRR5888735',
       'SRR5888736')
Week <- c(0,14,15,19,31,43,55,80,102,128,153,177,200)
Error <- mapply(error_notlog, experim, mod)
Significant = c('-','-','-',
       '-','-','+',
       '-','-','+',
       '+','-','-',
       '+')
df <- data.frame(ID, Week, Error, Significant)
bla <- capture.output(stargazer(df, out='table.html'))



HIV_x <- log10(c(4291, 4041,3572,4679,4480,3826,3856,3369,3673,3277,3603,3200,2727))
HIV_y <- mapply(error, experim, mod)

ggplot(appl_df) +
  geom_area(aes(x=sample_size, y=ub),size =1.5,color = "black",fill = "grey") +
  geom_line(aes(x=sample_size, y=mean),size =1.5,linetype = "dashed") +
  geom_point(aes(x = HIV_x[1], y = HIV_y[1], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[2], y = HIV_y[2], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[3], y = HIV_y[3], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[4], y = HIV_y[4], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[5], y = HIV_y[5], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[6], y = HIV_y[6], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[7], y = HIV_y[7], colour = "HIV"),shape=17,size=4) + #interesting
  geom_point(aes(x = HIV_x[8], y = HIV_y[8], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[9], y = HIV_y[9], colour = "HIV"),shape=17,size=4) + #interesting
  geom_point(aes(x = HIV_x[10], y = HIV_y[10], colour = "HIV"),shape=17,size=4) + #interesting
  geom_point(aes(x = HIV_x[11], y = HIV_y[11], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[12], y = HIV_y[12], colour = "HIV"),shape=17,size=4) +
  geom_point(aes(x = HIV_x[13], y = HIV_y[13], colour = "HIV"),shape=17,size=4) + #interesting
  scale_colour_manual(values=c("black", "black", "black", "black")) +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        plot.title = element_text(size=18,face="bold"),
        legend.position = "bottom") +
  scale_x_continuous(breaks=log10(c(10,100,1000,10000)),
                  labels=c("10", "100", "1000", "10000")) +
  labs(x="Sample size", y = "log(RSSE)")

ggplot(df, aes(Significant, Error)) + geom_boxplot() + theme_bw()