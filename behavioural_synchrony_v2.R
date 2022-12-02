#Install required packages
#install.packages('ggplot2')
#install.packages('tidyverse')
library(ggplot2)
library(ggridges)
library(tidyverse)
library(readr)
library(psych)
# Load data
library(readr)
df <- read_csv("~/LLs_MSG_BSynch2022.csv") #please change to your path

#OPTIONAL - change the column name
colnames(df)[10] <-"RC"
colnames(df)[16] <-"RS"

#Set seed so we have the same results each time we run randomization 
set.seed(123)

#extract info (columns) needed
c_ew_beh_col <- df[c('Group (1 = West; 2 = East)','VC','FC','LC','RC','PC')]
s_ew_beh_col <-df[c('Group (1 = West; 2 = East)','VS','FS','LS','RS','PS')]


################################################################################
# Functions
## This function create randomization of the Mixed group:
random_gen<- function(list_x,list_y,x,y){
  # list_x = list of behaviour data of group 1
  # list_x = list of behaviour data of group 2
  # x = number of monkey in group 1
  # y = numbber of monkey in group 2
  i<-0
  cat_list=list_x+list_y # combine the number of behaviours 
                         #in two groups to one dataset
  
  need_shuffle <- vector() # create a empty vector to catch the results
  for (beh in cat_list){ #this loop shuffles each rows of the dataset
    i = i+1
    if (beh!=0){
      need_shuffle= append(need_shuffle,rep(i,beh))
    }
    sf=sample(need_shuffle,length(need_shuffle),replace=FALSE)
    new_list_x=sf[1:x]
    new_list_y=sf[x+1:y]
  }
  return(list(new_list_x,new_list_y))
}


## This function computes Simpson diversity score
simpson <- function(list){ 
  # list = list of frequency of the behaviour
  syn=0
  for (i in 1:length(list)){
    sum_list=sum(list)
    syn_temp= (list[i]*(list[i]-1))/(sum_list*(sum_list-1))
    syn=syn+syn_temp
  }
  return(syn)
}


## This function turns row of a data frame to a list
row_collector <- function(df,row_idx,col_idx1,col_idx2){
  row=unname(unlist(df[row_idx,(col_idx1):(col_idx2)]))
  return(row)
}


################################################################################
# Analysis
## Cross-correlation randomization
beh_randomization_cross_corr <- function(df,no_beh){
  # df = data frame with only the number of behaviour occurs for each scan
  # no_beh = number of behaviours in the dataframe
  
  # create a data frame to load the new random values
  df_temp <- data.frame(matrix(ncol = 2*no_beh,nrow = 0)) 
  beh_x <- sprintf("beh_x%d", 1:no_beh) 
  beh_y <- sprintf("beh_y%d", 1:no_beh) 
  colnames(df_temp)<- c(beh_x,beh_y) # empty df created
  
  
  for (i in 1:nrow(df)){# loop through each row of the df
    
    # perform the randomization with the i-th row
    species_1 <- row_collector(df,i,1,no_beh)
    species_2 <-row_collector(df,i,no_beh+1,2*no_beh)
    num_list<- random_gen(species_1,species_2, sum(species_1),sum(species_2))
    
    # beh (e.g. 1=A,2=B etc. back to number of occurrence)
    list_x<-as.data.frame(table(unlist(num_list[[1]])))
    list_y<-as.data.frame(table(unlist(num_list[[2]])))
    
    
    #create dummy dataframe to deal with behaviour not present
    dummy <- data.frame (Var1  = c(seq(1:no_beh)), 
                         Freq = c(rep(0,no_beh)))
    
    list_x_1<-merge(x = dummy, y = list_x, by = "Var1", all.x = TRUE)
    list_y_1<-merge(x = dummy, y = list_y, by = "Var1", all.x = TRUE)
    
    list_x_1<-within(list_x_1, rm("Freq.x")) #remove the extra column
    list_y_1<-within(list_y_1, rm("Freq.x")) #remove the extra column 
    list_x_1[is.na(list_x_1)] = 0 #NA to 0
    list_y_1[is.na(list_y_1)] = 0
    list_x_1 <-list_x_1$Freq.y
    list_y_1 <-list_y_1$Freq.y
    #now we have length of list for each group is 5
    df_temp[i,1:no_beh] <- unlist(list_x_1) # add the list to the dataframe
    df_temp[i,(no_beh+1):(2*no_beh)] <- unlist(list_y_1)
    
  }
  cor_matrix <- round(cor(df_temp, method = "pearson"),4)
  cor_matrix<-as.data.frame(as.table(cor_matrix))
  cor_matrix<-na.omit(cor_matrix[2*no_beh*no_beh+1:nrow(cor_matrix),])#2*5*5=50 
  # this gives the correlation
  cor = cor_matrix[seq(1, nrow(cor_matrix), 2*no_beh+1), ]  
  rownames(cor) = 1:nrow(cor)
  return(cor)#return the final corr for this randomized dataset
  
}



#function for collect correlations
collect_corr <- function(df,iter,no_beh){
  # df = dataset
  # iter = iteration/number of correlation to collect
  # no_beh = number of behavious measured in the dataset
  corr <- beh_randomization_cross_corr(df,5)
  for (i in 1:(iter-1)){
    Freq <- beh_randomization_cross_corr(df,5)[,3]
    print(Freq)
    corr[,ncol(corr)+1] <- Freq # add the correlation results to 'corr'
    colnames(corr)[ncol(corr)] <- paste0("Freq", i+1)
    
  }
  return(corr)
}



#This function calculates the P-value
p_value_cal <- function(data,col_idx,act_value,side){
  # data= data.frame format data
  # act_value = a single value (float)
  # side = is the actual data on the left or the right hand side of 
  #        the distribution
  
  len <- length(data[,col_idx]) # the number of similated data
  if (side=='left'){
    total <- sum(data[,col_idx]>act_value) #count the number of the data greater
    #than the actual value
  } else if (side=='right'){
    total <- sum(data[,col_idx]<act_value) #count the number of the data greater
    #than the actual value
  }
  p_value=(len-total)/len
  return(p_value)
}



#This function is for running randomization for each row of data to answer the
#within species synchronization problem.
bsyn_scores_single <- function(df,g1,no_beh,prob_list){
  # df = data frame
  # g1 = index of the behaviour colomn starts 
  # no_beh = number of behaviours in the dataset
  # prob_list = probabilitie of the occurrence of each behaviour
  
  # create a data frame to load the new random values
  df_temp <- data.frame(matrix(ncol = no_beh,nrow = 0)) 
  beh_x <- sprintf("beh_x%d", 1:no_beh) 
  colnames(df_temp)<- c(beh_x) # empty df create here
  
  syn=vector() #empty vector to collect the sun scores
  for (i in 1:nrow(df)){# loop through each row of the df
    row_total=sum(row_collector(df,i,g1,g1+no_beh-1))# total no of group x 
    
    # perform randomization within the species, 
    # based on the probability of occurrence
    samp_temp<-sample(beh_x, size = row_total, replace=TRUE, prob= prob_list)
    tab_temp<-table(samp_temp) #add row to dataset
    tab_temp<-as.data.frame(tab_temp)
    colnames(tab_temp)<-c("Var1", "Freq")
    
    # create dataframe to account for rows without behaviour (nulls)
    dummy <- data.frame (Var1  = beh_x, Freq = c(rep(0,5))) 
    #merge dataframes
    tab_temp<-merge(x = dummy, y = tab_temp, by = "Var1", all.x = TRUE)
    tab_temp$Freq<-rowSums(tab_temp[,c("Freq.x", "Freq.y")], na.rm=TRUE)
    #tab_temp<-within(tab_temp, rm("Freq.x"))
    #tab_temp<-within(tab_temp, rm("Freq.y"))
    df_temp[i,1:no_beh] <- unlist(tab_temp$Freq)
    
    syn_score <- simpson(unlist(df_temp[i,], use.names = FALSE))
    syn[i] <- syn_score
    
  }
  return(mean(syn)) 
}
################################################################################
#The actual SDI score for Capuchin (West+East)
bs_c<-vector()
for (i in 1:nrow(df)){
  bs_temp <- simpson(df[i, 7:11])
  bs_c <-append(bs_c, bs_temp)
}
df$Bsyn_c<-as.numeric(unlist(bs_c)) #add to dataframe
summary(df$Bsyn_c)

#The actual SDI score for Squirrel monkey (West+East)
bs_s<-vector()
for (i in 1:nrow(df)){
  bs_temp <- simpson(df[i, 13:17])
  bs_s <-c(bs_s, bs_temp)
}

df$Bsyn_s<-as.numeric(unlist(bs_s)) #add to dataframe
summary(df$Bsyn_s)

#separate SDI scores for each EC, ES and WC, WS
##WS
mean(df[df$`Group (1 = West; 2 = East)`==1,]$Bsyn_s)
##WC
mean(df[df$`Group (1 = West; 2 = East)`==1,]$Bsyn_c)
##ES
mean(df[df$`Group (1 = West; 2 = East)`==2,]$Bsyn_s)
##EC
mean(df[df$`Group (1 = West; 2 = East)`==2,]$Bsyn_c)

################################################################################
no_beh=5
# The probabilities of each behaviour
#The (observed) probability of occurrence of each behaviour.

##CW=data frame only contains the data of West Capuchins.
cw<-c_ew_beh_col[c_ew_beh_col$`Group (1 = West; 2 = East)`==1,]
cw_prob<-unlist(colSums(as.data.frame(lapply(cw[,c(2:6)], 
                                  function(x) x / sum(colSums(cw[,c(2:6)]))))))

##CE
ce<-c_ew_beh_col[c_ew_beh_col$`Group (1 = West; 2 = East)`==2,]
ce_prob<-unlist(colSums(as.data.frame(lapply(ce[,c(2:6)], 
                                   function(x) x / sum(colSums(ce[,c(2:6)]))))))

##SW
sw<-s_ew_beh_col[s_ew_beh_col$`Group (1 = West; 2 = East)`==1,]
sw_prob<-unlist(colSums(as.data.frame(lapply(sw[,c(2:6)], 
                                   function(x) x / sum(colSums(sw[,c(2:6)]))))))

##SE
se<-s_ew_beh_col[s_ew_beh_col$`Group (1 = West; 2 = East)`==2,]
se_prob<-unlist(colSums(as.data.frame(lapply(se[,c(2:6)], 
                                   function(x) x / sum(colSums(se[,c(2:6)]))))))

################################################################################

#The actual SDI score and add to dataset
##CW
bs_cw<-vector()
for (i in 1:nrow(cw)){
  bs_temp <- simpson(cw[i, 2:6])
  bs_cw <-append(bs_cw, bs_temp)
}
cw$Bsyn<-as.numeric(unlist(bs_cw)) #add to dataframe
mean(cw$Bsyn)

##CE
bs_ce<-vector()
for (i in 1:nrow(ce)){
  bs_temp <- simpson(ce[i, 2:6])
  bs_ce <-append(bs_ce, bs_temp)
}
ce$Bsyn<-as.numeric(unlist(bs_ce)) #add to dataframe
mean(ce$Bsyn)

##SW
bs_sw<-vector()
for (i in 1:nrow(sw)){
  bs_temp <- simpson(sw[i, 2:6])
  bs_sw <-append(bs_sw, bs_temp)
}
sw$Bsyn<-as.numeric(unlist(bs_sw)) #add to dataframe
mean(sw$Bsyn)

##SE

bs_se<-vector()
for (i in 1:nrow(se)){
  bs_temp <- simpson(se[i, 2:6])
  bs_se <-append(bs_se, bs_temp)
}
se$Bsyn<-as.numeric(unlist(bs_se)) #add to dataframe
mean(se$Bsyn)

################################################################################

#Simulate Bsyn scores
##CW
system.time(Bsync_means_cw<-replicate(1000, #1000
                                      bsyn_scores_single(cw,2,5,cw_prob)))
##CE
system.time(Bsync_means_ce<-replicate(1000, #1000
                                      bsyn_scores_single(ce,2,5,ce_prob)))
##SW
system.time(Bsync_means_sw<-replicate(1000, #1000
                                      bsyn_scores_single(sw,2,5,sw_prob)))
##SE
system.time(Bsync_means_se<-replicate(1000, #1000
                                      bsyn_scores_single(se,2,5,se_prob)))



################################################################################
# The probabilities of each behaviour
#The (observed) probability of occurrence of each behaviour.

## West (Capuchin + Squirrel monkey)
west_cs <- data.frame(matrix(ncol = no_beh,nrow = 90)) 
beh <- c("V","F","L","RC","P") 
colnames(west_cs)<- beh

west_cs$V<- cw$VC+sw$VS
west_cs$F<- cw$FC+sw$FS
west_cs$L<- cw$LC+sw$LS
west_cs$RGS<- cw$RC+sw$RGSS
west_cs$P<- cw$PC+sw$PS

#Probabilities
west_cs_prob<-unlist(colSums(as.data.frame(lapply(west_cs, 
                                       function(x) x / sum(colSums(west_cs))))))


## East (Capuchin + Squirrel monkey)
east_cs <- data.frame(matrix(ncol = no_beh,nrow = 90)) 
colnames(east_cs)<- beh

east_cs$V<- ce$VC+se$VS
east_cs$F<- ce$FC+se$FS
east_cs$L<- ce$LC+se$LS
east_cs$RGS<-ce$RGSC+se$RS
east_cs$P<- ce$PC+se$PS

#Probabilities
east_cs_prob<-unlist(colSums(as.data.frame(lapply(east_cs, 
                                       function(x) x / sum(colSums(east_cs))))))
################################################################################

#The actual Bsyn score for Capuchin and Squirrel monkey (2 species) combined
##west
bs_west<-vector()
for (i in 1:nrow(west_cs)){
  bs_temp <- simpson(west_cs[i, 1:5])
  bs_west <-append(bs_west, bs_temp)
}
west_cs$Bsyn<-as.numeric(unlist(bs_west)) #add to dataframe
mean(west_cs$Bsyn)

##east
bs_east<-vector()
for (i in 1:nrow(east_cs)){
  bs_temp <- simpson(east_cs[i, 1:5])
  bs_east <-append(bs_east, bs_temp)
}
east_cs$Bsyn<-as.numeric(unlist(bs_east)) #add to dataframe
mean(east_cs$Bsyn)


#Simulate bsyn score for C+S (2 species) combined
##West
system.time(Bsync_means_cs_west<-replicate(1000, #1000
                                  bsyn_scores_single(west_cs,1,5,west_cs_prob)))

##East
system.time(Bsync_means_cs_east<-replicate(1000, #1000
                                  bsyn_scores_single(east_cs,1,5,east_cs_prob)))

################################################################################
#Actual Correlation
##West - correlation matrix
corr_west <-  round(cor(cw[2:6],sw[2:6], method = "pearson"),4) 
corPlot(corr_west)

##East - correlation matrix
cor_east <-  round(cor(ce[2:6],se[2:6], method = "pearson"),4) 
corPlot(cor_east)

##Capuchins - correlation matrix
corr_c <-  round(cor(cw[2:6],ce[2:6], method = "pearson"),4) 
corPlot(corr_c)

##Squirrel monkey - correlation matrix
cor_s <-  round(cor(sw[2:6],se[2:6], method = "pearson"),4) 
corPlot(cor_s)
################################################################################
#simulated correlation
#West vs east
west<-cbind(cw[2:6], sw[2:6])
east<-cbind(ce[2:6], se[2:6])

est_corr_west<-collect_corr(west,1000,5)
est_corr_east<-collect_corr(east,1000,5)

# mean of correlation of each simulation
est_corr_west_mean <- rowMeans(est_corr_west[,-c(1,2)]) 
est_corr_west_mean<-data.frame(est_corr_west_mean)

# mean of correlation of each simulation
est_corr_east_mean <- rowMeans(est_corr_east[,-c(1,2)]) 
est_corr_east_mean<-data.frame(est_corr_east_mean)

west_corr<-as.data.frame(t(est_corr_west[,-c(1,2)]))
east_corr<-as.data.frame(t(est_corr_east[,-c(1,2)]))

#Capuchins vs Squirrel monkey

capuchins<-cbind(cw[2:6], ce[2:6])
squirrel_monkey<-cbind(sw[2:6], se[2:6])

est_corr_c<-collect_corr(capuchins,1000,5)
est_corr_s<-collect_corr(squirrel_monkey,1000,5)##TODO

# mean of correlation of each simulation 
est_corr_c_mean <- rowMeans(est_corr_c[,-c(1,2)]) 
est_corr_c_mean<-data.frame(est_corr_c_mean)

# mean of correlation of each simulation 
est_corr_s_mean <- rowMeans(est_corr_s[,-c(1,2)]) 
est_corr_s_mean<-data.frame(est_corr_s_mean)

c_corr<-as.data.frame(t(est_corr_c[,-c(1,2)]))

################################################################################

#Correlation between number of behaviours and SDI  - correlation matrix
corr_bsdi_cw <-  round(cor(rowSums(cw[2:6]),cw[7], method = "pearson"),4)


################################################################################
#Correlation plot 
cor_matrix_act<-as.data.frame(as.table(cor_east))

ggplot(cor_matrix_act, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile()+
  labs(x="Behaviours of squrriel monkeys (east)",
       y="Behaviours of capuchins (east)", 
       fill="Correlation")+
  geom_text(aes(label = round(Freq,3)), color = "black", size = 4)+ 
  scale_fill_gradient2(low = "#ca0020",
                       mid="#f7f7f7",
                       high = "#0571b0",
                       midpoint = 0,
                       guide = "colorbar")+
  theme(legend.title=element_blank())+ 
  geom_segment(aes(x = 1.5, xend = 1.5, y = .5, yend = 1.5), colour = "black", size = 1) +
  geom_segment(aes(x = 2.5, xend = 2.5, y = 1.5, yend = 2.5), colour = "black", size = 1) +
  geom_segment(aes(x = 3.5, xend = 3.5, y = 2.5, yend = 3.5), colour = "black", size = 1) +
  geom_segment(aes(x = 4.5, xend = 4.5, y = 3.5, yend = 4.5), colour = "black", size = 1) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = 4.5, yend = 5.5), colour = "black", size = 1) +
  geom_segment(aes(x = .5, xend = 1.5, y = 1.5, yend = 1.5), colour = "black", size = 1) +
  geom_segment(aes(x = 1.5, xend = 2.5, y = 2.5, yend = 2.5), colour = "black", size = 1) +
  geom_segment(aes(x = 2.5, xend = 3.5, y = 3.5, yend = 3.5), colour = "black", size = 1) +
  geom_segment(aes(x = 3.5, xend = 4.5, y = 4.5, yend = 4.5), colour = "black", size = 1) +
  geom_segment(aes(x = 4.5, xend = 5.5, y = 5.5, yend = 5.5), colour = "black", size = 1) +
  theme_minimal()

################################################################################
#calculate p-value
#West
est_corr_west_df<-as.data.frame(t(est_corr_west))
print(paste0("The P-value is:",p_value_cal(est_corr_west_df,1,0.2600,'left')))
print(paste0("The P-value is:",p_value_cal(est_corr_west_df,2,-0.0209,'left')))
print(paste0("The P-value is:",p_value_cal(est_corr_west_df,3,0.0866,'left')))
print(paste0("The P-value is:",p_value_cal(est_corr_west_df,4,0.3576 ,'left')))
print(paste0("The P-value is:",p_value_cal(est_corr_west_df,5,0.2072,'left')))
#East
est_corr_east_df<-as.data.frame(t(est_corr_east))
print(paste0("The P-value is:",p_value_cal(est_corr_east_df,1,0.3092,'left')))
print(paste0("The P-value is:",p_value_cal(est_corr_east_df,2,0.0867,'left')))
print(paste0("The P-value is:",p_value_cal(est_corr_east_df,3,0.0973 ,'left')))
print(paste0("The P-value is:",p_value_cal(est_corr_east_df,4,0.0657 ,'left')))
print(paste0("The P-value is:",p_value_cal(est_corr_east_df,5,-0.0979,'left')))




























