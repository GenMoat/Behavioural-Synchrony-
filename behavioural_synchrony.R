#Install required packages
install.packages('ggplot2')
install.packages('tidyverse')
library(ggplot2)
library(ggridges)
library(tidyverse)
library(readr)

# Load data
library(readr)
df <- read_csv("~/Desktop/data.csv") #please change to your path


#Set seed so we have the same results each time we run randomization 
#please change this to a random number 
set.seed(123)



# Data visualization
head(df)


#Extract sub-dataset only contain the behaviour information of two species. 
#(This is the format required to make this script works)

#column name for Capuchins' behaviour 
c_beh_col <- c('VC','FC','LC','RC','PC') 
#column name for Squirrel monkey's behaviour
s_beh_col <- c('VS','FS','LS','RS','PS') 

#extract
df_beh <- df[,append(c_beh_col,s_beh_col)]
head(df_beh)

#Plot the frequency of the behaviour in the dataset

no_beh=5 #number of behaviours observed

# the number of columns of freq_beh should be even, such that it includes
# the behaviours of two species.
freq_beh <- data.frame(t(colSums(df_beh))) # column sum of the dataset
#create data frames for plotting purpose  

beh_col <- c("V", "F","L","R","P") # name of behaviours
fbl <- unname(unlist(freq_beh[1,])) #convert freq_beh to list


freq_total<- vector()
for (i in c(1:no_beh)){
  freq_total<- append(freq_total,c(fbl[[i+no_beh]],fbl[[i]], 
                                   sum(fbl[[i]],fbl[[i+no_beh]])))
} 
#freq_total=(# of beh1 of c,# of beh1 of s, total beh1, etc for every behaviour)


# reordered the number of squirrel monkey and capuchins to make the plot nicer
freq_beh_4_plot_t <- data.frame(behaviour=rep(beh_col, each=3),
                                species=rep(c("Squirrel monkey",
                                "Capuchins", "Combined"),5),
                                frequency=c(freq_total))


# fix the order of bar display in the ggplot 
freq_beh_4_plot_t$behaviour <- factor(freq_beh_4_plot_t$behaviour, 
                                      levels = unique(freq_beh_4_plot_t$behaviour))
freq_beh_4_plot_t$species <- factor(freq_beh_4_plot_t$species, 
                                    levels = unique(freq_beh_4_plot_t$species))


#plot the frequency of bhevaiours (and the total) 
ggplot(data=freq_beh_4_plot_t, aes(x=behaviour, y=frequency, fill=species)) +
  labs(x="Behaviour", y="Frequency", fill="Species")+
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=frequency), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

#The (observed) probability of occurrence of each behaviour.
prob_c<-vector() # capuchins
prob_s<-vector() # squirrel monkey
for (i in c(1:no_beh)){
  total_c <- sum(freq_beh[1:no_beh])
  total_s <- sum(freq_beh[(no_beh+1):(2*no_beh)])
  
  prob_c<- append(prob_c,freq_beh[[i]]/total_c)
  prob_s<- append(prob_s,freq_beh[[i+no_beh]]/total_s)
}

prob_c
prob_s

################################################################################
# Functions
## This function create randomization of the Mixed group:
random_gen<- function(list_x,list_y,x,y){
  # list_x = list of behaviour data of group 1
  # list_x = list of behaviour data of group 2
  # x = number of monkey in group 1
  # y = numbber of monkey in group 2
  i<-0
  cat_list=list_x+list_y # combine the # of behaviours in two groups to one dataset
  need_shuffle <- vector() # create a empty vector to catch the results
  for (beh in cat_list){
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
  # data = data frame with only the number of behaviour occurs for each scan
  
  # g1 = index of the behaviour colomn starts for species 1
  # g2 = index of the behaviour colomn starts for species 2
  # itr = how many iteration for each row's randomization
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
    
    
    #create dummy dataframe to deal with beh not present
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
  #return(df_temp)
  cor_matrix <- round(cor(df_temp, method = "pearson"),4)
  #print(cor_matrix)
  cor_matrix<-as.data.frame(as.table(cor_matrix))
  cor_matrix<-na.omit(cor_matrix[2*no_beh*no_beh+1:nrow(cor_matrix),])#2*5*5=50 
  cor = cor_matrix[seq(1, nrow(cor_matrix), 2*no_beh+1), ] # this gives the correlation 
  rownames(cor) = 1:nrow(cor)
  return(cor)#return the final corr for this randomized dataset
  
}



#function for collect correlations
collect_corr <- function(df,iter,no_beh){ # function to collect the pearson correlation 
  # df = dataset
  # iter = iteration/number of correlation to collect
  # no_beh = number of behavious measured in the dataset
  corr <- beh_randomization_cross_corr(df_beh,5)
  
  for (i in 1:(iter-1)){
    Freq <- beh_randomization_cross_corr(df_beh,5)[,3]
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
#The actural correlation (from the observed data)
cor_matrix_act <-  round(cor(df[,c(7:11,13:17)], method = "pearson"),4) # correlation matrix
cor_matrix_act<-as.data.frame(as.table(cor_matrix_act))
cor_act<-na.omit(cor_matrix_act[2*no_beh*no_beh+1:nrow(cor_matrix_act),])
cor_act <- cor_act[seq(1, nrow(cor_act), 2*no_beh+1), ] # correlation scores

summary(cor_act[,3])


#Plot the observed cross-correlations
#every possoble pair of the two list of behaviours of two species
corr_temp <- crossing(Var1 = s_beh_col, Var2 = c_beh_col)
corr_plot_df <- cor_matrix_act[(cor_matrix_act$Var1 %in% corr_temp$Var1) & (cor_matrix_act$Var2 %in% corr_temp$Var2) ,]

ggplot(corr_plot_df, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile()+
  labs(x="Behaviours of capuchins",y="Behaviours of squrriel monkeys", fill="Correlation")+
  geom_text(aes(label = round(Freq,3)), color = "white", size = 4)+ 
  scale_fill_gradient(low = "steelblue1",
                      high = "midnightblue",
                      guide = "colorbar")+
  theme(legend.title=element_blank())+ 
  geom_segment(aes(x = 1.5, xend = 1.5, y = .5, yend = 1.5), colour = "green", size = 1) +
  geom_segment(aes(x = 2.5, xend = 2.5, y = 1.5, yend = 2.5), colour = "green", size = 1) +
  geom_segment(aes(x = 3.5, xend = 3.5, y = 2.5, yend = 3.5), colour = "green", size = 1) +
  geom_segment(aes(x = 4.5, xend = 4.5, y = 3.5, yend = 4.5), colour = "green", size = 1) +
  geom_segment(aes(x = 5.5, xend = 5.5, y = 4.5, yend = 5.5), colour = "green", size = 1) +
  geom_segment(aes(x = .5, xend = 1.5, y = 1.5, yend = 1.5), colour = "green", size = 1) +
  geom_segment(aes(x = 1.5, xend = 2.5, y = 2.5, yend = 2.5), colour = "green", size = 1) +
  geom_segment(aes(x = 2.5, xend = 3.5, y = 3.5, yend = 3.5), colour = "green", size = 1) +
  geom_segment(aes(x = 3.5, xend = 4.5, y = 4.5, yend = 4.5), colour = "green", size = 1) +
  geom_segment(aes(x = 4.5, xend = 5.5, y = 5.5, yend = 5.5), colour = "green", size = 1) +
  theme_minimal()


################################################################################
#Run the randomization for 10,000 times to generate the cross-correlations 
system.time(corr_df<-collect_corr(df,1000,5))

#Summary of the simulated cross-correlations
corr_mean <- colMeans(corr_df[,-c(1,2)]) # mean of correlation of each simulation 
corr_mean<-data.frame(corr_mean)
rownames(corr_mean) <- 1:nrow(corr_mean) # mean corr of each simulation
summary(corr_mean$corr_mean)


corr_df$row_mean <- rowMeans(corr_df[,3:ncol(corr_df)])
corr_for_beh <- corr_df[,c("Var1","Var2","row_mean")]



#Plot the simulated cross-correlations
ggplot(corr_for_beh, aes(x = Var1, y = Var2, fill = row_mean)) +
  geom_tile()+
  #geom_tile(color = "black") +
  geom_text(aes(label = round(row_mean,3)), color = "white", size = 4) +
  scale_fill_gradient(low = "steelblue1",
                      high = "midnightblue",
                      guide = "colorbar")+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.3, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))



#Visualize the simulated cross-correlations
ggplot(corr_mean,aes(corr_mean))+geom_area(aes(y = ..count..),
                                           color="blue", 
                                           fill="lightsteelblue2", 
                                           stat = "bin") 


#Visualize the actual cross-correlations
ggplot(corr_mean,aes(corr_mean))+
  geom_area(aes(y = ..count..),color="blue",
            fill="lightsteelblue2", 
            stat = "bin")+ 
  geom_vline(xintercept = 0.10914, linetype="dotted", size = 0.3,colour="red")+
  geom_text(aes(x=0.11914, 
                label="Observed cross-corr=0.10914",y=100), 
            colour="red", angle=90)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#The P-value for the cross-correlation
print(paste0("The P-value is:",p_value_cal(corr_mean,1,0.10914,'left')))



################################################################################
## Compute Bsyn score within the species

### Capuchins
#The actual Bsyn score for Capuchins
bs_c<-vector()
for (i in 1:nrow(df)){
  bs_temp <- simpson(df[i, 7:11])
  bs_c <-append(bs_c, bs_temp)
}
df$Bsyn_c<-as.numeric(unlist(bs_c)) #add to dataframe
summary(df$Bsyn_c)

bs_act_c <- mean(df$Bsyn_c) # the actual mean of the Bsyn score for capuchins


#Run the ramdomization test 10,000 times for capuchins
system.time(Bsync_means_c<-replicate(1000, 
                                     bsyn_scores_single(df_beh,1,5,prob_c)))


#Summary of the simulated mean of Bsyn score for capuchins
summary(Bsync_means_c)

#Plot the result (simulated Bsyn scores)
Bsync_means_c_df <- as.data.frame(Bsync_means_c)
ggplot(Bsync_means_c_df,aes(Bsync_means_c))+geom_area(aes(y = ..count..),
                                                      color="blue", 
                                                      fill="lightsteelblue2", 
                                                      stat = "bin") 


#Plot the actual Bsyn score with the simulated data
ggplot(Bsync_means_c_df,aes(Bsync_means_c))+geom_area(aes(y = ..count..),
                                                      color="blue", 
                                                      fill="lightsteelblue2", 
                                                      stat = "bin")+
  geom_vline(xintercept =bs_act_c, linetype="dotted", size = 0.3,colour="red")+
  geom_text(aes(x=bs_act_c, label="Observed Bsyn score= 0.3983",y=100), 
            colour="red", angle=90)


#calculate the p-value
print(paste0("The P-value is:",p_value_cal(Bsync_means_c_df,1,bs_act_c,'right'))) 
# you might want to change 'left' or 'right'




### Squirrel monkey

#Run the ramdomization test 10,000 times for Squirrel monkey
system.time(Bsync_means_s<-replicate(1000, 
                                     bsyn_scores_single(df_beh,6,5,prob_s))) 

#The summary of the simulated dataset for squirrel monkey
summary(Bsync_means_s)

#Plot the result (simulated Bsyn scores)
Bsync_means_s_df <- as.data.frame(Bsync_means_s)
ggplot(Bsync_means_s_df,aes(Bsync_means_s))+geom_area(aes(y = ..count..),
                                                      color="blue", 
                                                      fill="lightsteelblue2", 
                                                      stat = "bin") 


#The actual Bsyn score for Squirrel monkey
bs_s<-vector()
for (i in 1:nrow(df)){
  bs_temp <- simpson(df[i, 13:17])
  bs_s <-c(bs_s, bs_temp)
}

df$Bsyn_s<-as.numeric(unlist(bs_s)) #add to dataframe
summary(df$Bsyn_s)

bs_act_s <- mean(df$Bsyn_s) # the actual mean of the Bsyn score for capuchins

#Plot the actual Bsyn score with the simulated data
ggplot(Bsync_means_s_df,aes(Bsync_means_s))+geom_area(aes(y = ..count..),
                                                      color="blue", 
                                                      fill="lightsteelblue2", 
                                                      stat = "bin")+
  geom_vline(xintercept =bs_act_s, linetype="dotted", size = 0.3,colour="red")+
  geom_text(aes(x=bs_act_s, label="Observed Bsyn score=0.5058",y=100), 
            colour="red", angle=90)


#calculate the p-value
print(paste0("The P-value is:",p_value_cal(Bsync_means_s_df,1,bs_act_s,'right')))


#Combine the two distributions into one plot
Bsync_means_c_df$species <- rep("capuchin",10)
Bsync_means_s_df$species <- rep("squirrel monkey",10)
colnames(Bsync_means_c_df) <- c("Bsyn_score","species")
colnames(Bsync_means_s_df) <- c("Bsyn_score","species")
#combine the two distributions
combined_cs <- rbind(Bsync_means_c_df, Bsync_means_s_df)
#data frame with the actual mean Bsyn scores
draw_act_mean <- data.frame(act_mean = c(bs_act_s, bs_act_c))
dummy_data<-data.frame(x = 0,y=0.6)

ggplot(combined_cs, aes(x =Bsyn_score , y = species)) + 
  geom_density_ridges(aes(fill = species))+
  labs(x="Bsyn score",y="Species")+
  scale_fill_manual(values = c("lightsteelblue1", "lightsteelblue3")) +
  geom_segment(data = draw_act_mean, aes(x=act_mean,
                                         xend=act_mean, 
                                         y=c(2,1), 
                                         yend=c(2.5,1.5)),
               color=c("royalblue4", "royalblue1"))+
  geom_segment(data=dummy_data,aes(x=0.468,xend=0.506,y=2.01,yend=2.01))+
  geom_text(data=draw_act_mean, aes(x=act_mean,xend=act_mean, 
                                    y=c(2.6,1.6), 
                                    yend=c(2.6,1.6), 
                                    label=c("0.506","0.398")))+
  scale_y_discrete(expand = c(0, 0))+
  theme_ridges(grid = FALSE, center = TRUE)+
  theme(legend.position="none")



################################################################################
#Bsyn score for combined data
#create a dummy data frame
df_combined <- data.frame(matrix(ncol = no_beh,nrow = 0)) 
beh <- c("V","F","L","R","P") 
colnames(df_combined)<- beh
#combine the original data 
for (i in 1:nrow(df_beh)){
  for (j in 1:no_beh){
    df_combined[i,j] <- sum(df_beh[i,j],df_beh[i,j+5])
  }
}

# calculate Bsyn score for the combine group 
bs_com<- vector()
for (i in 1:nrow(df_combined)){
  bs_tem <- simpson(df_combined[i,1:no_beh])
  bs_com <- append(bs_com,bs_tem)
}
# the mean Bsyn score of the observed data (combined two species)
mean(as.numeric(unlist(bs_com)))














































