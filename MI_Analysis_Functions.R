####Custom R Functions to Accompany MI Decay Paper####
#### Author: Elliot Howard-Spink ####
#### Date Created: 7.4.2022 ####
#### Date Last Edited: 7.4.2022 ####





#### Install and load Packages ####

#Please ensure you have the following packages installed and loaded from the library:
#data.table, tidyr, dplyr, minpack.lm, ggplot2

#install.packages("data.table") - This script was written with package version data.table_1.14.2
#install.packages("tidyr") - This script was written with package version tidyr_1.1.4
#install.packages("dplyr") - This script was written with package version dplyr_1.0.7
#install.packages("minpack.lm") - This script was written with package version minpack.lm_1.2-1
#install.packages("ggplot2") - This script was written with package version ggplot2_3.3.5

library(data.table)
library(tidyr)
library(dplyr)
library(minpack.lm)
library(ggplot2)





#### Import example data ####

#Please ensure that your directory is set to MI_Code_Master before trying to import data. 
#Alternatively, provide an extended file path in read.csv().
getwd() #"/Users/user/Desktop/MI_Code_Master"
#setwd() #If required.

example_data<-read.csv("example_data.csv") #Import data.

head(example_data) #Here we have an example dataset, taken from 1 individual (Fanle). The data features all sequences
#coded for fanle, concatenated together in order of collection. 
#Sequence data exists in the form of a dataframe, with the following columns:
#'X' Index of table
#'Time' Timestamp of action in video
#'Subject' E.g. Fanle
#'Behavior' E.g. grasp, place, drop
#'Modifier.1' Object being used
#'Hammer' Hammer stone number if recorded
#'Anvil' Anvil stone number if recorded
#'Complete' if sequence observation ran from start to finish
#'Nut_crack' if a nut was cracked within the sequence
#'Nut_Number' Number of nuts cracked in the sequence
#'Event' concatenated description of both behaviour, e.g. An action, grasp, and modifier object, Hammer, as grasp2Hammer
#Note: the number next to the modifier has no relation to quantity - it is merely a biproduct of the 
#keyboard key used to code the element in BORIS (i.e. when an action is coded, the number 2 is pressed to log that the
#action relates to a single Hammer, and BORIS records '2HAMMER').





#### MI Estimation ####

#To run MI estimation, the final function Permutation_MI_Adjust is used. This function runs MI estimations
#for sequential distances between 1 and a specified upper limit (Length). The function then adjusts this MI score
#relative to the mean MI scores taken from a number of permuted MI distributions. The number of permutations of the 
#sequence used to generate a permuted MI distribution can be specified by the input PermNo.

#To run MI estimation using Permutation_MI_Adjust, several other custom functions are drawn upon. These include:

#JoinAct - A function which runs through the sequence data (Events) and combines sequence elements which are separated by a specified distance.
#(e.g. Action now and action 1, 2, 3... elements along). The specific comparative distance is set at the start of the function.
#This prepares the data for calculations of joint entropy.

#Mutual_Info_Grass - Grassberger Entropy Estimation and MI Calculation. Calculates Grassberger entropy and MI estimates for
#elements at one specific distance apart (e.g. entropy and MI for sequential elements that are 3 units apart).

#Iterative_MI_Grass- Iterative Grassberger Estimation and MI Calculation - Requires JoinAct & Mutual_Info_Grass; this function
#calculates entropy and MI scores iteratively, for elements over a range of distances (e.g. entropy and MI scores 1,2,3...100 units apart).
#It returns all scores in a single dataframe.

#Permuted_Joint_MI - Requires JoinAct & Mutual_Info_Grass; this function creates a dataframe of MI values estimated using a permuted joint entropy distribution. 
#It permutes the joint entropy distribution once, and runs iterative MI estimations for unit elements in a specified range (e.g. entropy and MI scores 1,2,3...100 units apart).

#Permuted_Joint_MI_Many - Requires JoinAct, Mutual_Info_Grass & Permuted_Joint_MI; this function runs a specified number of permutations of the joint distribution, calculates MI for each in a specific element range
#(e.g. entropy and MI scores 1,2,3...100 units apart), calculates mean for each element range, and returns values in a dataframe.

#Each of these functions can be found below, and tested on their own with the example data.

#TO RUN FULL MI ESTIMATION AND ADJUSTMENT, HAVE ALL FUNCTIONS LOADED IN R (as some are dependent on one another) AND RUN Permutation_MI_Adjust 
#FOR FULL MI ESTIMATION AND PERMUTATION ADJUSTMENT.





#### MI Functions ####

#JoinAct - A function which runs through the sequence data (Events) and combines sequence elements which are separated by a specified distance.
#(e.g. Action now and action 1, 2, 3... elements along). The specific comparative distance is set at the start of the function.
#This prepares the data for calculations of joint entropy.

JoinAct<-function(x,distance){
  e1 <- vector()
  e2 <- vector()
  len<-length(x$Event)
  for(i in (distance+1):len){
    e1[i]<-x$Event[i-distance]
    e2[i]<-x$Event[i]
  }
  
  j <- data.frame(e1,e2)
  j$Joint<-paste(j$e1,j$e2) #Paste the two columns to get a column with both elements in
  names(j)[names(j) == "e1"] <- "ActNow" 
  names(j)[names(j) == "e2"] <- "ActNext"
  j <- na.omit(j) #Remove NA rows
  return(j)
}

#Try with the example data. The comparative distance here is set at 2.
JoinAct(example_data, distance = 2)
#ActNow is the anterior sequence element.
#ActNext is the sequence element 2 terms further down stream in the sequence (this is altered by altering the 'distance' term).
#Joint is the concatenation of these two sequence elements.




#Mutual_Info_Grass - Grassberger Entropy Estimation and MI Calculation. Calculates grassberger entropy and MI estimates for
#elements at one specific distance apart (e.g. entropy and MI for sequential elements that are 3 units apart).
#To run Mutual_Info_Grass, data must first be passed through JoinAct (see example below).
Mutual_Info_Grass<-function(x){
  #X entropy
  N_Now=sum(table(x$ActNow)) #Number of elements in the distribution
  Ni_Now=as.data.frame(table(x$ActNow)) #A dataframe of the frequency each action occurs
  K_Now = length(table(x$ActNow))
  sum_me_Now<-vector(mode = "list", length = K_Now) #Create an empty list to store iterations in before summing.
  for(i in 1:K_Now){
    sum_me_Now[[i]]<-(Ni_Now[i,'Freq']*digamma(Ni_Now[i,'Freq']))
  }
  
  summed_Now<-sum(unlist(sum_me_Now)) #Sum all values produced in iterative section.
  
  Xent<-log2(N_Now) - ((1/N_Now)*summed_Now)
  
  #Y entropy
  N_Next=sum(table(x$ActNext)) #Number of elements in the distribution
  Ni_Next=as.data.frame(table(x$ActNext)) #A dataframe of the frequency each action occurs
  K_Next = length(table(x$ActNext))
  sum_me_Next<-vector(mode = "list", length = K_Next) #Create an empty list to store iterations in before summing.
  for(i in 1:K_Next){
    sum_me_Next[[i]]<-(Ni_Next[i,'Freq']*digamma(Ni_Next[i,'Freq']))
  }
  
  summed_Next<-sum(unlist(sum_me_Next)) #Sum all values produced in iterative section.
  
  Yent<-log2(N_Next) - ((1/N_Next)*summed_Next)
  Yent
  
  #Joint
  N_J=sum(table(x$Joint)) #Number of elements in the distribution
  Ni_J=as.data.frame(table(x$Joint)) #A dataframe of the frequency each action occurs
  K_J = length(table(x$Joint))
  
  sum_me_J<-vector(mode = "list", length = K_J) #Create an empty list to store iterations in before summing.
  for(i in 1:K_J){
    sum_me_J[[i]]<-(Ni_J[i,'Freq']*digamma(Ni_J[i,'Freq']))
  }
  
  summed_J<-sum(unlist(sum_me_J)) #Sum all values produced in iterative section.
  
  Jent<-log2(N_J) - ((1/N_J)*summed_J)
  Jent
  
  MInfo = Xent + Yent - Jent
  df_results<- data.frame(Xent,Yent,Jent,N_J,MInfo)
  return(df_results)
}

#Use JoinAct to set distance over which entropy and MI will be calculated
  #E.g. for element distance of 4 elements.
joint_example_data<-JoinAct(example_data, distance = 4)
  #Then pass that data into Mutual_Info_Grass
Mutual_Info_Grass(joint_example_data)
  #Xent is marginal entropy at X (anterior elements of comparisons).
  #Yent is marginal entropy at Y (posterior elements of comparisons).
  #Jent is joint entropy number of elements in the joint distribution
  #MInfo is MI for that distance.





#Iterative_MI_Grass- Iterative Grassberger Estimation and MI Calculation - Requires JoinAct & Mutual_Info_Grass; this function
#calculates entropy and MI scores iteratively, for elements over a range of distances (e.g. entropy and MI scores 1,2,3...100 units apart).
#It returns all scores in a single dataframe.
#Unlike Mutual_Info_Grass, this function can be run directly from the concatenated dataframe: example_data.

Iterative_MI_Grass<-function(x,Length){  #Length is the max distance you want to make comparisons over
  empty_list <- vector(mode = "list", length = Length)
  for(i in 1:Length){
    empty_list[[i]]<-JoinAct(x, i)
  }
  ret<-empty_list
  l1<-lapply(ret, function(x) x<-Mutual_Info_Grass(x))
  d.info<-rbindlist(l1, idcol = "index")
  return(d.info)
} 

#Example of MI estimations for distances between 1 and 100 units apart (this may take a few seconds).
Iterative_MI_Grass(example_data, Length = 100)





#Permuted_Joint_MI - Requires JoinAct & Mutual_Info_Grass; this function creates a dataframe of MI values estimated using a permuted joint entropy distribution. 
#It permutes the joint entropy distribution once, and runs iterative MI estimations for unit elements in a specified range (e.g. entropy and MI scores 1,2,3...100 units apart).
#Much like in Iterative_MI_Grass, x is the concatenated sequence data, Length is the range of distances over which you would like MI estimated using the permuted sequence data, and seed is the
#seed number used to run a pseudorandom permutation of a sequence (ensuring reproducibility of results).

Permuted_Joint_MI<-function(x,Length,seed){ #x = data with actions and joint actions, Length is maximum distance to
  #compare over, and seed is the seed for pseudorandom sampling.
  Vec<-vector(mode = "list", length = Length)
  for(i in 1:Length){
    Vec[[i]]<-JoinAct(x,i)
  }
  for(i in 1:Length){
    set.seed(seed)
    Vec[[i]]$ActNow_sh<-sample(Vec[[i]]$ActNow)
    Vec[[i]]$ActNext_sh<-sample(Vec[[i]]$ActNext)
    Vec[[i]]$Joint<-paste(Vec[[i]]$ActNow_sh,Vec[[i]]$ActNext_sh)
  }
  unlist_me<-lapply(Vec, function(x) x<-Mutual_Info_Grass(x))
  return_me<-rbindlist(unlist_me, idcol = "index")
  return(return_me)
}

#Example concatenated data, max distance for sequential comparison is 100 elements, seed is 43.
Permuted_Joint_MI(example_data, Length = 100, seed = 43) 
#(This may take a few seconds to run)





#Permuted_Joint_MI_Many - Run a fixed number of permutations of the joint distribution, calculate mean MI at each sequential distance, and return the values as a dataframe.
#x is the data, LenMI is the total length to compare elements over, and PermNo is the number of permutations to run.
Permuted_Joint_MI_Many<-function(x, Length, PermNo){
  storage<-data.frame(matrix(ncol = PermNo + 1, nrow = Length)) #Create dataframe to store iterative MI calculations on permuted data. 
  storage[,1] <- c(1:Length)
  for(i in 1:PermNo){
    est<-Permuted_Joint_MI(x, Length, i) #x is the data, LenMI is the distance to calculate MI over, i is the seed for this permutation.
    storage[,i+1] <-est$MInfo
  }
  mean_MI_Perm<-rowSums(storage[,-1])/PermNo
  storage<-cbind(storage, mean_MI_Perm)
  names(storage)[names(storage) == 'X1'] <- 'Index'
  return(storage)
}
#Note: This function may take several minutes to run (approx 5 minutes). 
Permuted_Joint_MI_Many(example_data, Length = 100, PermNo = 100) 
#Max distance for sequential comparrison is 100 elements.
#Calculations are run over 100 different permutations of the joint distribution.
#Permutations will be run using the random seeds 1-PermNo. In this case, it is 1-100.
#In the output, index refers to the comparative distance across sequence elements (e.g. 1,2,3...100 units).
#X1,X2,X3 etc refer to each permutation run by the function, e.g. X1 is the first permutation, X2 is the second, etc.
#A final column, mean_MI_perm, creates an average from all permutations, for each sequential distance.





#Permutation_MI_Adjust - Adjust for chance MI estimation using permuted data.
#Can take a concatenated sequence dataset and estimate MI and MI from permuted joint entropy.
#It then returns all in one dataframe, including adjusted MI values.
#Requires functions Iterative_MI_Grass, Permuted_Joint_MI, Mutual_Info_Grass, and JoinAct.
Permutation_MI_Adjust<-function(x,Length,PermNo){
  Estimated<-Iterative_MI_Grass(x,Length)
  Permuted<-Permuted_Joint_MI_Many(x,Length,PermNo)
  Estimated$Permuted_MInfo<-Permuted$mean_MI_Perm
  Estimated$Adjusted_MI<-Estimated$MInfo - Permuted$mean_MI_Perm
  return(Estimated)
}

#Note - this function can take a very long time to run. When running at 1000 permutations, running one individual's data set takes ~ 30-40 minutes
example_MIAdj<-Permutation_MI_Adjust(example_data, Length = 100, PermNo = 1000)  
#Max distance for sequential comparison is 100 elements.
#Calculations are run over 100 different permutations of the joint distribution.
#Permutations will be run using the random seeds 1-PermNo. In this case, it is 1-100.
example_MIAdj
#In the output table:
#index refers to inter-element distance over which MI was estimated.
#Xent refers to the marginal entropy for the anterior position of element comparisons.
#(e.g. for the sequence A,B,C,A,A, and comparisons are 2 elements apart, comparison units are (A,C), (B,A), (C,A) - Xent is calculated from A,B,C)
#Yent refers to the marginal entropy at the posterior position of element comparisons (e.g. from C,A,A from previous example)
#Jent is the joint entropy for comparisons (e.g. calculated from AC, BA, CA).
#N_J is total number of units for estimating Joint entropy.
#MInfo is the MI estimated before adjusting to permutation.
#Permuted_MInfo is mean MI from permuted distributions.
#Adjusted_MI is MI estimate adjusted to chance combinations of elements.





#### Estimating when MIAdj is at 0 ####

#For the rationale of this section, see sections 4.5 and 4.8 of main text.

#Add a column to the dataset which states if MIAdj is positive or negative.
#Run a logistic regression on P(MIAdj > 0) to identify if the relationship is significant, and if so,
#whether PP(MIAdj > 0) = 0.5 at any point within the 1-100 element range.

#Put positive and negatives into data as a column to be used in logistic regression.
is_negative<-function(x){
  x$Pos<-0
  for(i in 1:length(x$Adjusted_MI)){
    if(x$Adjusted_MI[i] < 0){
      x$Pos[i]<- 0
    }
    else{ 
      x$Pos[i] <- 1 
    }
  }
  return(x)
}

#Run on example data.
example_MIAdj<-is_negative(example_MIAdj)
example_logreg<-glm(Pos ~ index, data = example_MIAdj, family = "binomial")
summary(example_logreg) #A non-significant effect of inter-element distance on P(MIAdj > 0).





#### Logistic Regression Confidence Interval ####

#For those models which had a significant interaction between element distance and P(MIAdj > 0),
#a confidence interval for the model can be calculated with the following function.
#Here, we will use the example data, despite it not showing a significant relationship. This is purely to demonstrate the function in operation.

#x is a glm() object with the family "binomial"
confint_glm<-function(x){
  ilink<-family(x)$linkinv
  fitted<-as_tibble(predict(x, se.fit = TRUE))
  fitted$Upper_CI_LINK<-fitted$fit + 1.96*fitted$se.fit
  fitted$Lower_CI_LINK<-fitted$fit - 1.96*fitted$se.fit
  
  fitted<-mutate(fitted,
                 fit_prob = ilink(fitted$fit),
                 Upper_CI_Prob = ilink(fitted$Upper_CI_LINK),
                 Lower_CI_Prob = ilink(fitted$Lower_CI_LINK))
  return(fitted)
}

#Trial with the logistic regression run on the example sequence data.
confint_glm(example_logreg)





#Set negative values to estimations of 0, and then mark where zeros are ready for plotting (see below).

#Set negative values to estimations of 0.
ind<-which(example_MIAdj$Adjusted_MI <= 0)
for(i in 1:length(ind)){
  example_MIAdj$Adjusted_MI[ind[i]]<-0
}
example_MIAdj

#Mark where zeros are ready for plotting.
zero_num<-which(example_MIAdj$Adjusted_MI == 0) 
example_MIAdj$is_zero<-'NO'
for(i in 1:length(zero_num)){
  example_MIAdj$is_zero[zero_num[i]] <- 'YES'
}
example_MIAdj





#### Plotting Data ####

#Recreating plot of MIAdj with fitted models.

#Exponential - pred_Exp
nls.exp<-nlsLM(Adjusted_MI~(A)*exp(-index*B),data=example_MIAdj, start = list(A=1,B=1)) #Fit Exponential Model.
Pred_Exp<-predict(nls.exp)
pred_exp<-data.frame(example_MIAdj$index,Pred_Exp)
summary(nls.exp) #To view coefficients

#Power Law - pred_PL
nls.pl<-nlsLM(Adjusted_MI~A*(index**B),data=example_MIAdj, start = list(A=0.5,B=-0.1),control = nls.control(maxiter = 10000, minFactor = 0)) #Fit Power-Law Model.
Pred_PL<-predict(nls.pl)
pred_PL<-data.frame(example_MIAdj$index,Pred_PL)
summary(nls.pl)#To view coefficients

#Mixed - pred_MX
nls.mx<-nlsLM(Adjusted_MI~A*exp(-index*B) + C*(index**D),data=example_MIAdj, trace = T,start = list(A=0.1,B=0.1,C=0.1,D=0.1),control = nls.control(maxiter = 10000, minFactor = 0.0001)) #Fit Composite Model.
Pred_MX<-predict(nls.mx)
pred_MX<-data.frame(example_MIAdj$index,Pred_MX)
summary(nls.mx)#To view coefficients

#Plot with models.
plot_example<-ggplot(data=example_MIAdj, aes(x=log10(index), y=log10(Adjusted_MI))) +
  geom_point(data = example_MIAdj, alpha = 0.4, aes(shape = is_zero, size = is_zero),show.legend = FALSE) + #aes tells the factor to separate shapes by, and to have different sizes.
  scale_shape_manual(values = c(16,17)) + #Dictates which shapes get applied to which factor type.
  scale_size_manual(values = c(1, 2)) +
  xlab(expression(Log[10]~Distance)) +
  ylab(expression(Log[10]~MI[Adj]~(bits))) +
  ggtitle("Example MI Decay Data and Models") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-5,1) +
  geom_line(data=pred_PL, 
            aes(x = log10(example_MIAdj.index), y=log10(Pred_PL),colour = "red")) +
  geom_line(data=pred_exp, 
            aes(x = log10(example_MIAdj.index), y=log10(Pred_Exp),colour = "springgreen3")) +
  geom_line(data=pred_MX, 
            aes(x = log10(example_MIAdj.index), y=log10(Pred_MX),colour = "dodgerblue2")) +
  theme_light() +
  scale_colour_manual(name = 'Decay Model', 
                      values =c('springgreen3'='springgreen3','red'='red','dodgerblue2' = 'dodgerblue2'),guide = 'legend', labels = c('Exponential','Power-Law',"Composite")) +
  theme() 
plot_example 





#### Composite Model Transition Points ####

#For method and context, see sections 2.4 and 4.7 in main text.

#Function to get values for second differential of the composite model in logspace.
#A,B,C and D refer to the coefficients of the composite model.

second_dif_val<-function(A,B,C,D){
  Y <- expression(log10(A * exp(-(10**x) * B) + C * ((10**x)^D)))
  x<-log10(1:100)
  vals<-eval(D(D(Y, 'x'), 'x'))
  return(vals)
}

#Function to plot second differential, as a function of inter-element comparative distance, in logspace.
second_dif_plot<-function(df_values,name,ymin, ymax,xline,xlab,ylab){
  x<-log10(1:100)
  df_data<-data.frame(x,df_values)
  names(df_data)<-c("x", "df_values")
  ggplot(df_data, aes(x = x, y = df_values)) +
    ggtitle(name) +
    geom_line(alpha = 0.7) +
    ylab(ylab) +
    xlab(xlab) +
    ylim(ymin, ymax) +
    theme_light() +
    geom_hline(yintercept = 0, alpha = 0.5, color = "red") +
    geom_vline(xintercept = log10(xline), color = 'black', linetype = 'dashed', alpha = 0.5) 
}

A<-1.94702
B<-1.09597
C<-0.32573
D<- -0.97475
dif_2<-second_dif_val(A,B,C,D)
second_dif_plot(dif_2,"Example Second Differential Plot of Fitted Composite Model", -4, 4,3,expression(Log[10]~Distance),"f ''")





#### Function to Condense Sequences ####

#To condense repeated elements in the sequence into individual elements, use the following function (See section 2.3 in main text for context).

condense_sequence<-function(x){
  indexes<-vector()
  a = 1
  for(i in 2:length(x$Event)){
    if(x$Event[i] == x$Event[i-1]){
      indexes[a]<-i
      a = a + 1
    }
  }
  x_condensed<-x[-indexes,]
  return(x_condensed)
}

#Trial with original sequence data.
#Original length of corpus.
length(example_data$Behavior)
#Condensed Sequence Data.
condensed_data<-condense_sequence(example_data)
length(condensed_data$Behavior)

#This code can be used to check if any 'Event' elements in the sequence are neighbouring elements of the same type
for(i in 2:length(condensed_data$Event)){
  if(condensed_data$Event[i] == condensed_data$Event[i-1]){
    print("TRUE")}
  else{
    print("FALSE")
  }
}

#Condensed sequences can be passed to the same functions as before, such as Permutation_MI_Adjust for MI estimation and adjustment.