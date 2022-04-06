####Custom R Functions to Accompany MI Decay Paper####

#### MI Estimation ####

#Data inputs take the form of a dataframe with the following columns:
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
#'Event' concatenated description of both behaviour, e.g. grasp, and modifier object, e.g. Hammer

#An example datastructure can be generated with the lines of code below:
X<-c(1,2,3,4,5,6)
Time<-c(657.1,661.2,665.3,667.4,667.8,668.2)
Subject<-c("Fanle","Fanle","Fanle","Fanle","Fanle","Fanle")
Behavior<-c("grasp","place","grasp","strikeonehand","drop","grasp")
Modifier.1<-c("1NUT","1NUT","2HAMMER","2HAMMER","2HAMMER","4KERNEL")
Hammer<-c("34","34","34","34","34","34")
Anvil<-c("32","32","32","32","32","32")
Complete<-c("N","N","N","N","N","N")
Nut_Crack<-c("N","N","N","N","N","N")
Nut_Number<-c("NA","NA","NA","NA","NA","NA")
Event<-c("grasp 1NUT","place 1NUT", "grasp 2HAMMER", "strikeonehand 2HAMMER", "drop 2HAMMER", "grasp 4KERNEL")

Example_Data<-data.frame(X, Time, Subject, Behavior, Modifier.1, Hammer, Anvil, Complete, Nut_Crack, Nut_Number, Event)
Example_Data

#Functions Included:

#JoinAct - Create a dataframe with the combined actions in (Action now and action 1, 2, 3... elements along). Comparative distance is set at the start of the function.

#Mutual_Info_Grass - Grassberger Entropy Estimation and MI Calculation 

#Iterative_MI_Grass- Iterative Grassberger Estimation and MI Calculation - Requires Grass Entropy & MI Function 

#Permuted_Joint_MI - Create dataframe with Permuted Joint Distributions

#permuted_Joint_MI_Many - Run a fixed number of permutations of the joint distribution, calculate average, and return

#Permutation_MI_Adjust: Estimate MI and adjust for chance MI estimation using permuted data

#TO RUN MI ESTIMATION, HAVE ALL FUNCTIONS LOADED IN R (as some are dependent on one another) AND RUN Permutation_MI_Adjust 
#FOR FULL MI ESTIMATION AND PERMUTATION ADJUSTMENT.



#JoinAct - Create a dataframe with the combined actions in (Action now and action 1, 2, 3... elements along). Comparative distance is set at the start of the function.
#x is dataframe containing actions. y is the distance apart you want to count the actions
JoinAct<-function(x,y){
  e1 <- vector()
  e2 <- vector()
  len<-length(x$Event)
  for(i in (y+1):len){
    e1[i]<-x$Event[i-y]
    e2[i]<-x$Event[i]
  }
  
  j <- data.frame(e1,e2)
  j$Joint<-paste(j$e1,j$e2) #Paste the two columns to get a column with both elements in
  names(j)[names(j) == "e1"] <- "ActNow" 
  names(j)[names(j) == "e2"] <- "ActNext"
  j <- na.omit(j) #Remove NA rows
  return(j)
}



#Mutual_Info_Grass - Insert a dataframe that has probabilities of Now, Next and Joint.
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



#Iterative_MI_Grass - Insert an unjoint dataframe, and the maximum distance you want to compare over.
Iterative_MI_Grass<-function(x,Length){  #Length is the max distance you want to make comparrisons over
  empty_list <- vector(mode = "list", length = Length)
  for(i in 1:Length){
    empty_list[[i]]<-JoinAct(x, i)
  }
  ret<-empty_list
  l1<-lapply(ret, function(x) x<-Mutual_Info_Grass(x))
  d.info<-rbindlist(l1, idcol = "index")
  return(d.info)
} 



#Permuted_Joint_MI - Create dataframe with Permuted Joint Distributions

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



#permuted_Joint_MI_Many - Run a fixed number of permutations of the joint distribution, calculate average, and return
#x is the data, LenMI is the total length to compare elements over, and PermNo is the number of permutations to run.
permuted_Joint_MI_Many<-function(x, LenMI, PermNo){
  storage<-data.frame(matrix(ncol = PermNo + 1, nrow = LenMI)) #Create dataframe to store iterative MI calculations on permuted data. 
  storage[,1] <- c(1:LenMI)
  for(i in 1:PermNo){
    est<-Permuted_Joint_MI(x, LenMI, i) #x is the data, LenMI is the distance to calculate MI over, i is the seed for this permutation.
    storage[,i+1] <-est$MInfo
  }
  mean_MI_Perm<-rowSums(storage[,-1])/PermNo
  storage<-cbind(storage, mean_MI_Perm)
  names(storage)[names(storage) == 'X1'] <- 'Index'
  return(storage)
}

#Permutation_MI_Adjust - Adjust for chance MI estimation using permuted data
#Can take a dataset that has undegone obact combine and estimate MI and MI from permuted joint entropy
#Then returns all in one dataframe
#Requires functions Iterative_MI_Grass, Permuted_Joint_MI, Mutual_Info_Grass, and JoinAct
Permutation_MI_Adjust<-function(x,Len,PermNo){
  Estimated<-Iterative_MI_Grass(x,Len)
  Permuted<-permuted_Joint_MI_Many(x,Len,PermNo)
  Estimated$Permuted_MInfo<-Permuted$mean_MI_Perm
  Estimated$Adjusted_MI<-Estimated$MInfo - Permuted$mean_MI_Perm
  return(Estimated)
}





#### Logistic Regression Confidence Interval ####

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



#### Composite Model Transition Points ####

#Function to get values for second differential of the composite function.
#A,B,C and D refer to the coefficients of the model model, estimated as outlined in the methods section.
second_dif_val<-function(A,B,C,D){
  Y <- expression(log10(A * exp(-(10**x) * B) + C * ((10**x)^D)))
  x<-log10(1:100)
  vals<-eval(D(D(Y, 'x'), 'x'))
  return(vals)
}
