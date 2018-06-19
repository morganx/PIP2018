library(ggplot2)
library(plyr)
library(vegan)
library(asbio)
library(reshape2)
##### functions for randomly sampling from each group by different independent variable(s)
# randomly select a number of samples (columns) as a subgroup
# make that a function to randomly select samples with two runs of test on a certain number of samples, use e. coli as a test for sampling depth of the selected number of samples per group
Random_samples_columns<-function(df,Samp_number){
  # if dataset is matrix, force it to be a dataframe
  subsamples<-df[,sample(ncol(df),Samp_number, replace = FALSE)]
  subsamples<-apply(subsamples, 2, function(x) x<-as.numeric(as.character(x)))
  # row.names(subsamples)<-row.names(df)
  # remove taxon that has 0 in all selected samples
  # subsamples<-subsamples[rowSums(subsamples)!=0,]
  
  # transpose data so each sample will be a row instead of a column (this is required for further PERMANOVA analyses)
  subsamples<-as.data.frame(t(subsamples), row.names<-colnames(subsamples) )
  return(subsamples)
}

# make a function to randomly select a number of samples from different groups, and calculate their PERMANOVA (permutation = 999). Note this function will need the function (Random_samples_columns) to run
# inde_variable needs a string (eg. inde_variable = "time")
### Consider to make the independent variable to multiple later on ###
PERMANOVA_from_different_groups<-function(df, Samp_number, metadf, inde_variable) {
  # separate samples based on the selected independent_variable (e.g. time)
  # defind groups by factor levels in the independent variable.
  lev<-levels(metadf[,inde_variable]) # record the factor levels
  #### Note if data has been pre-filtered and may not contain all factor levels, use the following command instead
  # lev<-unique(metadf[,inde_variable])
  x<-split(metadf, metadf[,inde_variable]) # this will separate the dataframe into a list of dataframes, one dataframe each factor level (same order as lev).
  
  # create a emplty list for saving the subsamples
  temp_list<-list()
  for (subgroup in 1:length(lev)) { # for each subgroup
    # Randomly select some number of samples  (Samp_number as the input of the function) from the taxonomy dataframe (df as the input of the function) as a subdataset, and assign each subdataset into a list one by one. So the list will have a few number of elements, that equals to the number of factor levels of the selected independent variable.
    temp_list[[subgroup]]<-as.data.frame(Random_samples_columns(df[,row.names(x[[subgroup]])],Samp_number))
  }
  # calculate the Bray-Curtis dissimilarities (default for Vegan) between subsamples groups
  dissim<-vegdist(do.call(rbind,temp_list))
  # set up the metadata as factor for later PERMANOVA test
  fac<-as.factor(metadf[row.names(do.call(rbind,temp_list)),inde_variable])
  # Permutational multivariate analysis of variance using Bray-Curtis distance matrices, with 999 permutations
  result<-adonis(dissim~fac)
  return(result)
  # print(lev)
  # print(result)
}

##### make functions to calculate PERMANOVA at different sampling depth

# calculate power for a number of runs with a certain sampling depth 
P_calculator<-function(df, Samp_number, metadf, inde_variable, Power_runs) {
  test<-replicate(Power_runs, {
    temp<-PERMANOVA_from_different_groups(df, Samp_number, metadf, inde_variable)
    return(temp$aov.tab[1,6])
  })
}

# now make a function that can hold a list numbers of sampling depth. e.g. sampling_depth<-seq(5,20,5) )
P_different_depth<-function(df, sampling_depth, metadf, inde_variable, Power_runs) {
  result<-list()
  for (depth in sampling_depth) {
    pos<-match(depth, sampling_depth)
    temp<-as.data.frame(P_calculator(df, depth, metadf, inde_variable, Power_runs))
    colnames(temp)<-"PERMANOVA_p_value"
    temp$sampling_depth<-depth
    result[[pos]]<-temp
  }
  return(result)
}

##### Power defined as the fraction of subsamples that are significant tested by PERMANOVA. Then plot the figure with power as y-axis and sampling depth as x-axis based on the functions we made before 

Power_calculation<-function(df, sampling_depth, metadf, inde_variable, Power_runs) {
  result<-P_different_depth(df, sampling_depth, metadf, inde_variable, Power_runs)
  combind_result<-do.call(rbind,result)
  # made a new binary variable "significance" for separating the P-values:  1 for p<=0.05, else 0 for the convinience of calculating the power
  combind_result$significance<-ifelse(combind_result$PERMANOVA_p_value<=0.05, 1, 0)
  #combind_result$sampling_depth<-as.factor(combind_result$sampling_depth)
  Pow<-ddply(combind_result,"sampling_depth", summarise, 
             Pw=sum(significance)/length(significance))
  return(Pow)
  
}

##### multiple rounds of power tests for mean and confidence interval
Power_multiple_rounds<-function(df, sampling_depth, metadf, inde_variable, Power_runs, rounds){
  test<-replicate(rounds, {Power_test<-Power_calculation(df, sampling_depth, metadf, inde_variable, Power_runs)
  return(Power_test)})
  result<-do.call(cbind,t(test))
  row.names(result)<-result[,1]
  result<-result[,-seq(1,rounds,1)]
  return(result)
}

# calculate CI for figure
Power_CI_figure<-function(df, sampling_depth, metadf, inde_variable, Power_runs, rounds){
  t<-Power_multiple_rounds(df, sampling_depth, metadf, inde_variable, Power_runs, rounds)
  t1<-as.data.frame(apply(t,1, mean))
  t2<-as.data.frame(apply(t,1, stan.error))
  t2<-cbind(t1,t2)
  colnames(t2)<-c("mean","standard_error")
  t2$CI_low<-t2$mean-qt(0.975, df=rounds-1)*t2$standard_error
  t2$CI_high<-t2$mean+qt(0.975, df=rounds-1)*t2$standard_error
  t2$Sampling_depth<-as.numeric(row.names(t2))
  t2$variable<-inde_variable
  return(t2)
}

##===============Monte-carlo on Alpha diversity=============##
# use some of the functions made for beta diversity and create a few for alpha diversity
Shannon_from_different_groups<-function(df, Samp_number, metadf, inde_variable) {
  # separate samples based on the selected independent_variable (e.g. time)
  # defind groups by factor levels in the independent variable.
  lev<-levels(metadf[,inde_variable]) # record the factor levels
  #### Note if data has been pre-filtered and may not contain all factor levels, use the following command instead
  # lev<-unique(metadf[,inde_variable])
  x<-split(metadf, metadf[,inde_variable]) # this will separate the dataframe into a list of dataframes, one dataframe each factor level (same order as lev).
  
  # create a emplty list for saving the subsamples
  temp_list<-list()
  for (subgroup in 1:length(lev)) { # for each subgroup
    # Randomly select some number of samples  (Samp_number as the input of the function) from the taxonomy dataframe (df as the input of the function) as a subdataset, and assign each subdataset into a list one by one. So the list will have a few number of elements, that equals to the number of factor levels of the selected independent variable.
    temp_list[[subgroup]]<-as.data.frame(Random_samples_columns(df[,row.names(x[[subgroup]])],Samp_number))
  }
  # calculate the Bray-Curtis dissimilarities (default for Vegan) between subsamples groups
  shannon<-as.data.frame(diversity(do.call(rbind,temp_list), index = "shannon"))
  colnames(shannon)<-"Shannon_diversity"
  # set up the metadata as factor for later test
  shannon$fac<-as.factor(metadf[row.names(shannon),inde_variable])

  # kruskal test to get the p value
  result<-kruskal.test(Shannon_diversity~fac, data = shannon)
  return(result)
  # print(lev)
  # print(result)
}

# calculate power for a number of runs with a certain sampling depth 
Shannon_P_calculator<-function(df, Samp_number, metadf, inde_variable, Power_runs) {
  test<-replicate(Power_runs, {
    temp<-Shannon_from_different_groups(df, Samp_number, metadf, inde_variable)
    return(temp[[3]]) # when use [n], gives you the whole nth item in the list, if you [[n]], gives the value(s) of the n th item in the list.
  })
}

# now make a function that can hold a list numbers of sampling depth. e.g. sampling_depth<-seq(5,20,5) )
Shannon_P_different_depth<-function(df, sampling_depth, metadf, inde_variable, Power_runs) {
  result<-list()
  for (depth in sampling_depth) {
    pos<-match(depth, sampling_depth)
    temp<-as.data.frame(Shannon_P_calculator(df, depth, metadf, inde_variable, Power_runs))
    colnames(temp)<-"Kruskal_p_value_for_shannon_diversity"
    temp$sampling_depth<-depth
    result[[pos]]<-temp
  }
  return(result)
}

##### Power defined as the fraction of subsamples that are significant tested by kruskal test. Then plot the figure with power as y-axis and sampling depth as x-axis based on the functions we made before 

Shannon_Power_calculation<-function(df, sampling_depth, metadf, inde_variable, Power_runs) {
  result<-Shannon_P_different_depth(df, sampling_depth, metadf, inde_variable, Power_runs)
  combind_result<-do.call(rbind,result)
  # made a new binary variable "significance" for separating the P-values:  1 for p<=0.05, else 0 for the convinience of calculating the power
  combind_result$significance<-ifelse(combind_result$Kruskal_p_value_for_shannon_diversity<=0.05, 1, 0)
  #combind_result$sampling_depth<-as.factor(combind_result$sampling_depth)
  Pow<-ddply(combind_result,"sampling_depth", summarise, 
             Pw=sum(significance)/length(significance))
  return(Pow)
  
}

##### multiple rounds of power tests for mean and confidence interval
Shannon_Power_multiple_rounds<-function(df, sampling_depth, metadf, inde_variable, Power_runs, rounds){
  test<-replicate(rounds, {Power_test<-Shannon_Power_calculation(df, sampling_depth, metadf, inde_variable, Power_runs)
  return(Power_test)})
  result<-do.call(cbind,t(test))
  row.names(result)<-result[,1]
  result<-result[,-seq(1,rounds,1)]
  return(result)
}

# calculate CI for figure
Shannon_Power_CI_figure<-function(df, sampling_depth, metadf, inde_variable, Power_runs, rounds){
  t<-Shannon_Power_multiple_rounds(df, sampling_depth, metadf, inde_variable, Power_runs, rounds)
  t1<-as.data.frame(apply(t,1, mean))
  t2<-as.data.frame(apply(t,1, stan.error))
  t2<-cbind(t1,t2)
  colnames(t2)<-c("mean","standard_error")
  t2$CI_low<-t2$mean-qt(0.975, df=rounds-1)*t2$standard_error
  t2$CI_high<-t2$mean+qt(0.975, df=rounds-1)*t2$standard_error
  t2$Sampling_depth<-as.numeric(row.names(t2))
  t2$variable<-inde_variable
  return(t2)
}
