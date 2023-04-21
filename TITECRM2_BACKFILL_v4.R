##################################################
##    TITE-CRM for comparison (Real doses)  ######
##################################################
#back-filling 
#v2 can input number of backfilling cohorts
#v3 backfill once escalates above AND efficacy signal seen on at least one patient

library(rjags)
source("data_generation_functions-V2.R")
##pre-amble functions
#for two parameter logistic


expit<-function (x) {
  return(exp(x)/(1 + exp(x)))
}

logit<-function(p){
  return(log(p/(1-p)))
}
TITE_CRM_Post_Bayes_2parlog<-function(beta,weights_obs,DLT_obs,dose_obs,doses_vec,prior_mean,prior_var){
  
  if(!((length(weights_obs)==length(DLT_obs))&(length(DLT_obs)==length(dose_obs)))){
    stop("non-conforming lengths")
  }
  n.obs<-length(DLT_obs)
  Post<-c()
  for (i in 1:n.obs){
    
    Post[i]<-((weights_obs[i]*(doses_vec[dose_obs[i]]^exp(beta)))^DLT_obs[i])*
      (1-(weights_obs[i]*(doses_vec[dose_obs[i]]^exp(beta))))^(1-DLT_obs[i])
    
  }
  
  return(prod(Post)*(dnorm(beta,mean=prior_mean,sd=sqrt(prior_var))))
  
}


TITE_CRM_Post_Bayes_vector_Np<-function(beta,weights_obs,DLT_obs,dose_obs,doses_vec,prior_mean,prior_var){
  return(sapply(beta, function(x) TITE_CRM_Post_Bayes_Np(x,weights_obs,DLT_obs,dose_obs,doses_vec,prior_mean,prior_var)))
  
}
int_Bayes_Np<-function(beta,weights_obs,DLT_obs,dose_obs,doses_vec,prior_mean,prior_var){
  beta*TITE_CRM_Post_Bayes_vector_Np(beta,weights_obs,DLT_obs,dose_obs,doses_vec,prior_mean,prior_var)
}

dose_recom<-function(betahat,doses_vec,target){
  return(which.min(abs(target-doses_vec^exp(betahat))))
  
}


current_patient_data_frame<-function(current_time,patient.dataframe,follow_up){
  #input the whole data and convert to format for TITECRM
  #each row is a patient, with entries c("entry.time","dose.level","DLT","DLT.time") 
  num_patients<-max(patient.dataframe$patient_ID)
  
  patient_dataframe<-matrix(NA,nrow=num_patients,ncol=4)
  for (i in 1:num_patients){
    patient_data_ind<-patient.dataframe[patient.dataframe$patient_ID==i,]
    patient_DLT<-max(patient_data_ind$DLT)
    patient_entry.time<-patient_data_ind$entry_time[1]
    patient_dose.level<-patient_data_ind$dose_level[1]
    if(patient_DLT==1){
      patient_DLT.time<-min(which(patient_data_ind$DLT==1))+patient_entry.time
    }else{
      patient_DLT.time<-NA
    }
    patient_dataframe[i,]<-c(patient_entry.time,patient_dose.level,patient_DLT,patient_DLT.time)
  }
  
  patient_dataframe<-data.frame(patient_dataframe)
  names(patient_dataframe)<-c("entry.time","dose.level","DLT","DLT.time") 
  follow_up_time<-function(x) min(follow_up,x)
  patient_follow_up<- sapply(current_time-patient_dataframe$entry.time,follow_up_time)
  patient_weights<-patient_follow_up/follow_up
  
  
  #only DLT if we have seen it at the current time
  current.DLT<-patient_dataframe$DLT
  current.DLT[which(current.DLT==1)]<-patient_dataframe$DLT.time[which(current.DLT==1)]<=current_time
  
  #only DLT time observed id current DLT is true
  current.DLT.time<-patient_dataframe$DLT.time
  current.DLT.time[which(current.DLT==0)]<-NA
  
  patient_weights[which(current.DLT==1)]<-1
  entry.time<-patient_dataframe$entry.time
  dose.level<-patient_dataframe$dose.level
  return(data.frame(entry.time,dose.level,current.DLT,current.DLT.time,patient_weights,patient_follow_up))
  
}


patient_data_frame<-function(patient.dataframe){
  #input the whole data and convert to format for TITECRM
  #each row is a patient, with entries c("entry.time","dose.level","DLT","DLT.time") 
  num_patients<-max(patient.dataframe$patient_ID)
  #browser()
  patient_dataframe<-matrix(NA,nrow=num_patients,ncol=4)
  for (i in 1:num_patients){
    patient_data_ind<-patient.dataframe[patient.dataframe$patient_ID==i,]
    patient_DLT<-max(patient_data_ind$DLT)
    patient_entry.time<-patient_data_ind$entry_time[1]
    patient_dose.level<-patient_data_ind$dose_level[1]
    if(patient_DLT==1){
      patient_DLT.time<-min(which(patient_data_ind$DLT==1))+patient_entry.time
    }else{
      patient_DLT.time<-NA
    }
    patient_dataframe[i,]<-c(patient_entry.time,patient_dose.level,patient_DLT,patient_DLT.time)
  }
  
  patient_dataframe<-data.frame(patient_dataframe)
  names(patient_dataframe)<-c("entry.time","dose.level","DLT","DLT.time") 
  return(patient_dataframe)
  
}


gibbs_sampler.tite.crm.2parlog<-function(data_list,iter,model_string){
  model1.spec<-textConnection(model_string)
  jags <- jags.model(model1.spec,data =data_list,n.chains=1,n.adapt=iter,quiet=T)
  update(jags, iter,progress.bar="none")
  tt<-jags.samples(jags,c('alpha','beta'),iter,progress.bar="none")
  return(tt)
}

#create hard safety matrix

hard.safety.mat.function<-function(perc,co.size,max.cohorts){
  hard.safety.mat<-matrix(nrow=2,ncol=max.cohorts)
  hard.safety.mat[2,]<-seq(from=co.size,to=max.cohorts*co.size, by=co.size)
  for(i in 1:max.cohorts){
    numi<- hard.safety.mat[2,i]
    try_vec<-c(1:numi)
    probs_vec<-1-pbeta(0.3,1+try_vec,1+numi-try_vec)
    hard.safety.mat[1,i]<- min( which(100*probs_vec>perc))
  }                     
  return(hard.safety.mat)
}


##INPUT:
#seed=seed for reproducability
#tru_cyc1=true underlying probability of DLT in cycle 1
#co_size=cohort size
#ncohorts=maximum number of cohorts in study
#target=target toxicity (all cycles)
#doses=real doses
#prior_mean= prior mean
#prior_var= prior variance
#nTTP_array= 5x5x5 array of nTTP values for each of the 125 combinations of type and grade of toxicity
#ncycles=number of cycles
#cycle_dec=decrease in P(DLT) by cycle
#dose.skipping=enforce dose skipping rule? T=no dose skipping (default)
#sufficient.information==enforce stopping for sufficient information? Input number of assignments needed before stopping when the next assignment is the same. (and e.g. 99 for no sufficient info rule)
#hard.safety.rule=enforce hard safety rule based on Beta(1.1)? 
#options: 85= 2/3,3/6,4/9. 90=2/3,4/6,5/9, 95=3/3,4/6,5/9, 0=no hard safety enforced
#safety.stopping.low.unsafe= (T=stop when P(p1>0.3)>0.8 (cycle 1))
#safety.stopping.high.toosafe= (T= stop when P(pJ>0.3)>0.8 (cycle 1))
#kfold.skipping: Is a kfold dose skipping rule implemented?
#kfold: the k-fold rise thsat is allowable in dose skipping
#precision.stopping: Is a precision stopping rule implemented? Input number of assigned cohorts needed before precision stopping can be triggered
#initial.one.cycle: Is the initial period based on one cycle at a time? (F=wait until all cycles completed before next dose in initial period)
#backfill: should doses deemed safe be backfilled?
#backfill.num: How many cohorts to add as backfilling?
#efficacy_vec: Probability of efficacy for at least one patient out of cohort (default is rep(1,10), so if no input is specified, function will do "full backfilling")



##OUTPUT:
#list with elements corresponding to:
#dose.rec : dose recommendation
#num.pat: number of patients
#dose.ass : number of cohorts per dose (vector)
#stop.code: stopping reason (4=sufficient information,3= max patients reached)
#num.DLT: number of DLTS (total)
#grade.mat: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses. max grade observed by each patient.
#duration: total trial duration (until all recruited patients are fully observed)
#max_admissible: max admissible dose at the end of the trial (has hard safety eliminated any?)
#backfill: Should all doses below the highest explored be expanded? (2 cohorts)
#output.full.data: Should the full dataset be output?

TITECRM_2parlog_backfill<-function(seed,tru_cyc1,co_size,ncohorts ,target,doses,efficacy_vec=rep(1,10),p_a_mean=0,p_b_mean=0,
                                   p_a_prec=0.1,p_b_prec=0.1,nTTP_array,ncycles,cycle_dec,dose.skipping=T,
                                   sufficient.information=3,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                                   safety.stopping.high.toosafe=T,kfold.skipping=T,kfold=2,precision.stopping=3,
                                   initial.one.cycle=T,backfill=T,backfill.num=2,output.full.data=F){
  set.seed(seed)
  rands<-runif(co_size*ncohorts)
  randsE<-runif(ncohorts)
  stop_vec<-rep(0,7)
  ndoses<-length(doses)
  expand_vec<-rep(0, ndoses)
  dosevec<-c()
  patient_ID1<-1
  current.time<-0
  efficacy_obs<-rep(0,ndoses)
  
  #define all doses as admissable before any are dropped for safety
  max_admissable<-ndoses
  #hard safety
  if(hard.safety.rule>50){
    hard.safety<-T
    hard.safety.mat<-hard.safety.mat.function(perc=hard.safety.rule,co.size = co_size,max.cohorts = ncohorts)
  }else{
    hard.safety<-F
  }
  
  
  dose_rec<-NA
  nextdose<-1
  stop<-0
  
  
  initial<-1
  model.tite.crm.string <-"
model {




for(j in 1:ndoses){

logit(pi[j]) = alpha + beta*doses[j]
}


for(i in 1:length(yDLT)){
G[i]=weight[i]*pi[patdoses[i]]
yDLT[i] ~ dbinom(G[i],1)
}




log_beta ~ dnorm(mu_beta, tau_beta)
beta = exp(log_beta)
alpha ~ dnorm(mu_alpha, tau_alpha)



}
"


while(stop==0){
  
  
  #DATA GENERATION
  
  if(current.time==0){
    #first cohort
    esc.cohort<-1
    all.data<-multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose,entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec)
    current.time<-current.time+1
    patient_ID1<-max(all.data$patient_ID)+1
    current.data<-all.data[all.data$time_of<=current.time,]
    dosevec[current.time]<-nextdose
    efficacy_obs[1]<-as.numeric(randsE[1]<efficacy_vec[1])
  }else{
    #subsequent cohorts
    #escalation
    esc.cohort<-esc.cohort+1
    all.data<-rbind(all.data,multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose,entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec))
    efficacy_obs[c(nextdose:ndoses)]<- efficacy_obs[c(nextdose:ndoses)]+rep(as.numeric(randsE[esc.cohort]<efficacy_vec[nextdose]),length(c(nextdose:ndoses)))
    
     patient_ID1<-max(all.data$patient_ID)+1
    if(backfill==T){
      #expansion
      if(nextdose>1){ #if there are doses below the next dose
        below.doses<-c(1:(nextdose-1)) #which doses are below? 
        expand.below<- expand_vec[below.doses] #have the below doses been expanded?
        # browser()
        # which.expand<-which[expand.below==0] #which doses should now be expanded?
        for(expand.doses in below.doses){ #for each dose that needs expanding (not been exapnded before, and efficacy signal seen)
          if((expand.below[expand.doses]==0)&(efficacy_obs[expand.doses]>0)){
            #generate data for backfill.num cohorts on the dose
            all.data<-rbind(all.data,multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,expand.doses,entry_time=current.time,nTTP_array = nTTP_array,num_patients = backfill.num*co_size,rands=rands,cycle_dec=cycle_dec))
            #update patient ID for next assignment
            patient_ID1<-max(all.data$patient_ID)+1
            #update the expansion vector to say this dose has been expanded
            expand_vec[expand.doses]<-1
          }
        }
      }
      
    }
    
    current.time<-current.time+1
    
    current.data<-all.data[all.data$time_of<=current.time,]
    dosevec[current.time]<-nextdose
    
    
    
  }
  patient.data<-patient_data_frame(all.data)
  
  current.patient.data<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=ncycles)
  
  ##TITE CRM has an initial period where we keep escalating until we see a DLT
  
  if(initial==1){
    
    nextdose<-min(max(patient.data$dose.level)+1,length(doses))
    
    ##sufficient information 
    #  if(sufficient.information==T){
    npats_doses<-co_size*tabulate(dosevec,nbins = ndoses)
    num_ass<-tabulate(dosevec,nbins = ndoses) #number of assignments
    if(!is.na(nextdose)){
      if(num_ass[nextdose]>=sufficient.information){
        stop<-4
        stop_vec[4]<-1
        dose_rec<-nextdose
        break
      }
    }
    
    if(initial.one.cycle==T){
      if(sum(current.patient.data$current.DLT)>0){
        initial<-0
      }
      
    }else{
      
      for (cyc in 1:(ncycles-1)){
        current.patient.data<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=ncycles)
        current.time<-current.time+1
        
        if(sum(current.patient.data$current.DLT)>0){
          initial<-0
          current.time<-current.time-1
        }
        
      }
    }
    
  }
  if(initial==0){
    #posterior 
    n<-rep(0,ndoses)
    yDLT<-c()
    weight<-c()
    patdoses<-current.patient.data$dose.level
    
    #formatting for gibbs sampler
    for(obser in 1:nrow(current.patient.data)){
      yDLT[obser]<-current.patient.data$current.DLT[obser]
      weight[obser]<-current.patient.data$patient_weights[obser]
      n[current.patient.data$dose.level[obser]]<-n[current.patient.data$dose.level[obser]]+1
      
    }
    
    
    
    
    
    
    currentdata<-list(ndoses=ndoses,weight=weight,yDLT=yDLT,doses=doses,
                      patdoses=patdoses,
                      mu_alpha=p_a_mean,mu_beta=p_b_mean,
                      tau_alpha=p_a_prec,tau_beta=p_b_prec)
    
    
    gibbs_out<-gibbs_sampler.tite.crm.2parlog(data_list=currentdata,iter=10000,model_string=model.tite.crm.string)
    
    
    nextdose<-which.min(abs(target-expit(mean(gibbs_out$alpha)+doses*mean(gibbs_out$beta))))
    
    ##catches errors in the gibbs sampler caused by poor prior parameter combination
    if(length(nextdose)==0){
      stop<-1
      stop_vec[1]<-1
      dose_rec<-NA
      break
    }    
    
    #k fold dose increase on EXPERIMENTED DOSES
    if(kfold.skipping==T){
      
      prevdoseval<-doses[max(dosevec,na.rm=T)]
      nextdoseval<-doses[nextdose]
      
      if(nextdoseval>(kfold*prevdoseval)){
        maxdoseval<-kfold*prevdoseval
        nextdose<-max(which(doses<=maxdoseval))
      }
    }
    
    #stopping based on beta-binomial
    npats_doses<-co_size*tabulate(dosevec,nbins = ndoses)
    
    if(((safety.stopping.low.unsafe==T)&(npats_doses[1]>0))|((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0))){
      
      current.patient.data_cyc1<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=1)
      
      
      #n1<-rep(0,ndoses)
      yDLT1<-c()
      weight1<-c()
      patdoses1<-current.patient.data_cyc1$dose.level
      
      #formatting for gibbs sampler
      for(obser in 1:nrow(current.patient.data_cyc1)){
        yDLT1[obser]<-current.patient.data_cyc1$current.DLT[obser]
        weight1[obser]<-current.patient.data_cyc1$patient_weights[obser]
        #n1[current.patient.data_cyc1$dose.level[obser]]<-n1[current.patient.data_cyc1$dose.level[obser]]+1
        
      }
      
      
      
      
      
      
      currentdata1<-list(ndoses=ndoses,weight=weight1,yDLT=yDLT1,doses=doses,
                         patdoses=patdoses1,
                         mu_alpha=p_a_mean,mu_beta=p_b_mean,
                         tau_alpha=p_a_prec,tau_beta=p_b_prec)
      
      
      gibbs_out1<-gibbs_sampler.tite.crm.2parlog(data_list=currentdata1,iter=10000,model_string=model.tite.crm.string)
      
      if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
        posterior1<- expit((gibbs_out1$alpha)+doses[1]*(gibbs_out1$beta))
        
        posterior1<-posterior1[1:length(posterior1)]
        cyc1_0.3g<-mean(posterior1>0.3,na.rm=T)
        
        if(cyc1_0.3g>0.8){
          stop<-6
          stop_vec[6]<-1
          nextdose<-NA
          # break
        }
      }
      
      if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
        posteriorJ<-expit((gibbs_out1$alpha)+doses[ndoses]*(gibbs_out1$beta))
        posteriorJ<-posteriorJ[1:length(posteriorJ)]
        cycJ_0.3l<-mean(posteriorJ<0.3,na.rm=T)
        if(cycJ_0.3l>0.8){
          stop<-7
          stop_vec[7]<-1
          nextdose<-NA
          #  break
        }
      }
      
    }
    
    
    
    
    #number of DLTs per dose level
    nDLTs_doses<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)]))
    nDLTs_doses_c1<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)&(current.data$cycle==1)]))
    
    npats_doses<-sapply(c(1:ndoses), function(x) sum(current.patient.data$dose.level==x))
    num_ass<-tabulate(dosevec,nbins = ndoses) #number of assignments
    
    
    
    ##hard safety
    if(hard.safety==T){
      
      explored<-c(1:ndoses)[npats_doses>0]
      for(do in explored){
        
        if(nDLTs_doses_c1[do]>=hard.safety.mat[1, which(hard.safety.mat[2,]==npats_doses[do])]){
          
          max_admissable<-min(c(max_admissable,do-1))
          
        }
      }
      if(!is.na(nextdose)){
        #  browser()
        if(nextdose>max_admissable){
          nextdose<-max_admissable
        }
      }
      if(max_admissable==0){
        stop<-5
        stop_vec[5]<-1
        dose_rec<-NA
        nextdose<-NA
        # break
      }
      
      
    }
    
    
    
    
    ##sufficient information 
    #  if(sufficient.information==T){
    # npats_doses<-co_size*tabulate(dosevec,nbins = ndoses)
    # browser()
    if(!is.na(nextdose)){
      if(num_ass[nextdose]>=sufficient.information){
        stop<-4
        stop_vec[4]<-1
        dose_rec<-nextdose
        break
      }
    }
    
    ##precision stopping
    
    if(sum(num_ass)>=precision.stopping){
      MTD_vector<-(logit(target)-(gibbs_out$alpha))/(gibbs_out$beta)
      
      CV<-mad(MTD_vector, center = median(MTD_vector), constant = 1.4826, 
              na.rm = FALSE, low = FALSE, high = FALSE)/median(MTD_vector)
      
      if(CV<0.3){
        dose_rec<-nextdose
        stop<-2
        stop_vec[2]<-1
      }
    }
    
  }
  if(nrow(patient.data)==(ncohorts*co_size)){ #max patients reached
    
    current.patient.data<-current_patient_data_frame(current.time+ncycles,all.data,follow_up = ncycles)
    n<-rep(0,ndoses)
    yDLT<-c()
    weight<-c()
    patdoses<-current.patient.data$dose.level
    for(obser in 1:nrow(current.patient.data)){
      yDLT[obser]<-current.patient.data$current.DLT[obser]
      weight[obser]<-current.patient.data$patient_weights[obser]
      n[current.patient.data$dose.level[obser]]<-n[current.patient.data$dose.level[obser]]+1
      
    }
    
    
    
    
    
    
    currentdata<-list(ndoses=ndoses,weight=weight,yDLT=yDLT,doses=doses,
                      patdoses=patdoses,
                      mu_alpha=p_a_mean,mu_beta=p_b_mean,
                      tau_alpha=p_a_prec,tau_beta=p_b_prec)
    
    
    gibbs_out<-gibbs_sampler.tite.crm.2parlog(data_list=currentdata,iter=10000,model_string=model.tite.crm.string)
    
    
    
    
    dose_rec<-which.min(abs(target-expit(mean(gibbs_out$alpha)+doses*mean(gibbs_out$beta))))
    
    stop<-3
    stop_vec[3]<-1
  }
  
}


#follow up for all patients
current.data<-all.data
current.time<-max(current.data$time_of)


#grade matrix out
grade.matrix.out<-matrix(0,ncol=ndoses,nrow=5)
for (pat in 1:(max(current.data$patient))){
  pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$max_grade)
  pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
  grade.matrix.out[pat_max_grade+1,pat_dose]<- grade.matrix.out[pat_max_grade+1,pat_dose]+1
}

if(output.full.data==F){
  #output list
  output<-list(
    dose.rec=dose_rec ,#dose recommendation
    num.pat=max(current.data$patient) ,#: number of patients
    dose.ass=tabulate(dosevec,nbins = ndoses) ,# : number of cohorts per dose (vector)
    stop.code=stop_vec, #: stopping reason 
    num.DLT= sum(current.data$DLT[current.data$patient_ID>0]),#: number of DLTS (total)
    grade.mat=grade.matrix.out ,#: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses.
    duration= current.time,#: total trial duration (until all recruited patients are fully observed)
    max_admissable=max_admissable, # max admissable dose at the end of the trial (has hard safety eliminated any?)
    escalation=dosevec
  )
}else if(output.full.data==T){
  #output list
  output<-list(
    dose.rec=dose_rec ,#dose recommendation
    num.pat=max(current.data$patient) ,#: number of patients
    dose.ass=tabulate(dosevec,nbins = ndoses) ,# : number of cohorts per dose (vector)
    stop.code=stop_vec, #: stopping reason 
    num.DLT= sum(current.data$DLT[current.data$patient_ID>0]),#: number of DLTS (total)
    grade.mat=grade.matrix.out ,#: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses.
    duration= current.time,#: total trial duration (until all recruited patients are fully observed)
    max_admissable=max_admissable, # max admissable dose at the end of the trial (has hard safety eliminated any?)
    escalation=dosevec,
    full.trial.data=current.data
  )
}

return(output)

}

