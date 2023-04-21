##example of simulations


#packages
library(doParallel)
registerDoParallel(cores=36)

#source code
source("TITECRM2_BACKFILL_v4.R")



##define scenarios
#different dose levels MTD (linear around MTD)
s1<-c(0.30, 0.40, 0.50 ,0.60, 0.70, 0.80)
s2<-c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70)
s3<-c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60)
s4<-c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
s5<-c(0.05, 0.10, 0.15, 0.20, 0.30 ,0.40)
s6<-c(0.02, 0.05, 0.10, 0.15, 0.20, 0.30)

#non-linear around MTD
s7<-c(0.15, 0.20, 0.25, 0.30, 0.45, 0.60)
s8<-c(0.05, 0.15, 0.30, 0.35, 0.40, 0.45)

#no dose has exactly target
s10<-c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55)
s11<-c(0.15, 0.20, 0.35, 0.40, 0.45, 0.50)
s12<-c(0.05, 0.10, 0.15, 0.20, 0.25, 0.40)

#all unsafe (but closer to MTD than D)
s9<-c(0.40, 0.45, 0.50, 0.55, 0.60, 0.65)

#defined scenarios
s13<-c(0.06, 0.07, 0.08, 0.09, 0.11, 0.12) #A
s14<-c(0.10, 0.14, 0.21, 0.30, 0.46, 0.58) #B
s15<-c(0.16, 0.30, 0.50, 0.70, 0.89, 0.95) #C
s16<-c(0.55, 0.91, 0.99, 1.00, 1.00, 1.00) #D
s17<-c(0.05, 0.05, 0.05, 0.80, 0.80, 0.80) #E

scen_list<-lapply(c(1:17),function(x) get(paste(c("s",x),collapse="")))



##define fixed inputs
co_size<-3
ncohorts<-18
doses<-c(1.5,2.5,3.5,4.5,6,7)
target1<-0.3
target3<-0.391
ncycles<-3
nsims<-5000
cycle_dec<-1/3

efficacy_vec<-c(0,0.15,0.3,0.45,0.6,0.75)


#execute simulations


for(sc in 1:length(scen_list)){
  
  
  cal_time <- system.time({
    
    
    ##partial backfill & 1 cycle
    ##define parameter values 
    p_a_mean<-log(0.5)
    p_b_mean<-0
    p_a_prec<-1/4
    p_b_prec<-1
    assign(paste(c("BF1e_simulation_s",sc),collapse = ""), foreach(i=1:nsims, combine = list) %dopar% {
      ##function
      
      TITECRM_2parlog_backfill(seed=i,tru_cyc1=scen_list[[sc]],co_size=co_size,ncohorts=ncohorts ,target=target1,doses=doses,
                               efficacy_vec = efficacy_vec, p_a_mean=p_a_mean,p_b_mean=p_b_mean,
                               p_a_prec=p_a_prec,p_b_prec=p_b_prec,nTTP_array=nTTP_array,ncycles=1,cycle_dec=cycle_dec,dose.skipping=T,
                               sufficient.information=3,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                               safety.stopping.high.toosafe=T,kfold.skipping=T,kfold=2,precision.stopping=3,initial.one.cycle=T,backfill=T)
      
      
    })
    save.image("backfill_sims_e.RData")
    
    
    
    ##partial backfill & 3 cycle
    
    ##define parameter values 
    p_a_mean<-log(0.5)
    p_b_mean<-0
    p_a_prec<-1/4
    p_b_prec<-1
    assign(paste(c("BF3e_simulation_s",sc),collapse = ""), foreach(i=1:nsims, combine = list) %dopar% {
      ##function
      
      TITECRM_2parlog_backfill(seed=i,tru_cyc1=scen_list[[sc]],co_size=co_size,ncohorts=ncohorts ,target=target3,doses=doses,
                               efficacy_vec = efficacy_vec, p_a_mean=p_a_mean,p_b_mean=p_b_mean,
                               p_a_prec=p_a_prec,p_b_prec=p_b_prec,nTTP_array,ncycles=3,cycle_dec=cycle_dec,dose.skipping=T,
                               sufficient.information=3,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                               safety.stopping.high.toosafe=T,kfold.skipping=T,kfold=2,precision.stopping=3,initial.one.cycle=T,backfill=T)
      
      
    })
    save.image("backfill_sims_e.RData")
    
    
  })
  
  
}

