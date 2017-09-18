##' This function simulates genetic drift of variants at a neutral polymorphism 
##' in populations under selection based on different initial variant frequencies. 
##' Implemented in research article: Chen S, Kaeppler SM, Vogel KP, Casler MD (2016) Selection Signatures in Four Lignin Genes from Switchgrass Populations Divergently Selected for In Vitro Dry Matter Digestibility. PLoS ONE 11(11): e0167005. doi:10.1371/journal.pone.0167005

##' Input: ifreq is the initial variant frequency
##'        n.allel is the ploid level
##'        cycle and cyclediv are the number of selection cycles for divergent directions, 
##'  example: cycle=3, 3 cycles of selection for high level of lignin content; 
##'  cyclediv=3, 3 cycles of selection of low level of lignin content
##'        Csize is a vector of populations sizes at each selection cycle starting from base population at cycle 0
##'        selecsize is a vector of the number of individuals selected at each cycle starting from cycle 1. 
##' Note: assuming that the divergent directions select the same number of individuals and keep the same population size at each cycle
##'        sim.rep is the number of simulations
##' Note: The following parameters account for the errors occuring at genotype sampling and sequence sampling. 
##'        samplesize and samplesizediv are vectors representing the number of individuals were sampled for estimating the variant frequencies
##'        seqsize and seqsizediv are vectors representing the number of reads obtained at each polymorphism for variant frequency estimation 

##' Output: all simulated variant frequencies of all selection cycles for divergent directions

simpoly<-function(ifreq,n.allel,cycle,cyclediv,Csize,selecsize,seqsize,seqsizediv,samplesize,samplesizediv,sim.rep){
  SimFreq<-matrix(,nrow=sim.rep,ncol=((cycle+1)*length(ifreq)))
  SimFreqdiv<-matrix(,nrow=sim.rep,ncol=((cyclediv+1)*length(ifreq)))
  for (infq in 1: length(ifreq)){ #ifreq was estimated from the C0 population sequences 
    #print(infq)
    alleltype=c("a","b")
    # Generations with drift
   
    pop<-matrix(,nrow=max(Csize),ncol=n.allel)
    samplefreq<-numeric((cycle+1))
    samplefreqdiv<-numeric((cyclediv+1))
    simpop <-character(Csize[1]*n.allel)
    allelprob<-c(ifreq[infq],1-ifreq[infq])
    
    
    for (simid in 1:sim.rep){
      #print(simid)
      pop<-matrix(,nrow=max(Csize),ncol=n.allel)
      popdiv<-matrix(,nrow=max(Csize),ncol=n.allel)
      # C0 population with ESTMATED initial allele frequency ifreq[]
      simpop=sample(alleltype,Csize[1]*n.allel,prob=allelprob,replace=TRUE)
      # Assign the alleles to individuals to the C0 population with n.allel number of variables and # of individuals as nrow
      for (id in 1:Csize[1]){  
        Indtemp=simpop[(1+(id-1)*n.allel):(id*n.allel)]
        pop[id,]<-Indtemp
      }
      rownames(pop)=NULL
      popdiv<-pop



# The allele frequency of the C0 population is NOT!!! normally distributed with mean of p and SE of sqrt(pq/n)
#samplefreq[1]<-rnorm(1,ifreq[infq],sqrt(ifreq[infq]*(1-ifreq[infq])/seqsize[1])) #seqsize[] is the sequencing number for each population 

samplefreq[1]<-prop.table(table(sample(pop,seqsize[1],replace=TRUE)))["a"]
samplefreqdiv[1]<-samplefreq[1]
      
      for (cyl in 1: cycle){ #one direction of selection
        
         # this is for debugging the sample() vector size NA error
        a<-cyl+1
        b<-seqsize[a]
        
        selecID<-sample(1:Csize[cyl],selecsize[cyl],replace=FALSE) # Build selection nursery
        D<-samplesize[cyl]+Csize[cyl+1] 
        E<-Csize[cyl+1]
        G<-Csize[cyl+1]+1
        temp<-matrix(,nrow=max(Csize),ncol=n.allel) # The matrix temp[,] to store progenies of the selection nursery
        for (n in 1:D){
          parentsID<-sample(selecID,2,replace=FALSE) # Randomly choose two individuals as parents
          parents<-pop[parentsID,]
          gamete1<-c(sample(parents[1,1:4],0.25*n.allel,replace=FALSE),sample(parents[1,5:8],0.25*n.allel,replace=FALSE))
          gamete2<-c(sample(parents[2,1:4],0.25*n.allel,replace=FALSE),sample(parents[2,5:8],0.25*n.allel,replace=FALSE))
          progeny<-c(gamete1,gamete2) 
          temp[n,]<-progeny
        } #the end of polycross
        pop<-matrix(,nrow=max(Csize),ncol=n.allel)
        pop<-temp[1:E,]
        pop<-na.omit(pop)
        Sample<-matrix(,nrow=samplesize[a],ncol=n.allel)  
        Sample<-temp[G:D,] # Individual sample error. Matrix Sample[,] contains the evaluated individuals
        
        Sample<-sample(Sample,b,replace=TRUE)
        samplefreq[cyl+1]<-prop.table(table(Sample))["a"] #sequencing sample error
          
        } # the end of cycles  
      SimFreq[simid,(1+(infq-1)*(cycle+1)):(infq*(cycle+1))]<-samplefreq[1:(cycle+1)]
       
      for (cyl in 1: cyclediv){ #divergent selection
        
        selecID<-sample(1:Csize[cyl],selecsize[cyl],replace=FALSE)
        D<-samplesizediv+Csize[cyl+1]
        E<-Csize[cyl+1]
        G<-Csize[cyl+1]+1
        temp<-matrix(,nrow=max(Csize),ncol=n.allel)
        for (n in 1:D){
          parentsID<-sample(selecID,2,replace=FALSE)
          parents<-popdiv[parentsID,]
          gamete1<-c(sample(parents[1,1:4],0.25*n.allel,replace=FALSE),sample(parents[1,5:8],0.25*n.allel,replace=FALSE))
          gamete2<-c(sample(parents[2,1:4],0.25*n.allel,replace=FALSE),sample(parents[2,5:8],0.25*n.allel,replace=FALSE))
          
          progeny<-c(gamete1,gamete2) 
          temp[n,]<-progeny
        } # the end of polycross
        popdiv<-matrix(,nrow=max(Csize),ncol=n.allel)
        popdiv<-temp[1:E,]
        popdiv<-na.omit(pop)
        Sample<-temp[G:D,] # Individual sample error for divergent selection

        Sample <- sample(as.vector(Sample),seqsizediv,replace=TRUE)
        samplefreqdiv[cyl+1]<-prop.table(table(Sample))["a"] # Sequencing sample error for divergent selection
      } # the end of divergent cycles
      SimFreqdiv[simid,(1+(infq-1)*(cyclediv+1)):(infq*(cyclediv+1))]<-samplefreqdiv[1:(cyclediv+1)]
      

    } # the end of simrep
    
    SimFreq[is.na(SimFreq)]<-0
    SimFreqdiv[is.na(SimFreqdiv)]<-0
  } # the end of initial freqencies
  result<-list(SimFreq,SimFreqdiv)
  return(result)
} # the end of the function
