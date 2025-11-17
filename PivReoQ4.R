PivReoQ4<-function(mcombinedchains,
                   piv){
  library(combinat)
  library(label.switching)
  
  if (!is.atomic(piv) || length(piv) != 1) {
    stop("Input piv is not a scalar. Execution stopped.")
  }
  if (!is.matrix(mcombinedchains) ||
      dim(mcombinedchains)[1] != 40000 ||
      dim(mcombinedchains)[2] != 177) {
    stop("Input mcombinedchains is not a matrix or has
         incorrect dimensions (it should be 40000 X 177).
         Execution stopped.")
  }
  
  perm=permn(1:3)# all permutation of 3 labels
  # extract distribution of means, vars and weights
  mcmc.pars=array(NA,c(dim(mcombinedchains)[1],3,3))#
  for(i in 1:dim(mcombinedchains)[1]){
    mcmc.pars[i,,1]=mcombinedchains[i,c(1:3)]# means var 1
    mcmc.pars[i,,2]=mcombinedchains[i,c(13:15)] #weights
    mcmc.pars[i,,3]=mcombinedchains[i,c(16:18)]# variances var 1
  }
  #
  # extract pivot values
  pivo=array(NA,c(3,3))
  pivo[,1]=mcombinedchains[piv,c(1:3)]
  pivo[,2]=mcombinedchains[piv,c(13:15)]
  pivo[,3]=mcombinedchains[piv,c(16:18)]
  #
  #actual algorithm
  run=pra(mcmc=mcmc.pars, pivot=pivo)
  #
  #
  Ord_mcombinedchainsNEW=mcombinedchains
  for(i in 1:dim(mcombinedchains)[1]){
    #
    # reorder ith iteration MU[1:3,1] by the selected permutation
    reord1=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,1:3]))
    reord1=rbind(reord1,reord1[2,c(as.vector(reord1[1,1:3]))])
    #
    # reorder ith iteration MU[1:3,2] by the selected permutation
    reord2=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,4:6]))
    reord2=rbind(reord2,reord2[2,c(as.vector(reord2[1,1:3]))])
    # 
    # reorder ith iteration MU[1:3,3] by the selected permutation
    reord3=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,7:9]))
    reord3=rbind(reord3,reord3[2,c(as.vector(reord3[1,1:3]))])
    #
    # reorder ith iteration MU[1:3,4] by the selected permutation
    reord4=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,10:12]))
    reord4=rbind(reord4,reord4[2,c(as.vector(reord4[1,1:3]))])
    # 
    # reorder ith iteration PSIs by selected permutation
    reord5=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,13:15]))
    reord5=rbind(reord5,reord5[2,c(as.vector(reord5[1,1:3]))])
    #
    # reorder ith iteration sigma[1:3,1] by the selected permutation
    reord6=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,16:18]))
    reord6=rbind(reord6,reord6[2,c(as.vector(reord6[1,1:3]))])
    #
    # reorder ith iteration sigma[1:3,2] by the selected permutation
    reord7=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,19:21]))
    reord7=rbind(reord7,reord7[2,c(as.vector(reord7[1,1:3]))])
    # 
    # reorder ith iteration sigma[1:3,3] by the selected permutation
    reord8=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,22:24]))
    reord8=rbind(reord8,reord8[2,c(as.vector(reord8[1,1:3]))])
    #
    # reorder ith iteration sigma[1:3,4] by the selected permutation
    reord9=rbind((as.vector(run$permutations[i,])),(mcombinedchains[i,25:27]))
    reord9=rbind(reord9,reord9[2,c(as.vector(reord9[1,1:3]))])
    #
    # order ith iteration Zs by selected permutation.
    reordZ=cbind((as.vector(run$permutations[i,])),1:3)
    for(j in 28:150){
      Ord_mcombinedchainsNEW[i,j]=reordZ[mcombinedchains[i,j],1]
    }
    #
    # save the new ordered 
    Ord_mcombinedchainsNEW[i,1:3]=reord1[3,1:3]
    Ord_mcombinedchainsNEW[i,4:6]=reord2[3,1:3]
    Ord_mcombinedchainsNEW[i,7:9]=reord3[3,1:3]
    Ord_mcombinedchainsNEW[i,10:12]=reord4[3,1:3]
    Ord_mcombinedchainsNEW[i,13:15]=reord5[3,1:3]
    Ord_mcombinedchainsNEW[i,16:18]=reord6[3,1:3]
    Ord_mcombinedchainsNEW[i,19:21]=reord7[3,1:3]
    Ord_mcombinedchainsNEW[i,22:24]=reord8[3,1:3]
    Ord_mcombinedchainsNEW[i,25:27]=reord9[3,1:3]
  }
  #
  Ord_mcombinedchainsNEW1=as.mcmc(Ord_mcombinedchainsNEW[1:20000,])  
  Ord_mcombinedchainsNEW2=as.mcmc(Ord_mcombinedchainsNEW[20001:40000,]) 
  Ord_mcombinedchainsNEW=as.mcmc.list(list(Ord_mcombinedchainsNEW1,Ord_mcombinedchainsNEW2))
  
  Ord_combinedchains<<-Ord_mcombinedchainsNEW
}