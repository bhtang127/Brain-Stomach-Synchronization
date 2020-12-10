## Hilbert transformation (to analytic signal)
HT = function(sig){
  ndata = length(sig)
  h = rep(0, ndata)
  if(ndata %% 2 == 0){
    h[c(1, ndata/2+1)] = 1 
    h[2:(ndata/2)] = 2 
  }
  else{
    h[1] = 1
    h[2:((ndata + 1)/2)] = 2 
  }
  fft(h * fft(sig), inverse = TRUE)/ndata
}

## PLV value (time averaged)
PLV = function(sig1, sig2){
  arg1 = Arg(HT(sig1))
  arg2 = Arg(HT(sig2))
  Mod(mean(exp(1i*(arg1-arg2))))
}


library(R.matlab)

## Load full data
full_data_ds = readMat("GICA_EGGRSN_analytic_ds.mat")

## Extract de-spiked Hilbert transformed signal
analytic_ds = full_data_ds$EGGRSN.analytic[[2]]

## Extract region name
region_name = c()
for(i in 1:18){
  region_name[i] = c(full_data_ds$EGGRSN.analytic[[1]][[2]][[i]][[1]])
}
region_rank = c(full_data_ds$EGGRSN.analytic[[1]][[3]])
#### Correct visual region name
region_name[region_name=="VIS_b"] = "VIS_a"
region_name[region_name=="VIS_c"] = "VIS_b"


## Function to read signal with specific 
## session, run and component
## Component 0 for EGG and 1-18 for 18 resting states fmri
readMRIEGG = function(session, run, component=0, data=analytic_ds){
  if(component == 0){
    data[[session]][[run*2-1]][,1]
  }
  else{
    data[[session]][[run*2]][,component]
  }
}
