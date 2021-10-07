//
// Two location Cormack Jolly Saber
data {
  int<lower=0> N; // number of individuals in sample
  int<lower=0> T; // number of sampling periods 
  vector<lower=0,upper=1>[5] y[N,T]; // capture histry for each individual 
}

transformed data{
  // some helpful vectors for calcualting
  // probaiblities of survival 
  vector[5] survival; 
  vector[5] unobserved; 
  
  // time of last observaiton 
  int<lower=0,upper=T> final[N];
  // total catch each time step 
  vector[T] total1; // location 1
  vector[T] total2; // location 2
  
  // compute final observation 
  for(i in 1:N){
    for(j in 1:T){
      if(sum(y[i,j]) > 0){
        final[i] = j;
      }
    }
  }

  // add ones to possible states
  survival[1] = 0;
  survival[2] = 1.0;
  survival[3] = 1.0;
  survival[4] = 0.0;
  survival[5] = 0.0;
  
  unobserved[1] = 1.0;
  unobserved[2] = 1.0;
  unobserved[3] = 1.0;
  unobserved[4] = 0.0;
  unobserved[5] = 0.0;
  
  // loop over years and sum accross y[N,T]
  for(i in 1:T){
    total1[i] = 0;
    total2[i] = 0;
    for(j in 1:N){
      total1[i] += y[j,i,4]; // capture in location 1?
      total2[i] += y[j,i,5]; // capture in location 2?
    }
  }
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector<lower=0,upper=1>[T-1] phi1; // survival probability location 1
  vector<lower=0,upper=1>[T-1] phi2; // survival probability location 1
  vector<lower=0,upper=1>[T-1] sigma12; // migraiton location 1 to location 2
  vector<lower=0,upper=1>[T-1] sigma21; // migraiton location 2 to location 1
  vector<lower=0,upper=1>[T] p1; // probabilithy of capture location 1
  vector<lower=0,upper=1>[T] p2; // probability of capture locaiton 2
}

transformed parameters {
  matrix<lower=0,upper=1>[5,5] transition[T-1];
  vector<lower=0, upper=1>[5] p_first_catch[T-1];
  vector<lower=0,upper=1>[N] probs;
  for(i in 1:(T-1)){
    // row 1
    transition[i,1,1] = 1.0;
    transition[i,1,2] = 1.0 - phi1[i];
    transition[i,1,3] = 1.0 - phi2[i];
    transition[i,1,4] = 1.0 - phi1[i];
    transition[i,1,5] = 1.0 - phi2[i];
    // row 2 
    transition[i,1,1] = 0.0;
    transition[i,1,2] = phi1[i]*(1-sigma12[i])*(1-p1[i+1]);
    transition[i,1,3] = phi2[i]*sigma21[i]*(1-p1[i+1]);
    transition[i,1,4] = phi1[i]*(1-sigma12[i])*(1-p1[i+1]);
    transition[i,1,5] = phi2[i]*sigma21[i]*(1-p1[i+1]);
    // row 3
    transition[i,1,1] = 0.0;
    transition[i,1,2] = phi1[i]*sigma12[i]*(1-p2[i+1]);
    transition[i,1,3] = phi2[i]*(1-sigma21[i])*(1-p2[i+1]);
    transition[i,1,4] = phi1[i]*sigma12[i]*(1-p2[i+1]);
    transition[i,1,5] = phi2[i]*(1-sigma21[i])*(1-p2[i+1]);
    // row 4
    transition[i,1,1] = 0.0;
    transition[i,1,2] = phi1[i]*(1-sigma12[i])*p1[i+1];
    transition[i,1,3] = phi2[i]*sigma21[i]*p1[i+1];
    transition[i,1,4] = phi1[i]*(1-sigma12[i])*p1[i+1];
    transition[i,1,5] = phi2[i]*sigma21[i]*p1[i+1];
    // row 5
    transition[i,1,1] = 0.0;
    transition[i,1,2] = phi1[i]*sigma12[i]*p2[i+1];
    transition[i,1,3] = phi2[i]*(1-sigma21[i])*p2[i+1];
    transition[i,1,4] = phi1[i]*sigma12[i]*p2[i+1];
    transition[i,1,5] = phi2[i]*(1-sigma21[i])*p2[i+1];
    
    // probability of first catch 
    p_first_catch[i,1] = 0;
    p_first_catch[i,2] = 0;
    p_first_catch[i,3] = 0;
    p_first_catch[i,4] = p1[i];
    p_first_catch[i,5] = p2[i];
  }
  
  for(i in 1:N){
    real p;
    int counter;
    vector[5] state;
    p = 1;
    counter = 0;
    for(j in 1:T){
      
      // initialize state after first observaiton
      // increment counter 
      if(sum(y[i,j]) > 0.0 && counter == 0){
        state = transition[j]*y[i,j]; // update state 
        counter += 1;
        
      // increment probability and update state
      // for later observaitons 
      } else if(sum(y[i,j]) > 0.0 && counter > 0 ){
        p *= dot_product(state, y[i,j]);
        state = transition[j]*y[i,j]; // update state 
        
      // if not observed in present but observed in the 
      // future update state and add proability of survival
      } else if(counter > 0 && j < final[i]){
        p*= dot_product(state, survival);
        state = transition[j]*state; // update state
        
      // after last observaiton add probability 
      // that the individual is left the popuatlion
      // or that they are not captured
      } else if(j > final[i]){
        p*= dot_product(state, unobserved);
        state = transition[j]*state; // update state
      }
      
    }
    probs[i] = p;
  }
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for(i in 1:N){
    target += log(probs[i]);
  }
}

generated quantities {
  vector[T] Nt1; // abundnace location 1
  vector[T] Nt2; // abundnace location 2
  for(i in 1:T){
    Nt1[i] = total1[i] / p1[i]; // Total catch divided 
    Nt2[i] = total2[i] / p2[i]; // by  prob of catch.
  }
}



