// Generalized normal (aka power-exponential) unnormalized log-probability
// Vectorized version
real gnormal_lpdf(vector y, real mu, real sigma, real shape) {
  vector[num_elements(y)] LL;
  
  for(i in 1:num_elements(LL))
    LL[i] = -pow(abs(y[i] - mu)/sigma, shape);
  
  return(sum(LL));
}
