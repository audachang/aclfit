functions {
    real shifted_Wald_lpdf(real X, real gamma, real alpha, real theta){
        real tmp1;
        real tmp2;
        tmp1 = alpha / (sqrt(2 * pi() * (pow((X - theta), 3))));
        tmp2 = exp(-1 * (pow((alpha - gamma * (X-theta)),2)/(2*(X-theta))));
        return log(tmp1*tmp2) ;
    }
}
data{
    int t;
    vector<lower=0>[t] RT;
}
parameters {
    real<lower=0> gamma;
    real<lower=0> alpha;
    real<lower=0> theta;
}
transformed parameters{
    real<lower=0> E;
    real<lower=0> SD;
    E = alpha/gamma+theta;
    SD = sqrt(alpha / pow(gamma,3));
}
model{
    for(i in 1 : t)
        target+=shifted_Wald_lpdf(RT[i]|gamma,alpha,theta);
}
