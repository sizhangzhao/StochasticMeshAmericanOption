rm(list = ls())

#This one is for the control variate and the antithetic variate

#=========  Parameter Specification  =========
S_0 = 110       #initial price
K = 100        #strike price
TT = 3        #maturity
n = 5          #number of assets
N = 500        #simulation times
b = 50        #mesh size
r = 0.05       #annual interest rate
sigma = 0.2    #volatility
dt = 1         #time interval
delta = 0.1    #dividend rate
sigma_ = sigma * sqrt(dt)                 #for density function calculation
mu_ = (r - delta - 0.5 * sigma^2) * dt    #for density function calculation


#=========  Mesh Construction  ============

simu = function(S_0, TT, n, b, r, sigma, dt, delta, two){
  if(two == FALSE){
    S_1toT = list()
    tem = matrix(rnorm(b*n,mean=0, sd = sqrt(dt)), b, n)
    S_1toT[[1]] = S_0 + (r - delta)*S_0*dt + sigma*S_0*tem
    while(any(S_1toT[[1]] <= 0)){
      tem = matrix(rnorm(b*n,mean=0, sd = sqrt(dt)), b, n)
      S_1toT[[1]] = S_0 + (r - delta)*S_0*dt + sigma*S_0*tem
    }
    for (i in 2:TT){
      tem = matrix(rnorm(b*n,mean=0, sd = sqrt(dt)), b, n)
      S_1toT[[i]] = S_1toT[[(i-1)]] * matrix((1 + (r - delta)*dt + sigma*tem), b, n)
      while(any(S_1toT[[i]] <= 0)){
        tem = matrix(rnorm(b*n,mean=0, sd = sqrt(dt)), b, n)
        S_1toT[[i]] = S_1toT[[(i-1)]] * matrix((1 + (r - delta)*dt + sigma*tem), b, n)
      }
    }
  }
  else{
    S_1toT = list()
    path1 = list() 
    path2 = list()
    tem = matrix(rnorm(b*n,mean=0, sd = sqrt(dt)), b, n)
    path1[[1]] = S_0 * (1 + (r - delta) * dt + sigma * tem)
    path2[[1]] = S_0 * (1 + (r - delta) * dt - sigma * tem)
    while(any(path1[[1]] <= 0) || any(path2[[1]] <= 0)){
      tem = matrix(rnorm(b*n,mean=0, sd = sqrt(dt)), b, n)
      path1[[1]] = S_0 * (1 + (r - delta) * dt + sigma * tem)
      path2[[1]] = S_0 * (1 + (r - delta) * dt - sigma * tem)
    }
    for (i in 2:TT){
      tem = matrix(rnorm(b*n,mean=0, sd = sqrt(dt)), b, n)
      path1[[i]] = path1[[(i-1)]] * matrix((1 + (r - delta)*dt + sigma*tem), b, n)
      path2[[i]] = path2[[(i-1)]] * matrix((1 + (r - delta)*dt - sigma*tem), b, n)
      while(any(path1[[i]] <= 0) || any(path2[[i]] <= 0)){
        tem = matrix(rnorm(b*n,mean=0, sd = sqrt(dt)), b, n)
        path1[[i]] = path1[[(i-1)]] * matrix((1 + (r - delta)*dt + sigma*tem), b, n)
        path2[[i]] = path2[[(i-1)]] * matrix((1 + (r - delta)*dt - sigma*tem), b, n)
      }
    }
    S_1toT[[1]] = path1
    S_1toT[[2]] = path2
  }
  S_1toT
}

#mesh = simu(S_0, TT, n, 1, r, sigma, dt, delta, TRUE)

#=========  Weight Calculation  ============
dens = function(prev, mu, sigma, x){
  mu = log(prev) + mu
  dens = exp(-(log(x) - mu)^2/(2 * sigma^2))/(x * sigma * sqrt(2 * pi))
  dens
}

weight = function(curr, xt1, b, n, mu, sigma){
  curr_f = c()
  for(i in 1:n){
    curr_f = cbind(curr_f, sapply(xt1[,i], dens, mu = mu, sigma = sigma, x = curr[i], simplify = TRUE))
  }
  curr_p = apply(curr_f, 1, prod)
  denom = sum(curr_p) / b
  weight = curr_p / denom
  weight
}

weightall = function(mesh, b, n, T, mu, sigma){
  weightall = list()
  for(i in 1:(T-1)){
    curr = mesh[[(i+1)]]
    prev = mesh[[i]]
    tem = c()
    for (j in 1:b){
      tem = cbind(tem, weight(curr[j,], prev, b, n, mu, sigma)) 
    }
    weightall[[i]] = tem
  }
  weightall
}

#wt = weightall(mesh, b, n, TT, mu_, sigma_) ###Just for test

#=============   Mesh estimator  ===============
intrinsic = function(state, k){
  value = max(max(state) - k, 0) 
  value
}

Q_mat = function(mesh, b, K, TT, wt){
  Q_1toT = matrix(nrow = b, ncol = TT)
  for (i in 1:b){
    Q_1toT[i,TT] = intrinsic(mesh[[TT]][i,], K)
  }
  for (t in (TT-1):1){
    for (i in 1:b){
      Q_1toT[i,t] = max(intrinsic(mesh[[t]][i,], K), 1/b *  Q_1toT[,(t+1)]%*% wt[[t]][i,] * exp(-r * dt) ) # * exp(-r * dt) 
    }
  }
  Q_1toT
}

mesh_est = function(Q_1toT, mesh, b, K){
  mesh_est = max(intrinsic(S_0,K), 1/b * sum(Q_1toT[,1])* exp(-r * dt) )
  mesh_est
}

#mesh_est(Q_mat, mesh, b, K) #just for test


#===========   Path Estimator    =================
path_estimator = function(mesh, q, n, TT, S_0, r, sigma, dt, delta, k, mesh_est, path){
  sigma_ = sigma * sqrt(dt)
  mu_ = (r - delta - 0.5 * sigma^2) * dt
  find = FALSE
  time = TT
  if(intrinsic(S_0, k) >= mesh_est){
    find = TRUE
    time = 0
  }
  if(find == FALSE){
    for(i in 1:(TT-1)){
      tomorrow = mesh[[i+1]]
      today = rbind(mesh[[i]], path[[i]])
      intrin = intrinsic(path[[i]], k)
      wt = c()
      for (j in 1:b){
        wt = cbind(wt, weight(tomorrow[j,], today, b, n, mu_, sigma_)) 
      }
      wts = wt[(b+1),]
      nextq = q[,(i+1)]
      value = wts %*% nextq * exp(-r * dt) / b
      if(intrin >= value){
        find = TRUE
        time = dt * i
        break
      }
    }
  }
  if(!find){
    intrin = intrinsic(path[[TT]], k)
  }
  intrin * exp(-r * time)
}

#========= Simulation for both mesh and path ============

meshestimator = c()
pathestimator = c()
two = FALSE

for (i in 1:N){
	mesh = simu(S_0, TT, n, b, r, sigma, dt, delta, FALSE)
	wt = weightall(mesh, b, n, TT, mu_, sigma_)
	q = Q_mat(mesh, b, K, TT, wt)
  meshest = mesh_est(q, mesh, b, K)
  meshestimator = c(meshestimator, meshest)
  if(two == TRUE){
    path = simu(S_0, TT, n, 1, r, sigma, dt, delta, two)
  	path_est1 = path_estimator(mesh, q, n, TT, S_0, r, sigma, dt, delta, K, meshest, path[[1]])
    path_est2 = path_estimator(mesh, q, n, TT, S_0, r, sigma, dt, delta, K, meshest, path[[2]])
    path_est = (path_est1+path_est2) / 2
    pathestimator = c(pathestimator, path_est)
  }
  else{
    path = simu(S_0, TT, n, 1, r, sigma, dt, delta, two)
    path_est = path_estimator(mesh, q, n, TT, S_0, r, sigma, dt, delta, K, meshest, path)
    pathestimator = c(pathestimator, path_est)
  }
}

varp = var(pathestimator)
varm = var(meshestimator)
meanmesh = mean(meshestimator)
meanpath = mean(pathestimator)
sdmesh = sd(meshestimator)
sdpath = sd(pathestimator)
widmesh = sdmesh * 1.96 / sqrt(N-1)
widpath = sdpath * 1.96 / sqrt(N-1)
