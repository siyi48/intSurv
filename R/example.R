# data preparation
q.c1 <- 3
q.c0 <- 3
q.lambda <- 3
nfolds <- 5
c1.mat <- matrix(dat$x1, ncol = 1)
c1.discrete <- matrix(dat$x2, ncol = 1)
c0.mat <- matrix(dat$x1, ncol = 1)
c0.discrete <- matrix(dat$x2, ncol = 1)
lambda.mat <-  matrix(dat$x1, ncol = 1)
lambda.discrete <- cbind(1, dat$x2)
tau.mat <- cbind(1, dat$x1, dat$x1^2)
s <- dat$s
delta <- dat$delta
t <- dat$t
a <- dat$a

# run the function
res.hte <- surv.hte(tau.mat = tau.mat, c1.mat = c1.mat,
                    c1.discrete = c1.discrete, q.c1 = q.c1,
                    c0.mat = c0.mat, c0.discrete = c0.discrete, q.c0 = q.c0,
                    lambda.mat = lambda.mat,
                    lambda.discrete = lambda.discrete, q.lambda = q.lambda,
                    a = a, s = s, t = t, delta = delta,
                    nfolds = nfolds, type = c("int", "rct"))
res.hte
