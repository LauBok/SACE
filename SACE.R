SACE <- function(Z, S, Y, X, A, subset = NULL){
    #### Test Methods ####
    private.logical.test <- function(vec, name = NULL, interrupt = T){
        if (is.null(vec))
            return (T)
        if (!is.logical(vec) && !is.numeric(vec)){
            if (interrupt)
                stop(paste(name, "must be logical or numeric!", sep = " "))
            else return (F)
        }
        if (!is.logical(vec) && length(which(vec != 0 & vec != 1)) != 0)
            if (interrupt)
                warning(paste("Non 0/1 numbers in", name, "will be converted to TRUE.", sep = " "))
        return (T)
    }
    private.matrix.test <- function(mat, name = NULL, interrupt = T){
        if (is.null(mat))
            return (T)
        if (!is.matrix(mat) && !is.vector(mat)){
            if (interrupt)
                stop(paste(name, "should be a vector or a matrix", sep = " "))
            else return (F)
        }
        if (is.vector(mat))
            if (interrupt)
                warning(paste("Vector", name, "will be converted into a single column matrix.", sep = " "))
        return (T)
    }
    private.na.test <- function(mat, name = NULL, interrupt = T){
        if (is.null(mat))
            return (T)
        if (length(which(is.na(mat))) != 0){
            if (interrupt)
                stop(paste(name, "should not have NAs.", sep = " "))
            else return (F)
        }
        return (T)
    }
    
    #### Convert Methods ####
    private.logical.convert <- function(vec, name = NULL, interrupt = T){
        if (is.null(vec))
            return (NULL)
        private.logical.test(vec, name, interrupt)
        return (as.logical(vec))
    }
    private.matrix.convert <- function(mat, name = NULL, interrupt = T){
        if (is.null(mat))
            return (NULL)
        private.matrix.test(mat, name, interrupt)
        return (as.matrix(mat))
    }
    
    #### Get Methods ####
    Z.get <- function(){ private.Z }
    S.get <- function(){ private.S }
    Y.get <- function(){ private.Y }
    X.get <- function(){ private.X }
    A.get <- function(){ private.A }
    subset.get <- function(){ private.subset }
    
    private.mat.subset.get <- function(mat){
        if (is.null(mat))
            return (NULL)
        if (is.null(private.subset))
            return (mat)
        return (mat[private.subset,,drop = F])
    }
    
    Z.subset.get <- function(){ 
        if (private.checked) private.Z.sub.check
        else private.mat.subset.get(private.Z)
    }
    S.subset.get <- function(){ 
        if (private.checked) private.S.sub.check
        else private.mat.subset.get(private.S)
    }
    Y.subset.get <- function(){ 
        if (private.checked) private.Y.sub.check
        else private.mat.subset.get(private.Y)
    }
    X.subset.get <- function(){ 
        if (private.checked) private.X.sub.check
        else private.mat.subset.get(private.X)
    }
    A.subset.get <- function(){ 
        if (private.checked) private.A.sub.check
        else private.mat.subset.get(private.A)
    }
    
    #### Set Methods ####
    Z.set <- function(Z){ private.checked <<- FALSE; private.Z <<- private.matrix.convert(private.logical.convert(Z, "Z"), interrupt = F); }
    S.set <- function(S){ private.checked <<- FALSE; private.S <<- private.matrix.convert(private.logical.convert(S, "S"), interrupt = F); }
    Y.set <- function(Y){ private.checked <<- FALSE; private.Y <<- private.matrix.convert(Y, "Y"); }
    X.set <- function(X){ 
        private.checked <<- FALSE
        if (is.null(X)) private.X <<- NULL
        else private.X <<- private.matrix.convert(X, "X") 
        return (private.X)
    }
    A.set <- function(A){ 
        private.checked <<- FALSE
        if (is.null(A)) private.A <<- NULL
        else private.A <<- private.matrix.convert(A, "A") 
        return (private.A)
    }
    subset.set <- function(subset){ private.checked <<- FALSE; private.subset <<- private.logical.convert(subset, "subset"); }
    subset.unset <- function(){ private.subset <<- NULL; private.checked <<- FALSE }
    
    #### Check Method ####
    private.check <- function(interrupt){
        if (private.checked)
            return()
        ## size check ##
        if (ncol(private.Z) != 1 || ncol(private.S) != 1){
            stop("Z and S should be single-valued for each subject.")
        }
        if (is.null(private.Z) || is.null(private.S) || is.null(private.Y))
            stop("Z, S and Y are not nullable variables.")
        if (is.null(private.subset)){
            if (nrow(private.Z) != nrow(private.S))
                stop("Numbers of subjects in Z and S do not match.")
            if (nrow(private.Z) != nrow(private.Y))
                stop("Numbers of subjects in Z and Y do not match.")
            if (!is.null(private.X) && nrow(private.Z) != nrow(private.X))
                stop("Numbers of subjects in Z and X do not match.")
            if (!is.null(private.A) && nrow(private.Z) != nrow(private.A))
                stop("Numbers of subjects in Z and A do not match.")
        }
        else {
            if (nrow(private.Z) != length(private.subset))
                stop("The length of subset does not match the number of subjects.")
        }
        
        ## NA check ##
        private.na.test(Z.subset.get(), "Z")
        private.na.test(S.subset.get(), "S")
        private.na.test(X.subset.get(), "X")
        private.na.test(A.subset.get(), "A")
        if (!private.na.test(Y.subset.get()[S.subset.get() == 1], interrupt = F))
            stop("Y should not have missing values where S=1.")
        
        ## if successfully passed check ##
        private.Z.sub.check <<- Z.subset.get()
        private.S.sub.check <<- S.subset.get()
        private.Y.sub.check <<- Y.subset.get()
        private.X.sub.check <<- X.subset.get()
        private.A.sub.check <<- A.subset.get()
        private.checked <<- TRUE
        private.n.check <<- nrow(Z.subset.get())
        private.k.check <<- ncol(Y.subset.get())
        if (is.null(X.subset.get()))
            private.p.check <<- 0
        else
            private.p.check <<- ncol(X.subset.get())
        if (is.null(A.subset.get()))
            private.q.check <<- 0
        else
            private.q.check <<- ncol(A.subset.get())
    }
    
    #### MIE Method ####
    mie <- function(thres = 1e-6, optim.method = "BFGS", max.step = 1000, singular.ok = TRUE, need.variance = TRUE){
        private.check()
        if (is.null(A.subset.get())){
            warning("A is not provided, naive method (OLS) is used.")
            return(private.naive.mie())
        }
        private.model1.mie(thres, optim.method, max.step, singular.ok, need.variance)
    }
    
    private.naive.mie <- function(){
        W <- cbind(Z.subset.get(), X.subset.get())
        reg <- lm(Y.subset.get() ~ W)
        if (private.k.check == 1)
            sace <- reg$coefficients[2]
        else
            sace <- reg$coefficients[2,]
        return (list(CALL = match.call(), model = "naive",
                     n = private.n.check, sace = sace, reg = reg))
    }
    
    private.model1.mie <- function(thres = 1e-6, optim.method = "BFGS", max.step = 1000, singular.ok = TRUE, need.variance = TRUE){
        alpha1.estimate <- matrix(0, nrow = private.p.check + private.q.check + 1, ncol = private.k.check)
        alpha2.estimate <- matrix(0, nrow = private.p.check + 2, ncol = private.k.check)
        beta.estimate <- gamma.estimate <- numeric(private.p.check + private.q.check + 1)
        tildeX <- cbind(rep(1,private.n.check), X.subset.get())
        tildeW <- cbind(tildeX, A.subset.get())
        id_s1z1 <- S.subset.get() == T & Z.subset.get() == T
        id_s0z1 <- S.subset.get() == F & Z.subset.get() == T
        id_s1z0 <- S.subset.get() == T & Z.subset.get() == F
        id_s0z0 <- S.subset.get() == F & Z.subset.get() == F
        
        ## internal functions ##
        expit <- function(x) exp(x) / (1 + exp(x))
        evec <- function(vec) expit(tildeW %*% vec)
        sz.vec <- function(vec_s1z1, vec_s0z1, vec_s1z0, vec_s0z0) {
            vec_s1z1 * id_s1z1 + vec_s0z1 * id_s0z1 + vec_s1z0 * id_s1z0 + vec_s0z0 * id_s0z0
        }
        getbeta <- function(vec) vec[1:(length(vec)/2)]
        getgamma <- function(vec) vec[-(1:(length(vec)/2))]
        likelihood <- function(beta_gamma){
            ebeta <- evec(getbeta(beta_gamma))
            egamma <- evec(getgamma(beta_gamma))
            return(sum(sz.vec(log(ebeta),
                              log(1-ebeta),
                              log(ebeta*egamma),
                              log(1-ebeta*egamma))))
        }
        gradient <- function(beta_gamma){
            ebeta <- evec(getbeta(beta_gamma))
            egamma <- evec(getgamma(beta_gamma))
            zeros <- 0 * ebeta
            return(cbind(t(sz.vec(1 - ebeta,
                                  -ebeta,
                                  1 - ebeta,
                                  (1 - ebeta) / (1 - 1 / (ebeta * egamma)))) %*% tildeW,
                         t(sz.vec(zeros,
                                  zeros,
                                  1 - egamma,
                                  (1 - egamma) / (1 - 1 / (ebeta * egamma)))) %*% tildeW))
        }
        
        ## Estimation ##
        #---------------- beta & gamma ---------------#
        opt <- optim(c(beta.estimate, gamma.estimate), likelihood, gradient,
                     method = optim.method, hessian = need.variance,
                     control = list(fnscale = -1, maxit = 2 * max.step))
        beta.estimate <- getbeta(opt$par)
        gamma.estimate <- getgamma(opt$par)
        if (opt$convergence != 0) { 
            warning(paste("Optimization of beta and gamma didn't converge in", max.step, "steps!",sep = " "))
        }
        #------------------- alpha1 ------------------#
        lm.y.z0 <- lm(Y.subset.get() ~ tildeW - 1, subset = id_s1z0, singular.ok = singular.ok)
        alpha1.estimate <- lm.y.z0$coef
        names(alpha1.estimate) <- NULL
        rownames(alpha1.estimate) <- NULL
        #------------------- alpha2 ------------------#
        XL <- cbind(tildeX, evec(gamma.estimate))
        lm.y.z1 <- lm(Y.subset.get() ~ XL - 1, subset = id_s1z1, singular.ok = singular.ok)
        alpha2.estimate <- lm.y.z1$coef
        names(alpha2.estimate) <- NULL
        rownames(alpha2.estimate) <- NULL
        
        ## Prediction ##
        prob_W_raw <- evec(beta.estimate) * evec(gamma.estimate)
        prob_W <- t(prob_W_raw)/sum(prob_W_raw)
        #-------------------- mu0 --------------------#
        mu0_W <- tildeW %*% alpha1.estimate
        mu0 <- as.vector(prob_W %*% mu0_W)
        #-------------------- mu1 --------------------#
        mu1_W <- cbind(tildeX, rep(1, private.n.check)) %*% alpha2.estimate
        mu1 <- as.vector(prob_W %*% mu1_W)
        #------------------- sace --------------------#
        sace <- mu1 - mu0
        
        ## Variance ##
        if (need.variance){
            #----- !TODO -----#
        }
        results <- list(CALL = match.call(),
                        model = "model1",
                        n = private.n.check,
                        sace = sace,
                        options = list(optim.method = optim.method,
                                       need.variance=need.variance),
                        estimate = list(alpha1.estimate = alpha1.estimate,
                                        alpha2.estimate = alpha2.estimate,
                                        beta.estimate = beta.estimate,
                                        gamma.estimate = gamma.estimate,
                                        beta_gamma.convergence = opt$convergence,
                                        mu0.estimate = mu0,
                                        mu1.estimate = mu1))
        class(results) <- c("list", "mie")
        return (results)
    }
    
    #### Confident Interval ####
    private.boot.ci <- function(alpha = 0.05, nboot = 1000, max.step = 1000, singular.ok = TRUE, print.progress = TRUE, guarantee = T){
        options(warn = -1)
        private.check()
        i <- nskip <- 0
        sace.record <- matrix(nrow = 0, ncol = private.k.check)
        while(i < nboot) {
            i <- i + 1
            if(print.progress) cat("This is the", i, "th step.\n")
            bt.index <- sample(1:private.n.check, private.n.check, replace = TRUE)
            mie <- try(SACE(Z.subset.get()[bt.index,],
                            S.subset.get()[bt.index,],
                            Y.subset.get()[bt.index,],
                            X.subset.get()[bt.index,],
                            A.subset.get()[bt.index,])$MIE(max.step = max.step, singular.ok = singular.ok, need.variance = FALSE))
            if ('try-error' %in% class(mie) || is.na(sum(mie$sace))) { nskip <- nskip + 1; if (guarantee) i <- i - 1; next }
            else {
                sace.record <- rbind(sace.record, mie$sace)
            }
        }
        options(warn = 0)
        return(list(nskip = nskip, ci = apply(sace.record, 2, quantile, c(alpha / 2, 1 - alpha / 2), na.rm = T),
                    record = sace.record))
    }
    
    conf.int <- function(method = c("boot.ci"), alpha = 0.05, ...){
        if (method == "boot.ci")
            return (private.boot.ci(alpha,...))
        else if (method == "variance"){
            #----- !TODO -----#
        }
    }
    
    #### Generalized S3 Methods ####
    
    #### Variables ####
    if (missing(X)) X <- NULL
    if (missing(A)) A <- NULL
    private.Z <- Z.set(Z)
    private.S <- S.set(S)
    private.Y <- Y.set(Y)
    private.X <- X.set(X)
    private.A <- A.set(A)
    private.subset <- subset.set(subset)
    
    #### Status Variables ####
    private.checked <- FALSE
    
    #### Cache Variables (Should only be used when private.checked = TRUE) ####
    private.n.check <- private.k.check <- private.p.check <- private.q.check <- 0
    private.Z.sub.check <- private.S.sub.check <- private.Y.sub.check <- private.X.sub.check <- private.A.sub.check <- NULL
    
    return(list(GetZ = Z.get, GetS = S.get, GetY = Y.get, GetX = X.get, GetA = A.get, GetSubset = subset.get,
                GetZSubset = Z.subset.get, GetSSubset = S.subset.get, GetXSubset = X.subset.get, GetYSubset = Y.subset.get, GetASubset = A.subset.get,
                SetZ = Z.set, SetS = S.set, SetY = Y.set, SetX = X.set, SetA = A.set, SetSubset = subset.set, UnsetSubset = subset.unset,
                MIE = mie, ConfInt = conf.int))
}


expit <- function(x) { exp(x) / (1 + exp(x)) }
size <- 5000
X1 <- sample(c(1,-1),size,prob=c(0.5,0.5),replace = T)
mu <- c(1,-1)
sigma <- matrix(c(1,0.5,0.5,1),2,2)
X23 <- mvrnorm(size, mu, sigma)
X <- cbind(X1,X23)
A <- rbinom(size, 1, expit(X %*% c(0.5,0.5,0.5)))
Z <- rbinom(size, 1, expit(X %*% c(0.5,0.5,0.5) + 1))
S1 <- rbinom(size, 1, expit(X %*% c(0.5,0.5,0.5) + A + 2))
S0 <- rbinom(size, 1, expit(X %*% c(-1.5,0.5,0.5) + A))
S0[S1 == 0] <- 0
S <- S1 * Z + S0 * (1 - Z)
Y1LL <- rnorm(size, X %*% c(0.5,0.5,0.5), 0.5^2)
Y0LL <- rnorm(size, X %*% c(0.5,0.5,0.5) - 1, 0.5^2)
Y1LD <- rnorm(size, X %*% c(0.5,0.5,0.5) + 1, 0.5^2)
Y1 <- S0 * Y1LL + (1 - S0) * Y1LD
Y0 <- Y0LL
Y <- Y1 * Z + Y0 * (1 - Z)
Y[S == 0] <- NA

sace <- SACE(Z,S,Y,X,A)
print(sace$MIE())
print(sace$ConfInt(nboot=500,print.progress = F)$ci)
