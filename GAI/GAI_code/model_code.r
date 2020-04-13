# This file contains the functions required fit the GAI model approach


# Logit link functions
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

library(splines)

# Likelihood function
ll_func <- function(parm,irep=1,Nguess=NULL){
	if(a.type == "N" | a.type == "SO"){
		par.index <- 0
		# mu can be constant or varying with a spatial covariate
		mu.int <- parm[par.index+1]; par.index <- par.index + 1
		if(mu.type == "cov"){mu.slope <- parm[par.index+1]; par.index <- par.index + 1}
		switch(mu.type,
			common = {mu.est1 <- rep(exp(mu.int),nS)},
			cov = {mu.est1 <- exp(mu.int + mu.slope*mu.cov)})
	
		switch(B,
			"1" = {
				sigma.est <- rep(exp(parm[par.index+1]),2); par.index <- par.index + 1},
			"2" = {
				# mu_d can be constant or varying with a spatial covariate
				mu.diff.int <- parm[par.index+1]; par.index <- par.index + 1
				if(mu.diff.type == "cov"){mu.diff.slope <- parm[par.index+1]; par.index <- par.index + 1} 
				switch(mu.diff.type,
					common = {mu.diff.est <- rep(exp(mu.diff.int),nS)},
					cov = {mu.diff.est <- exp(mu.diff.int + mu.diff.slope*mu.diff.cov)})	
				# if B = 2, sigma can be the same or different for each brood
				switch(sigma.type,
					hom = {sigma.est <- rep(exp(parm[par.index+1]),2); par.index <- par.index + 1},
					het = {sigma.est <- exp(parm[(par.index+1):(par.index+2)]); par.index <- par.index + 2})
				# w can be constant or varying with a spatial covariate
				w.int <- parm[par.index+1]; par.index <- par.index + 1
				if(w.type == "cov"){w.slope <- parm[par.index+1]; par.index <- par.index + 1}
				switch(w.type,
					common = {w.est <- rep(expit(w.int),nS)},
					cov = {w.est <- expit(w.int + w.slope*w.cov)})})
		
		switch(a.type,
			N = {
				switch(B,
					"1" = {
						afunc <- matrix(dnorm(rep(1:nT,each=nS),mu.est1,sigma.est),nrow=nS,ncol=nT)},
					"2" = {
						afunc <- matrix(rep(w.est,nT)*dnorm(rep(1:nT,each=nS),mu.est1,sigma.est[1]) + (1-rep(w.est,nT))*dnorm(rep(1:nT,each=nS),mu.est1+mu.diff.est,sigma.est[2]),nrow=nS,ncol=nT)})},
			SO = {
				# We assume constant survival phi, but this could be adapted
				phi.est <-  expit(parm[par.index+1]); par.index <- par.index + 1
				switch(B,
					"1" = {	
						betta.est <- matrix(c(pnorm(1,mean=mu.est1,sd=sigma.est[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.est1,length(2:(nT-1))),sd=sigma.est[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.est1,length(1:(nT-2))),sd=sigma.est[1]),1-pnorm(nT-1,mean=mu.est1,sd=sigma.est[1])),nrow=nS,ncol=nT)},
					"2" = {
						betta.est <- matrix(rep(w.est,each=nT)*c(pnorm(1,mean=mu.est1,sd=sigma.est[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.est1,length(2:(nT-1))),sd=sigma.est[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.est1,length(1:(nT-2))),sd=sigma.est[1]),1-pnorm(rep(nT-1,each=nS),mean=mu.est1,sd=sigma.est[1])) + (1-rep(w.est,each=nT))*c(pnorm(1,mean=mu.est1 + mu.diff.est,sd=sigma.est[2]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu.est1+ mu.diff.est,length(2:(nT-1))),sd=sigma.est[2])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu.est1+ mu.diff.est,length(1:(nT-2))),sd=sigma.est[2]),1-pnorm(nT-1,mean=mu.est1+ mu.diff.est,sd=sigma.est[2])) ,nrow=nS,ncol=nT)})
						
				afunc <- betta.est				
				for(j in 2:nT){
					for(b in 1:(j-1)){
						afunc[,j] <- afunc[,j] + betta.est[,b]*phi.est^length(b:(j-1))
						}
					}})
		afunc[is.na(y)] <- NA
		
	} else if (a.type == "S"){
		if(dist.type == "P"){alpha.est <- parm} else {alpha.est <- parm[1:(length(parm)-1)]; par.index <- length(parm)-1}
		bsbasis <- bs(1:nT,df=degf,degree=deg,intercept=TRUE)
		afunc <- exp(matrix(bsbasis%*%alpha.est,nrow=nS,ncol=nT,byrow=TRUE))
		afunc <- afunc/rowSums(afunc)
		afunc[is.na(y)] <- NA
		}
	
	if(dist.type == "NB"){
		r.est <- exp(parm[par.index+1])}
	if(dist.type == "ZIP"){
		psi.est <- expit(parm[par.index+1])}
	
	# Concentrated likelihood formulation
	switch(dist.type, 
		P = { 
			llik <- dpois(y,lambda=afunc*rep(apply(y,1,sum,na.rm=TRUE)/apply(afunc,1,sum,na.rm=TRUE),nT),log=TRUE)},
		NB = { 
			if(irep > 1){
			llik <- dnbinom(y,mu=afunc*matrix(Nguess,nrow=nS,ncol=nT),size=r.est,log=TRUE)
			} else {
			llik <- dnbinom(y,mu=afunc*rep(apply(y,1,sum,na.rm=TRUE)/apply(afunc,1,sum,na.rm=TRUE),nT),size=r.est,log=TRUE)}}, 
		ZIP = {
			if(irep > 1){
			llik <- dpois(y,lambda=afunc*matrix(Nguess,nrow=nS,ncol=nT),log=FALSE)
			llik[y==0 & !is.na(y)] <- log((1-psi.est) + psi.est*llik[y==0 & !is.na(y)])
			llik[y!=0 & !is.na(y)] <- log(psi.est*llik[y!=0 & !is.na(y)])
			} else {
			llik <- dpois(y,lambda=afunc*rep(apply(y,1,sum,na.rm=TRUE)/apply(afunc,1,sum,na.rm=TRUE),nT),log=FALSE)
			llik[y==0 & !is.na(y)] <- log((1-psi.est) + psi.est*llik[y==0 & !is.na(y)])
			llik[y!=0 & !is.na(y)] <- log(psi.est*llik[y!=0 & !is.na(y)])}})
	-1*sum(llik,na.rm=TRUE) 
	}
    

#Starting values   
start_val_func <- function(){
	psi.st <- NULL
	if(dist.type == "ZIP"){
		psi.st <- logit(sample(seq(0.5,0.9,0.1),1))}	
	r.st <- NULL
	if(dist.type == "NB"){
		r.st <- log(sample(1:5,1))}
	
	if(a.type == "S"){
		parm <- c(psi.st,sample(seq(-2,2,.1),degf),r.st)
	} else {
		samp <- sort(sample(5:15,2))
		mu1.int.st <- log(samp[1])
		mu1.slope.st <- NULL
		if(mu.type == "cov") mu1.slope.st <- 0
		
		mu.diff.int.st <- mu.diff.slope.st <- NULL
		if(B == 2){
			mu.diff.int.st <- log(samp[2]-samp[1])
			if(exp(mu.diff.int.st)<7)mu.diff.int.st <- log(7)
			if(mu.diff.type == "cov") mu.diff.slope.st <- 0
			}
			
		if(B == 1){
			sigma.st <- log(sample(2:3,1))
			w.int.st <- w.slope.st <- NULL
		} else {
			switch(sigma.type,
				hom = {sigma.st <- log(sample(2:3,1))},
				het = {sigma.st <- rep(sample(2:3,1),2)})
			w.int.st <- logit(sample(seq(.2,.8,.1),1))
			w.slope.st <- NULL
			if(w.type == "cov") w.slope.st <- 0
			}
		switch(a.type,
			N = {
				parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st)},
			SO = {
				phi.st <- logit(sample(seq(.3,.9,.1),1))
				parm <- c(mu1.int.st,mu1.slope.st,mu.diff.int.st,mu.diff.slope.st,sigma.st,w.int.st,w.slope.st,r.st,psi.st,phi.st)})
		}		
	return(parm)
	}		
	

# Equation 4 in the paper
der_func <- function(x,model,isite,bmat){
	switch(dist.type,
		NB = {
			sum(model$modelfit$y[isite,]/x - (model$r.out+model$modelfit$y[isite,])*model$afunc.out[isite,]/(model$r.out+x*model$afunc.out[isite,]),na.rm=TRUE)},
		ZIP = {
			sum((1-bmat[isite,])*(-model$afunc.out[isite,]*model$psi.est[1]*exp(-model$afunc.out[isite,]*x))/(1-model$psi.est[1]+model$psi.est[1]*exp(-model$afunc.out[isite,]*x)) - bmat[isite,]*model$afunc.out[isite,] + bmat[isite,]*model$modelfit$y[isite,]/x,na.rm=TRUE)})
	}	
	
# Wrapper for fitting GAI the model for multiple starts, including the iterative approach for NB and ZIP	
fit_it_model <- function(){
	if(!(dist.type %in% c("P","NB","ZIP")))
		stop("Distribution must be P, NB or ZIP")
	if(!(a.type %in% c("N","SO","S")))
		stop("Function for {a} must be N, SO or S")

	fit_k <- list(); fit_k.ll <- rep(NA,nstart)
	for(k in 1:nstart){	
		st <- proc.time()
		irep <- 1
		fit1 <- try(fit_model(irep=irep),silent=FALSE)
		# If dist.type is "NB" or "ZIP" the iterative procedure is required
		if(dist.type %in% c("NB","ZIP")){
			if(dist.type == "ZIP"){
				# A matrix b indicating where y_{i,j} > 0
				bmat <- matrix(1,nrow=nS,ncol=nT)	
				bmat[is.na(y)] <- NA
				bmat[!is.na(y) & y==0] <- 0
				}
			lld <- 1
			fit <- list()
			fit[[1]] <- fit1
			ll <- NA
			ll[1] <- fit1$ll.val
			uppvals <- lowvals <- NULL
			# Iterate until convergence (here defined when the difference in likelihoods is sufficiently small)
			while(lld > 0.001){
				irep <- irep + 1
				Nest <- rep(NA,nS)
				for(isite in 1:nS){
					if(max(y[isite,],na.rm=TRUE) == 0){
						low <- 0
					} else {
						low <- 0.1
						}
					lowvals <- c(lowvals,low)
					uppvals <- c(uppvals,2500)
					# Find each N_i numerically
					temp <- try(uniroot(der_func,lower=low,upper=2500,model=fit1,isite=isite,bmat=bmat)$root,silent=TRUE)
					dtemp <- 1
					# Different options for the interval in the 1-d root finding
					lowval <- c(0,0,0,0,.1,.1,.1,.1,.01,.01,.01,.01,.01)
					uppval <- c(1000,2500,5000,10000,1000,2500,5000,10000,1000,2500,5000,10000)
					while(class(temp) == "try-error" & dtemp < 9){
						temp <- try(uniroot(der_func,lower=lowval[dtemp],upper=uppval[dtemp],bmat=bmat,model=fit1,isite=isite)$root,silent=TRUE)
						dtemp <- dtemp + 1
						lowvals <- c(lowvals,lowval)
						uppvals <- c(uppvals,uppval)
						}
					Nest[isite] <- unlist(temp)
					}
				vals <- fit1$modelfit$allval
				fit[[irep]] <- try(fit_model(irep=irep,Nguess=Nest,vals=vals),silent=FALSE)
				ll[irep] <- fit[[irep]]$ll.val
				fit1 <- fit[[irep]]
				lld <- abs(ll[irep]-ll[irep-1])
				}
			fit1 <- fit[[irep]]
			fit1$iterations <- list(fit)
			}
			et <- proc.time()
			fit1$time <- (et-st)[3]
			fit_k[[k]] <- fit1
			fit_k.ll[k] <- fit_k[[k]]$ll.val	
		}
	output <- list(fit_k[[min(c(1:nstart)[fit_k.ll==max(fit_k.ll,na.rm=T)])]],fit_k,fit_k.ll)	
	return(output)		
	}
	
# Fit the GAI model
fit_model <- function(irep=1,Nguess=NULL,vals=NULL){
	
	if(irep==1){parm <- start_val_func()} else {parm <- vals}	
		
	#if(a.type == "S") { meth <- "BFGS" } else {meth <- "Nelder-Mead"}
	meth <- "Nelder-Mead"

	this.fit <- try(optim(par=parm,
					fn=ll_func,Nguess=Nguess,irep=irep,hessian=TRUE,method=meth,
					control=list(maxit=100000)),silent=TRUE)
			
	if(is.list(this.fit) & class(try(solve(this.fit$hessian),silent=TRUE)) != "try-error"){	
		# Model output
		N.out <- psi.out <- r.out <- mu1.out <- mu1.int.out <- mu1.slope.out <- mu.diff.out <- mu.diff.int.out <- mu.diff.slope.out <- w.out <- w.int.out <- w.slope.out <- sigma.out <- alpha.out <- phi.out <- betta.out <- bsbasis <- NULL
	
		out.index <- 0
	
		if(a.type == "N" | a.type == "SO"){
			alpha.out <- NULL
			switch(mu.type,
				common = {mu1.out <- rep(exp(this.fit$par[out.index+1]),nS); mu1.int.out <- this.fit$par[out.index+1]; out.index <- out.index + 1},
				cov = {mu1.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*mu.cov); mu1.int.out <- this.fit$par[out.index+1]; mu1.slope.out <- this.fit$par[out.index+2]; out.index <- out.index + 2})
			
			switch(B,
				"1" = {
					sigma.out <- rep(exp(this.fit$par[out.index + 1]),2); out.index <- out.index + 1}, 
				"2" = {
					switch(mu.diff.type,
						common = {mu.diff.out <- rep(exp(this.fit$par[out.index + 1]),nS); mu.diff.int.out <- this.fit$par[out.index+1]; out.index <- out.index + 1},
						cov = {mu.diff.out <- exp(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*mu.diff.cov); mu.diff.int.out <- this.fit$par[out.index+1]; mu.diff.slope.out <- this.fit$par[out.index+2]; out.index <- out.index + 2})
			
					switch(sigma.type,
						hom = {sigma.out <- rep(exp(this.fit$par[out.index + 1]),2); out.index <- out.index + 1},
						het = {sigma.out <- exp(this.fit$par[(out.index + 1):(out.index + 2)]); out.index <- out.index + 2})
					
					switch(w.type,
						common = {w.out <- rep(expit(this.fit$par[out.index + 1]),nS); w.int.out <- this.fit$par[out.index+1];out.index <- out.index + 1},
						cov = {w.out <- expit(this.fit$par[out.index + 1] + this.fit$par[out.index + 2]*w.cov); w.int.out <- this.fit$par[out.index+1]; w.slope.out <- this.fit$par[out.index+2];out.index <- out.index + 2})})
					
			if(a.type == "SO"){
				phi.out <-  expit(this.fit$par[out.index+1]); out.index <- out.index + 1
				}			
		} else if(a.type == "S"){
			if(dist.type == "P"){alpha.out <- this.fit$par} else {alpha.out <- this.fit$par[1:(length(this.fit$par)-1)]; out.index <- length(this.fit$par)-1}
			bsbasis <- bs(1:nT,df=degf,degree=deg,intercept=TRUE)
			afunc.out <- exp(matrix(bsbasis%*%alpha.out,nrow=nS,ncol=nT,byrow=TRUE))
			afunc.out <- afunc.out/rowSums(afunc.out)	
		}

		if(dist.type == "NB"){
			r.out <- exp(this.fit$par[out.index+1])}
		if(dist.type == "ZIP"){
			psi.out <- expit(this.fit$par[out.index+1])
			}
	
		switch(a.type,
			N = {
				switch(B,
					"1" = {
						afunc.out <- matrix(dnorm(rep(1:nT,nS),mu1.out,sigma.out),nrow=nS,ncol=nT,byrow=TRUE)},
					"2" = {
						afunc.out <- matrix(rep(w.out,nT)*dnorm(rep(1:nT,nS),mu1.out,sigma.out[1]) + (1-rep(w.out,nT))*dnorm(rep(1:nT,nS),mu1.out + mu.diff.out,sigma.out[2]),nrow=nS,ncol=nT,byrow=TRUE)})},
			SO = {
				switch(B,
					"1" = {	
						betta.out <- matrix(c(pnorm(1,mean=mu1.out,sd=sigma.out[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu1.out,length(2:(nT-1))),sd=sigma.out[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu1.out,length(1:(nT-2))),sd=sigma.out[1]),1-pnorm(nT-1,mean=mu1.out,sd=sigma.out[1])),nrow=nS,ncol=nT)},
					"2" = {
						betta.out <- matrix(rep(w.out,nT)*c(pnorm(1,mean=mu1.out,sd=sigma.out[1]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu1.out,length(2:(nT-1))),sd=sigma.out[1])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu1.out,length(1:(nT-2))),sd=sigma.out[1]),1-pnorm(nT-1,mean=mu1.out,sd=sigma.out[1])) + (1-rep(w.out,nT))*c(pnorm(1,mean=mu1.out + mu.diff.out,sd=sigma.out[2]),pnorm(rep(2:(nT-1),each=nS),mean=rep(mu1.out+ mu.diff.out,length(2:(nT-1))),sd=sigma.out[2])-pnorm(rep(1:(nT-2),each=nS),mean=rep(mu1.out+ mu.diff.out,length(1:(nT-2))),sd=sigma.out[2]),1-pnorm(nT-1,mean=mu1.out+ mu.diff.out,sd=sigma.out[2])) ,nrow=nS,ncol=nT)})
						
				afunc.out <- betta.out				
				for(j in 2:nT){
					for(b in 1:(j-1)){
							afunc.out[,j] <- afunc.out[,j] + betta.out[,b]*phi.out^length(b:(j-1))
						}
					}})
			
		afunc.outNA <- afunc.out
		afunc.outNA[is.na(y)] <- NA
		N.out <- apply(y,1,sum,na.rm=TRUE)/apply(afunc.outNA,1,sum,na.rm=TRUE)

	
		output <- list(ll.val=-this.fit$value,
					npar=length(this.fit$par),
					N.est=N.out,
					w.est=w.out,
					w.int=w.int.out,
					w.slope=w.slope.out,
					mu1.est=mu1.out,
					mu1.int=mu1.int.out,
					mu1.slope=mu1.slope.out,
					mu.diff.est=mu.diff.out,mu.diff.int=mu.diff.int.out,
					mu.diff.slope=mu.diff.slope.out,
					sigma.out=sigma.out,
					r.out=r.out,
					psi.est=psi.out,
					phi.out=phi.out,
					afunc.out=afunc.out,
					afunc.outNA=afunc.outNA,
					betta.out=betta.out,
					bsbasis=bsbasis,
					alpha.out=alpha.out,
					modelfit=list(Hessian=this.fit$hessian,
								  starts=parm,
								  allval=this.fit$par,
								  nS=nS,nT=nT,
								  convergence=this.fit$convergence,
								  y=y))
						  
		output$Fitted <- output$afunc.out*output$N.est
		if(dist.type == "ZIP") output$Fitted <- output$psi.est*output$Fitted
		output$dev <- switch(dist.type,
						P={
							2*(sum((output$modelfit$y*log(output$modelfit$y/output$Fitted)-(output$modelfit$y-output$Fitted))[!is.na(output$modelfit$y) & output$modelfit$y!=0])+sum(output$Fitted[output$modelfit$y== 0 & !is.na(output$modelfit$y)]))},
						NB={
							2*(sum((output$modelfit$y*log(output$modelfit$y/output$Fitted)-((output$r.out+output$modelfit$y)*log((output$r.out+output$modelfit$y)/(output$r.out+output$Fitted))))[output$modelfit$y>0 & !is.na(output$modelfit$y)])-sum((output$r.out*log(output$r.out/(output$r.out+output$Fitted)))[!is.na(output$modelfit$y) & output$modelfit$y==0]))},
						ZIP={NA})
		output$D <- output$dev/(length(output$modelfit$y[!is.na(output$modelfit$y)])-output$npar)
		output
	} else {NA}
	}	
		
	
	