## pipe.design version 0.1
## MJS 18/03/2015

if(getRversion() >= "2.15.1") globalVariables(c("y","z"))


## N - sample size
## S - Number of trial simulations (default 1)
## c - cohort size
## theta - TTL
## pi - true p(DLT) matrix. If not specified only next recommended dose is returned.
## prior.med - prior median p(DLT) matrix
## prior.ss - prior sample size matrix
## strategy - dose escalation strategy. Options are 
	## "ss" - Sample size strategy. Choose from the admissible set of doses based on inverse of the sample size
	## "ss-random" - Random sample size strategy. Choose from the admissible set of doses based on randomisation, weighted by the inverse of the sample size
## admis - Choice of admissible doses around the MTC. Choices are
	## "adjacent" - All doses that are adjacent the MTC are admissible (even corner doses)
	## "closest" - All doses that are closest to the MTC
## constraint - dose-skipping constraint (defaults to no constraint). Options are
	## "neighbouring" - any dose combination can be chosen up to one dose level above OR BELOW current drug A and drug B levels (i.e. neighbouring current dose)
	## "no.dose.skip" - any dose combination can be chosen up to one dose level above ANY previously experimented drug A and drug B levels
## epsilon - should there be a safety constraint imposed if the posterior probability(dose > MTC)>epsilon. This averages over posterior distribution of MTC contours
	## Defaults to NULL
	## Any number is the posterior tail probability for which any dose that exceeds this is not experimented on.
## test - Testing mode. Options are
	## "none" - no testing (use for simulations) [the default]
	## "nodlt" - Every patient is DLT free
	## "alldlt"  - Every patient suffers a DLT
## data (Optional) A named data frame giving information about dose and toxicity from previously recruited patients. If missing, then it is assumed that no data have thus far been collected. Contains the following variables:
    ##patient Recruited patient numbers, 1,...,n
    ##doseA Dose levels of Drug A 
    ##doseB Dose levels of Drug B
    ##tox An indicator variable for each patient (1=toxicity, 0=no toxicity)
## a - Matrix of values for prior Beta(a,b) distributions [Alternative to specifying prior.med and prior.ss]
## b - Matrix of values for prior Beta(a,b) distributions [Alternative to specifying prior.med and prior.ss]
## alternate - deescalate if above the MTC, escalate if below - defaults to FALSE
## uppertox.constraint - should an upper toxicity safety constraint be imposed?
	## Defaults to NULL
	## Any number is the upper toxicity level from which no dose should be experimented or recommended. 
	## N.B. This constraint is NOT documented in the Mander and Sweeting 2015 paper and its performance is often substandard to using the safety constraint that weights over all possible MTCs, due to rigidity. It should therefore be used with caution. 
## stop - an alternative stopping rule, whereby we stop if posterior prob. of being > TTL at lowest dose is > stop
	## Defaults to null
	## Any number >0 and <1
	## N.B. Not documented in Mander and Sweeting 2015.


pipe.design<-function(N=dim(data)[1]+1,S=1,c,theta,pi=NULL,prior.med,prior.ss,strategy,admis,constraint=NULL,epsilon=NULL,mode="sim",data=matrix(nrow=0,ncol=0),a=NULL,b=NULL,alternate=FALSE,uppertox.constraint=NULL,stop=NULL){

	# Dimensions of two-agent design
	I=dim(prior.med)[1]
	J=dim(prior.med)[2]
	
	# Do some checks of inputs
	if(!(strategy %in% c("ss","ss-random"))) stop("strategy must be one of `ss', `ss-random' or `p'")
	if(!(admis %in% c("adjacent","closest"))) stop("admis must be one of `adjacent', `closest' or `both'")
	if(!is.null(constraint) & !(constraint %in% c("neighbouring","no.dose.skip"))) stop("constraint must be one of `neighbouring' or `no.dose.skip'")
	if(!(mode %in% c("sim","nodlt","alldlt"))) stop("mode must be one of `sim', `nodlt', or `alldlt'")
	if(mode=="sim" & is.null(pi)){
		cat("`pi' must be specified to conduct a simulation study, only next recommended dose will be given \n")
		if(!is.null(data)){
			N=dim(data)[1]+1
		} else {
			N=1
		}
	}
	if(!is.null(uppertox.constraint)){
		if((theta>uppertox.constraint | uppertox.constraint>1)) stop("uppertox.constraint must be a number between theta and 1")
	}
	if(!is.null(epsilon)){
		if(epsilon<0 | epsilon>1) stop("epsilon must be a number between 0 and 1")
	}
	if(ceiling(N/c)!=floor(N/c)) stop("Total sample size `N' must be divisible by cohort size `c'")
	if (nrow(data)>0) {
        if (any(!(c("patient", "doseA", "doseB","tox") %in% names(data)))) 
            stop("data must have variables named 'patient', 'doseA', 'doseB' and 'tox'")
        data <- data[order(data$patient), ]
        if (any(data$patient != 1:dim(data)[1])) 
            stop("'patient' variable in data must be an ascending vector of positive integers")
        if (any(!(data$tox %in% c(0, 1)))) 
            stop("'tox' variable in data must be a vector of zeros (no toxicity) and ones (toxicity)")
        if (any(!(data$doseA %in% 1:I))) 
            stop(paste("'doseA' variable in data must contain the dose levels (1 to ",I,")", sep = ""))
        if (any(!(data$doseB %in% 1:J))) 
            stop(paste("'doseB' variable in data must contain the dose levels (1 to ",J,")", sep = ""))
	}            


	## No. of dose combinations
	k=I*J

	## Monotonic matrices
	matrices<-monotonic.matrices(I,J)

	# Number of doses being recommended for next cohort 
	# Currently PIPE only recommends 1 dose-combination for the next cohort
	doses<-1

	# Set up list of No. DLTs and Number Patients, recommended doses and number RPII doses per simulation
	r.sim<-n.sim<-list()
	rec.i.sim<-rec.j.sim<-array(NA,dim=c(S,N/c,doses))
	rec<-matrix(0,ncol=J,nrow=I)
	n.rpII<-vector()

	## Find a and b parameters from beta distribution with median = prior.med and prior strength given by prior.ss
	if(is.null(a) & is.null(b)){
		prior <- beta.med(prior.med,prior.ss)
		a<-prior$a
		b<-prior$b
	}

	# Set up more lists
	mat.list=uppermat.list=uppermat2.list=n.list=r.list=list()
	for(s in 1:S){ # Loop over simulations
		# Initialise No DLTs and No. patients for each dose combination
		r=matrix(0,nrow=I,ncol=J)
		n=matrix(0,nrow=I,ncol=J)
				
		## Prior probability that each dose is less than or equal to theta
		p<-pbeta(theta,a,b)
		# If uppertox.constraint is specified find probability that each dose is less than or equal to uppertox
		if(!is.null(uppertox.constraint)){
			pconstraint<-pbeta(uppertox.constraint,a,b)
		} else {
			pconstraint<-NULL
		}
		# Initialise recommended doses for dimension i and j
		rec.i=rec.j=matrix(nrow=0,ncol=doses)

		## Where does the first cohort get dosed?
		create<-mtc.create(matrices,p,constraint,pconstraint,epsilon,admis,rec.i,rec.j,n)
		mat<-create$mat			
		if(strategy=="ss" | strategy=="ss-random"){
			## CHOOSE NEXT DOSE AS ONE LEAST EXPERIMENTED ON
			pi.theta<-1/(a+b)
		} 
		# Using admissible doses, and prior sample sizes choose the next dose
		nxt<-mtc(create$dominant,create$admissible,strategy,rec.i,rec.j,pi.theta,mat,p,alternate)
		
		# Store results if only one trial run
		if(S==1){
			# This is most likely monotonic contour
			mat.list[[1]]=mat
			# This is monotonic matrix corresponding to upper toxicity constrain
			uppermat.list[[1]]=create$matupper
			# This is monotonic matrix corresponding to weighted Posterior p(>MTC) for each dose combination
			uppermat2.list[[1]]=create$matupper2
			n.list[[1]]=n
			r.list[[1]]=r
		}
	
		rec.i=nxt$rec.i
		rec.j=nxt$rec.j

		for(m in 1:(N/c)){ # Loop over cohorts
			if(doses==1){		
				if(!is.null(data) & nrow(data)>=m*c){ # If data is already specified then add it	
					if(length(unique(data$doseA[(m*c-c+1):(m*c)]))>1 | length(unique(data$doseB[(m*c-c+1):(m*c)]))>1){
						stop("Data given does not have all patients in the same cohort on the same dose combination")				
					} 
					for(pt in (m*c-c+1):(m*c)){
						r[data$doseA[pt],data$doseB[pt]]<-r[data$doseA[pt],data$doseB[pt]]+data$tox[pt]
						n[data$doseA[pt],data$doseB[pt]]<-n[data$doseA[pt],data$doseB[pt]]+1
					}
					rec.i[m,1]<-data$doseA[m*c]
					rec.j[m,1]<-data$doseB[m*c]
				} else if(mode=="alldlt" | mode=="nodlt") { # If nodlt or alldlt specified, add 0 or c respectively, else generate from binomial with true p(DLT)=pi
					r[rec.i[m,1],rec.j[m,1]]<-r[rec.i[m,1],rec.j[m,1]]+ifelse(mode=="nodlt",0,c)
					n[rec.i[m,1],rec.j[m,1]]<-n[rec.i[m,1],rec.j[m,1]]+c
				} else if(!is.null(pi)) {
					r[rec.i[m,1],rec.j[m,1]]<-r[rec.i[m,1],rec.j[m,1]]+rbinom(1,c,pi[rec.i[m,1],rec.j[m,1]])
					n[rec.i[m,1],rec.j[m,1]]<-n[rec.i[m,1],rec.j[m,1]]+c
				} else {
					break
				}
			} else { # This section is not currently used and corresponds to if more than one dose recommended after each cohort
				r[rec.i[m,1],rec.j[m,1]]<-r[rec.i[m,1],rec.j[m,1]]+ifelse(mode=="nodlt",0,ifelse(mode=="alldlt",c/2,rbinom(1,c/2,pi[rec.i[m,1],rec.j[m,1]])))
				n[rec.i[m,1],rec.j[m,1]]<-n[rec.i[m,1],rec.j[m,1]]+c/2
				r[rec.i[m,2],rec.j[m,2]]<-r[rec.i[m,2],rec.j[m,2]]+ifelse(mode=="nodlt",0,ifelse(mode=="alldlt",c/2,rbinom(1,c/2,pi[rec.i[m,2],rec.j[m,2]])))
				n[rec.i[m,2],rec.j[m,2]]<-n[rec.i[m,2],rec.j[m,2]]+c/2
			}
		
			# Calculate new posterior probailities
			p<-pbeta(theta,a+r,b+n-r)
			if(!is.null(uppertox.constraint)){
				pconstraint<-pbeta(uppertox.constraint,a+r,b+n-r)
			} else {
				pconstraint<-NULL
			}
			create<-mtc.create(matrices,p,constraint,pconstraint,epsilon,admis,rec.i,rec.j,n)
			mat<-create$mat		
				if(strategy=="ss" | strategy=="ss-random"){
				## CHOOSE NEXT DOSE AS ONE LEAST EXPERIMENTED ON
				## OR WEIGHTED RANDOMISATION BASED ON INVERSE SAMPLE SIZE (INCLUDING PRIOR SS)
				pi.theta<-1/(a+b+n)
			} 
			
			if(S==1){ # For a single trial (i.e. no simulation) store information after each cohort
				mat.list[[m+1]]=mat
				uppermat.list[[m+1]]=create$matupper
				uppermat2.list[[m+1]]=create$matupper2
				n.list[[m+1]]=n
				r.list[[m+1]]=r
			}
			
			## IF NO DOSES ARE ADMISSIBLE THEN STOP THE TRIAL FOR SAFETY
			if(all(!create$admissible)){
				rec.i<-rbind(rec.i,matrix(0,nrow=N/c+1-m,ncol=doses))
				rec.j<-rbind(rec.j,matrix(0,nrow=N/c+1-m,ncol=doses))
				break
			}
			if(!is.null(stop)){
				## If lowest dose combination has posterior probability of being greater than the TTL of > stop then stop the trial
				if(1-p[1,1]>stop){
					rec.i<-rbind(rec.i,matrix(0,nrow=N/c+1-m,ncol=doses))
					rec.j<-rbind(rec.j,matrix(0,nrow=N/c+1-m,ncol=doses))
					break
				}
			}
		
			# Find next dose
			nxt<-mtc(create$dominant,create$admissible,strategy,rec.i,rec.j,pi.theta,mat,p,alternate)
			rec.i=nxt$rec.i
			rec.j=nxt$rec.j
		}
		r.sim[[s]]<-r
		n.sim[[s]]<-n
		## Recommendations made for each cohort throughout the trial
		rec.i.sim[s,,]<-rec.i[-(m+1),]
		rec.j.sim[s,,]<-rec.j[-(m+1),]
		
		## Recommended PII dose combinations
		## THESE ARE DOSE COMBINATIONS THAT HAVE BEEN EXPERIMENTED ON, ARE CLOSEST TO ESTIMATED MTC_THETA
		## AND ARE LOWER THAN UPPER CONSTRAINT CONTOUR, OR p(dose>MTC)< epsilon
		## CAN ONLY RECOMMEND PII DOSES IF TRIAL IS NOT STOPPED EARLY
		if(any(create$admissible)){
			create.rpII<-mtc.create(matrices,p,constraint="none",pconstraint=pconstraint,epsilon,admis="closest",rec.i,rec.j,n)
			rpIIs<-create.rpII$dominant & create.rpII$mat==0 & n!=0 & create$matupper==0 & create$matupper2==0 
			rpII.i<-row(mat)[rpIIs]
			rpII.j<-col(mat)[rpIIs]
			for(i in 1:length(rpII.i)){
				rec[rpII.i[i],rpII.j[i]]<-rec[rpII.i[i],rpII.j[i]]+1
			}
			n.rpII[s]<-length(rpII.i)
		} else {
		## If trial has stopped early
			n.rpII[s]<-0
		}
		cat(s,"\n")
	}
	exp<-Reduce('+',n.sim)/sum(Reduce('+',n.sim))
	no.not.treated<-N*S-sum(Reduce('+',n.sim))
	dlts<-sapply(1:S,function(l){sum(r.sim[[l]])/sum(n.sim[[l]])})
	results<-list(r.sim=r.sim,n.sim=n.sim,rec.i.sim=rec.i.sim,rec.j.sim=rec.j.sim,exp=exp,rec=rec/sum(rec),dlts=dlts,mat.list=mat.list,uppermat.list=uppermat.list,uppermat2.list=uppermat2.list,r.list=r.list,n.list=n.list,n.rpII=n.rpII,no.not.treated=no.not.treated,pi=pi,theta=theta)
	if(S>1){
		class(results)<-"pipe.sim"
	} else {
		class(results)<-"pipe"
	}
	return(results)
}


## Function that returns all monotonic matrices of dimension IxJ
monotonic.matrices<-function(I,J){
	comb.col<-combinations(2, J, c(0,1), repeats.allowed=TRUE)
	n.com1<-dim(comb.col)[1]
	comb.row<-combinations(n.com1,I,repeats.allowed=TRUE)
	n.com2<-dim(comb.row)[1]
	matrices<-sapply(1:n.com2,function(i){comb.col[comb.row[i,],]},simplify=F)
	return(matrices)
}

## Obtain doses closest to MTC
closest<-function(mat){
	I<-nrow(mat)
	J<-ncol(mat)
	dominantu<-mat==1 & rbind(0,mat[-I,]) %in% c(0,2) & cbind(0,mat[,-J]) %in% c(0,2)
	dominantl<-mat==0 & rbind(mat[-1,],1) %in% c(1,2) & cbind(mat[,-1],1) %in% c(1,2)
	dominant<-dominantl | dominantu
	dominant
}


## Create admissible dose matrix
mtc.create<-function(matrices,p,constraint,pconstraint,epsilon,admis,rec.i,rec.j,n){

	m<-dim(rec.i)[1]

	## Assess all possible MTC contours to find most likely
	mtc.num<-unlist(lapply(matrices,function(l){prod((1-p)[l==1])*prod(p[l==0])}))
	mtc.lik<-mtc.num/sum(mtc.num)
	mtc.mode<-which.max(mtc.lik)
	mat<-matrices[[mtc.mode]]
	
	I<-dim(mat)[1]
	J<-dim(mat)[2]

	## Find upper toxicity constraint contour if uppertox.constraint is used
	if(!is.null(pconstraint)){
		upper.lik<-which.max(unlist(lapply(matrices,function(l){prod((1-pconstraint)[l==1])*prod(pconstraint[l==0])})))
		matupper<-matrices[[upper.lik]]
	} else {
		matupper<-matrix(0,nrow=I,ncol=J)
	}

	## Find doses that do not satisfy constraint if weightedMTC.constraint is used
	if(!is.null(epsilon)){
		## Posterior p(>MTC) for each dose combination
		mat.lik<-sapply(1:length(mtc.lik),function(l){matrices[[l]]*mtc.lik[[l]]},simplify=F)
		weight.pMTC<-Reduce('+',mat.lik)
	} else {
		weight.pMTC<-NULL
	}
	if(!is.null(epsilon)){
		matupper2<-weight.pMTC>=epsilon
	} else {
		matupper2<-matrix(0,nrow=I,ncol=J)
	}
	
	if(constraint=="neighbouring"){
		## IF NEIGHBOURING CONSTRAINT AND MORE THAN ONE DOSE RECOMMENDATION PER COHORT THEN USE UNION OF BOTH ADMISSIBLE REGIONS
		if(dim(rec.i)[2]>1){
			admissible1<- row(mat)<=max(rec.i[m,1],0)+1 & col(mat)<=max(rec.j[m,1],0)+1 & row(mat)>=max(rec.i[m,1],0)-1 & col(mat)>=max(rec.j[m,1],0)-1
			admissible2<- row(mat)<=max(rec.i[m,2],0)+1 & col(mat)<=max(rec.j[m,2],0)+1 & row(mat)>=max(rec.i[m,2],0)-1 & col(mat)>=max(rec.j[m,2],0)-1
			admissible<- admissible1 | admissible2
		} else {
			admissible<- row(mat)<=max(rec.i[m,1],0)+1 & col(mat)<=max(rec.j[m,1],0)+1 & row(mat)>=max(rec.i[m,1],0)-1 & col(mat)>=max(rec.j[m,1],0)-1
		}
	} else if(constraint=="no.dose.skip"){
		## IF NO.DOSE.SKIP CONSTRAINT AND MORE THAN ONE DOSE RECOMMENDATION PER COHORT THEN USE UNION OF BOTH ADMISSIBLE REGIONS
		if(dim(rec.i)[2]>1){
			admissible1<- row(mat)<=max(rec.i[,1],0)+1 & col(mat)<=max(rec.j[,1],0)+1
			admissible2<- row(mat)<=max(rec.i[,2],0)+1 & col(mat)<=max(rec.j[,2],0)+1
			admissible<- admissible1 | admissible2
		} else {
			admissible<- row(mat)<=max(rec.i[,1],0)+1 & col(mat)<=max(rec.j[,1],0)+1 
		}
	} else {
		## IF NO NEIGHBOURING CONSTRAINT THEN ALL DOSE COMBINATIONS ARE ADMISSIBLE	
		admissible<- matrix(TRUE,nrow=I,ncol=J)		
	}
	if(!is.null(pconstraint) | !is.null(epsilon)){
		admissible<-admissible & matupper==0 & matupper2==0
		## IF THERE ARE NO ADMISSIBLE DOSES LEFT (THAT IS ALL NEIGHBOURING DOSES ARE NOW UNSAFE)
		## CHOOSE CLOSEST DOSE THAT IS SAFE (TO FIRST COHORT DOSE)
		if(all(!admissible)){
			test<-abs(rec.i[m,1]-row(mat))+abs(rec.j[m,1]-col(mat))
			admissible<- test==min(c(test[matupper==0 & matupper2==0],-Inf)) & matupper==0 & matupper2==0
		}
	}
	
	separate<-FALSE
	if(admis=="adjacent"){
  		## ANY DOSE COMBINATION ADJACENT TO THE MTC IS ADMISSIBLE
		if(I<2 | J<2) stop("Admissible doses can only be calculated when both drugs have more than one level")
		admat<-mat
		dominantu<-admat==1 & (rbind(0,admat[-I,])==0 | cbind(0,admat[,-J])==0 | rbind(0,cbind(0,admat[,-J])[-I,])==0)
		dominantl<-admat==0 & (rbind(admat[-1,],1)==1 | cbind(admat[,-1],1)==1 | rbind(cbind(admat[,-1],1)[-1,],1)==1)
		dominant<-dominantl | dominantu
		## If dominant and admissible regions are separate choose the "closest" dose in the admissible region
		if(!any(dominant & admissible)){
			separate<-TRUE
		}
	} 
	if(admis=="closest" | separate==TRUE){
		## ONLY DOSE COMBINATIONS CLOSEST TO THE MTC ARE ADMISSIBLE
		admat <- mat
		## SET ALL DOSES OUTSIDE ADMISSIBLE RANGE TO 2 AND ALLOW ANY TOUCHING CONSTRAINT TO BE DOMINANT (IF SATISFY OTHER CRITERIA)
		admat[admissible==FALSE]=2
							
		dominant<-closest(admat)
	}
	return(list(dominant=dominant,admissible=admissible,mat=mat,matupper=matupper,matupper2=matupper2,weight.pMTC=weight.pMTC))
}

## Choose next dose from admissible doses that use either smallest sample size or weighted randomisation of sample size
mtc<-function(dominant,admissible,strategy,rec.i,rec.j,pi.theta,mat,p,alternate){

	m<-dim(rec.i)[1]
	I<-dim(pi.theta)[1]
	J<-dim(pi.theta)[2]
	k<-I*J

	## Strategy "ss": Select the dominant dose with smallest sample size
	## If there are more than two dose combinations that are equal choose from them at random
	if(strategy=="ss"){
		if(alternate==T){
			## If more than one admissible & dominant dose then go below MTC if last dose is above, or vice-versa
			if(sum(dominant & admissible)>1){
				if(mat[rec.i[[m]],rec.j[[m]]]==1){
					if(any(mat[dominant & admissible]==0)) dominant[mat==1]<-FALSE
				} else{
					if(any(mat[dominant & admissible]==1)) dominant[mat==0]<-FALSE
				}		
			}
		}
		## If still more than one dose comb. then choose one with smallest ss then choose at random
		test<- pi.theta==max(pi.theta[dominant & admissible]) & dominant & admissible
		chosen=ifelse(sum(test)>1,sample(sum(test),1),1)
		rec.i<-rbind(rec.i,row(pi.theta)[test][chosen])
		rec.j<-rbind(rec.j,col(pi.theta)[test][chosen])
	} else if(strategy=="ss-random"  | strategy=="weighted_mtc"){
		## Strategy "ss-random": Select the dominant dose with probability weighted by inverse sample size or weighted by the (weighted) probability of being a closest dose to the true MTC
		pi.theta[!(dominant & admissible)]=0
		chosen=sample(k,1,prob=pi.theta)
		rec.i<-rbind(rec.i,row(pi.theta)[chosen])
		rec.j<-rbind(rec.j,col(pi.theta)[chosen])
	} else if(strategy=="p"){
		## Strategy "p": Select the dominant dose with posterior probability of being less than p closest to 0.5 (i.e. most uncertain)
		p[!(dominant & admissible)]=2
		test<- abs(p-0.5)==min(abs(p-0.5)) & dominant & admissible
		## If more than one dose comb. with same ss then choose at random
		chosen=ifelse(sum(test)>1,sample(sum(test),1),1)
		rec.i<-rbind(rec.i,row(pi.theta)[test][chosen])
		rec.j<-rbind(rec.j,col(pi.theta)[test][chosen])
	} 
	return(list(rec.i=rec.i,rec.j=rec.j))	
}




## Obtain a and b parameters for beta prior from median and sample size using numerical optimisation
beta.med<-function(prior.med,prior.ss){
	opt.function<-function(par,prior.med,prior.ss){
		a<-pmin(exp(par),prior.ss)
		b<-prior.ss-a
		suppressWarnings(sum(abs(pbeta(prior.med,a,b)-0.5)))
	}
	init.par<-log(prior.ss/2)
	a.med<-exp(optim(init.par,opt.function,prior.med=prior.med,prior.ss=prior.ss,method="BFGS")$par)
	b.med<-prior.ss-a.med
	return(list(a=a.med,b=b.med))
}


## Function to plot dose-escalation steps for a PIPE design
## Arguments:
## x: a PIPE object as obtained from pipe
## type: What data should be shown on the figure. Options are
## 			"b" (default): Show both the number of DLTs (numerator) and number of patients recruited (denominator) at each dose combination
## 			"r": Only show number of DLTs (numerator) at each dose combination
## 			"n": Only show number of patients recruited (denominator) at each dose combination
## pi: True p(DLT) matrix
## theta: Target toxicity level
## epsilon.line: Plot safety constraint line formed by a weighted average of all possible MTCs
## uppertox.constraint.line: Plot safety constraint line based on most likely MTC at upper toxicity constraint
plot.pipe<-function(x,type="b",pi=x$pi,theta=x$theta,epsilon.line=TRUE,uppertox.constraint.line=FALSE,...){
	mat<-x$mat.list
	c<-length(mat)
	I<-dim(mat[[1]])[1]
	J<-dim(mat[[1]])[2]
	cohort.size=sum(x$n.sim[[1]])/(c-1)
	xlevels=rep(1:(I+1)-0.5,c)
	ylevels=c(sapply(1:c,function(k){c(apply(mat[[k]],1,function(i){min(which(i==1),J+1)-0.5}),0.5)}))
	cohort=rep(1:c,each=(I+1))
	if(!is.null(pi)){
		mat.true<-pi>theta
		x.true<-1:(I+1)-0.5
		y.true<-c(apply(mat.true,1,function(i){min(which(i==1),J+1)-0.5}),0.5)
	}
			
	df<-data.frame(x=xlevels,y=ylevels,cohort=cohort)
	ncol<-2 ## Minimum number of coloured tiles to be shown in plot
	## Data to be shown
	if(type=="n"){
		df2<-data.frame(x=rep(rep(1:I,J),c-1),y=rep(rep(1:J,each=I),c-1),z=factor(unlist(x$n.list[-1]),levels=unique(sort(unlist(x$n.list[-1]),decreasing=T))),cohort=rep(1:(c-1),each=I*J))
		ncol<-length(levels(df2$z))
	} else if(type=="r"){
		df2b<-data.frame(x=rep(rep(1:I,J),c-1),y=rep(rep(1:J,each=I),c-1),z=unlist(x$r.list[-1]),cohort=rep(1:(c-1),each=I*J))
		df2b<-df2b[df2b$z!=0,]
		df2b$z<-factor(df2b$z,levels=unique(sort(df2b$z,decreasing=T)))
	} else if(type=="b"){
		df2<-data.frame(x=rep(rep(1:I,J),c-1),y=rep(rep(1:J,each=I),c-1),z=factor(unlist(x$n.list[-1]),levels=unique(sort(unlist(x$n.list[-1]),decreasing=T))),cohort=rep(1:(c-1),each=I*J))
		ncol<-length(levels(df2$z))
		df2b<-data.frame(x=rep(rep(1:I,J),c-1),y=rep(rep(1:J,each=I),c-1),z=unlist(x$r.list[-1]),cohort=rep(1:(c-1),each=I*J))
		df2b<-df2b[df2b$z!=0,]
		df2b$z<-factor(df2b$z,levels=unique(sort(df2b$z,decreasing=T)))
	}
	if(!is.null(pi)){
		df3<-data.frame(x=x.true,y=y.true,cohort=rep(c,I+1))
	}
	df4<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=factor(c(cohort.size*as.numeric(x$rec!=0)),levels=c(cohort.size,0)),cohort=rep(c,I*J))
	v1<-ggplot()+
		geom_step(aes(x=x,y=y),data=df,size=2)+facet_wrap(~cohort)+
		xlab("Drug A level")+ylab("Drug B level")
	if(!is.null(pi)){
		v1<-v1+geom_step(aes(x=x,y=y),data=df3,size=1.5,colour="green",linetype=4)
	}
	if(type=="n" | type=="b"){
		v1<-v1+geom_tile(aes(x=x,y=y,fill = z),alpha=0.5,data=df2)
	}
	if(type=="r" | type=="b"){
		if(dim(df2b)[1]>0){
			v1<-v1+geom_point(aes(x=x,y=y,shape = z),size=2,data=df2b)+scale_shape(name="Number DLTs")	
		}
	}
	v1<-v1+geom_tile(aes(x=x,y=y,fill = z),alpha=0.5,data=df4)+scale_fill_manual(name="Number pts.",values=c(rainbow(ncol-1),"#FFFFFFFF"))
	if(length(x$uppermat.list)>0 & uppertox.constraint.line){
		uppermat<-x$uppermat.list
		x.high<-rep(1:(I+1)-0.5,c)
		y.high<-c(sapply(1:c,function(k){c(apply(uppermat[[k]],1,function(i){min(which(i==1),J+1)-0.5}),0.5)}))
		df5<-data.frame(x=x.high,y=y.high,cohort=cohort)
		v1<-v1+geom_step(aes(x=x,y=y),data=df5,size=1,colour="red")
	}
	if(length(x$uppermat2.list)>0 & epsilon.line){
		uppermat2<-x$uppermat2.list
		x.high<-rep(1:(I+1)-0.5,c)
		y.high<-c(sapply(1:c,function(k){c(apply(uppermat2[[k]],1,function(i){min(which(i==1),J+1)-0.5}),0.5)}))
		df5<-data.frame(x=x.high,y=y.high,cohort=cohort)
		v1<-v1+geom_step(aes(x=x,y=y),data=df5,size=1,colour="red4",linetype=4)
	}
	print(v1)
}

## Function to plot dose-escalation steps for a PIPE design
## Arguments:
## x: a PIPE object as obtained from pipe
## pi: True p(DLT) matrix
## theta: Target toxicity level
## plot: What operating characteristics should be plotted? Options are:
##		"exp": Experimentation proportions heat map
##		"rec": Recommendation proportions heat map
plot.pipe.sim<-function(x,pi=x$pi,theta=x$theta,plot="both",...){
	exp<-x$exp
	rec<-x$rec
	I<-dim(x$n.sim[[1]])[1]
	J<-dim(x$n.sim[[1]])[2]

	mat.true<-pi>theta
	x.true<-1:(I+1)-0.5
	y.true<-c(apply(mat.true,1,function(i){min(which(i==1),J+1)-0.5}),0.5)
	
	s<-length(x$n.sim)
	
	if(plot=="exp"){
		# Experimentation proportions plot
		df<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=100*c(exp))
		df2<-data.frame(x=x.true,y=y.true)
		v1<-ggplot()+geom_tile(aes(x=x,y=y,fill = z),data=df)+scale_fill_gradient(name="Experimentation percentages",low="white",high="red")+
			geom_step(aes(x=x,y=y),data=df2,size=1.5,colour="green",linetype=4)+
			xlab("Drug A level")+ylab("Drug B level")
		print(v1)
	}

	if(plot=="rec"){
		# Recommendation proportions
		df<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=100*c(rec))
		df2<-data.frame(x=x.true,y=y.true)
		v1<-ggplot()+geom_tile(aes(x=x,y=y,fill = z),data=df)+scale_fill_gradient(name="Recommendation percentages",low="white",high="red")+
			geom_step(aes(x=x,y=y),data=df2,size=1.5,colour="green",linetype=4)+
			xlab("Drug A level")+ylab("Drug B level")
		print(v1)
	}

	if(plot=="both"){
		# Experimentation and Recommendation proportions plot
		df.exp<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=100*c(exp),type="Experimentation")
		df.rec<-data.frame(x=rep(1:I,J),y=rep(1:J,each=I),z=100*c(rec),type="Recommendation")
		df<-rbind(df.exp,df.rec)
		df2<-data.frame(x=x.true,y=y.true)
		v1<-ggplot()+geom_tile(aes(x=x,y=y,fill = z),data=df)+facet_grid(~type)+scale_fill_gradient(name="Percent",low="white",high="red")+
			geom_step(aes(x=x,y=y),data=df2,size=1.5,colour="green",linetype=4)+
			xlab("Drug A level")+ylab("Drug B level")
		print(v1)
	}
}


print.pipe<-function(x,...){
	I=dim(x$r.sim[[1]])[1]
	J=dim(x$r.sim[[1]])[2]
	n<-x$n.sim[[1]]
	r<-x$r.sim[[1]]
	mat<-x$mat.list[[length(x$mat.list)]]
	matupper<-x$uppermat.list[[length(x$uppermat.list)]]
	matupper2<-x$uppermat2.list[[length(x$uppermat2.list)]]
	tab1<-t(n)[J:1,]
	rownames(tab1)<-paste("Level",J:1)
	colnames(tab1)<-paste("Level",1:I)
	names(dimnames(tab1))<-c("Drug B","Drug A")
	tab2<-t(r)[J:1,]
	rownames(tab2)<-paste("Level",J:1)
	colnames(tab2)<-paste("Level",1:I)
	names(dimnames(tab2))<-c("Drug B","Drug A")
	tab3<-t(mat)[J:1,]
	rownames(tab3)<-paste("Level",J:1)
	colnames(tab3)<-paste("Level",1:I)
	names(dimnames(tab3))<-c("Drug B","Drug A")
	tab4<-t((matupper | matupper2)*1)[J:1,]
	rownames(tab4)<-paste("Level",J:1)
	colnames(tab4)<-paste("Level",1:I)
	names(dimnames(tab4))<-c("Drug B","Drug A")
	cat("\n Number of patients dosed:\n")
	print(tab1)
	cat("\n Toxicities observed:\n")
	print(tab2)
	if(length(x$r.list)==length(x$rec.i.sim)){
		cat("\n Next recommended dose level: \n Dose A: ",x$rec.i.sim[1,length(x$rec.i.sim),1],"\n Dose B: ",x$rec.j.sim[1,length(x$rec.j.sim),1],"\n")
	}
	cat("\n MTC:\n")
	print(tab3)
	cat("\n Upper toxicity constraint (1 indicates doses not allowed):\n")
	print(tab4)
}

###### Function to print experimentation and recommendation percentages from a simulated PIPE design
## Arguments:
## x: a PIPE object as obtained from pipe
## pi: True p(DLT) matrix
## cut.points: cut points on the true DLT range to present the operating characteristics
## digits: Number of decimal places to be used for reporting
## print: (defaul=TRUE). Should the output be printed?
print.pipe.sim<-function(x,pi=x$pi,cut.points=c(0,15,25,35,45.1,100)/100,digits=1,print=TRUE,...){
	exp<-x$exp
	rec<-x$rec

	cuts<-cut(pi,cut.points,right=F)
	# Experimentation proportions table
	exp.table<-sapply(levels(cuts),function(i){sum(exp[cuts==i])})

	# Recommendation proportions table
	rec.table<-sapply(levels(cuts),function(i){sum(rec[cuts==i])})

	if(print){
		cat("\n Experimentation percentages by true toxicity: \n")
		print(round(100*exp.table,digits))

		cat("\n Recommendation percentages by true toxicity: \n")
		print(round(100*rec.table,digits))		
	}
	return(list(exp.table=exp.table,rec.table=rec.table))
}


