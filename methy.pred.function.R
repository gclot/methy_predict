
require("e1071")

methy.prediction=function(data,which.predictor="all",impute.missings=FALSE,
	path="/home/guillem/analisis/marti/v2/analisis/prepare easy code/upload to github/methy.pred.RData"){

	# data: cpgs in columns and samples in rows
	# which.predictor: one of c("main","ALL","CLL","DLBCL","MCL")
	# path: path to methy.pred.RData (the dfault path is to a github repository)

	# if which predictor have length two or more, use only the first
	if (length(which.predictor)>=2){
		which.predictor=which.predictor[1]
		warning(paste0("Only first value of which.predictor object (",which.predictor,") is used"))
	}

	# check if which.predictor is correctly specified
	if (!which.predictor %in% c("all","main","ALL","CLL","DLBCL","MCL")) stop(paste0("which.predictor value (",which.predictor,") is not valid"))

	# load predictor and cpgs
	load(path)

	# subset the cpg list according to which.predictor
	if (!which.predictor %in% "all"){
		list.cpgs=list.cpgs[list.cpgs$predictor %in% which.predictor,]
	}

	# check for missing cpgs in data object
	if (!all(list.cpgs$cpg %in% colnames(data))){
		w=list.cpgs$cpg[!list.cpgs$cpg %in% colnames(data)]
		stop(paste0("The following cpgs are missing in data object: ",paste(w,collapse=",")))	
	}

	# subset data to the relevant cpgs

	data=as.matrix(data[,list.cpgs$cpg,drop=F])

	# if impute.missings=FALSE, stop if missings are found

	if ((!impute.missings)&any(is.na(data))){
		w=which(is.na(data),T)
		w=paste(paste0(rownames(w)," in ",colnames(data)[w[,"col"]]),collapse="\n")
		stop(paste0("Missing values in data object and option impute.missings=FALSE. Found at:\n",w))
	}

	# check for missing values in cpgs of the CLL subtypes (not supported)

	if (which.predictor %in% c("all","CLL")){
		if (any(is.na(data[,list.cpgs[list.cpgs$predictor %in% "CLL","cpg"]]))){
			stop("Missing values in CLL subtype cpgs are not supported (cg03462096,cg17014214,cg11472422,cg00869668)")
		}
		
	}
		
	# impute missing values in other classes

	if (any(is.na(data))){

		w=which(is.na(data),T)
		w=paste(paste0(rownames(w)," in ",colnames(data)[w[,"col"]]),collapse="\n")

		orig=data

		to.impute=sweep(is.na(data),MARGIN=2,STATS=list.cpgs[,"impute"],FUN="*")
		to.impute[is.na(to.impute)]=0

		data[is.na(data)]=0
		data=data+to.impute

		warning(paste0("Missing values were imputed at:\n",w))
	}

	# calc signatures

	if (!which.predictor %in% "CLL"){
		temp=list.cpgs[!list.cpgs$predictor %in% "CLL",]

		signatures=by(data=temp[,c("cpg","sign")],INDICES=temp$class,FUN=function(x,meth){
			rowMeans(sweep(meth[,x$cpg,drop=F],MARGIN=2,STATS=x$sign,FUN=function(x,y) abs(y-x)))
		},meth=data)
		signatures=do.call("cbind",signatures)
		if (which.predictor %in% "all"){
			signatures=cbind(signatures,data[rownames(signatures),list.cpgs[list.cpgs$predictor %in% "CLL","cpg"]])
		}
	} else {
		temp=list.cpgs[list.cpgs$predictor %in% "CLL",]
		signatures=data[,temp$cpg]
	}

	# apply predictors
	
	output.list=list()

	if (which.predictor %in% c("all","main")){
		svm.pred=predict(fit.main,newdata=signatures[,unique(list.cpgs[list.cpgs$predictor %in% "main","class"]),drop=F],probability=T)
		res=data.frame(attributes(svm.pred)$probabilities,check.names=F)
		suggested.classes=apply(res>0.35,1,function(x,y){paste(y[x],collapse="|")},y=colnames(res))
		suggested.classes[apply(res,1,max)<0.5]="unclassified"
		colnames(res)=paste0("P(",colnames(res),")")
		res$svm.prediction=as.character(svm.pred)
		res$suggested.classes=suggested.classes
		output.list$main.predictor=res
	}
	if (which.predictor %in% c("all","ALL")){
		svm.pred=predict(fit.ALL,newdata=signatures[,unique(list.cpgs[list.cpgs$predictor %in% "ALL","class"]),drop=F],probability=T)
		res=data.frame(attributes(svm.pred)$probabilities,check.names=F)
		suggested.classes=apply(res>0.35,1,function(x,y){paste(y[x],collapse="|")},y=colnames(res))
		suggested.classes[apply(res,1,max)<0.5]="unclassified"
		colnames(res)=paste0("P(",colnames(res),")")
		res$svm.prediction=as.character(svm.pred)
		res$suggested.classes=suggested.classes
		output.list$ALL.predictor=res
	}
	if (which.predictor %in% c("all","CLL")){
		svm.pred=predict(fit.CLL,newdata=signatures[,unique(list.cpgs[list.cpgs$predictor %in% "CLL","cpg"]),drop=F],probability=T)
		res=data.frame(attributes(svm.pred)$probabilities,check.names=F)
		suggested.classes=apply(res>0.35,1,function(x,y){paste(y[x],collapse="|")},y=colnames(res))
		suggested.classes[apply(res,1,max)<0.5]="unclassified"
		colnames(res)=paste0("P(",colnames(res),")")
		res$svm.prediction=as.character(svm.pred)
		res$suggested.classes=suggested.classes
		output.list$CLL.predictor=res
	}
	if (which.predictor %in% c("all","DLBCL")){
		svm.pred=predict(fit.DLBCL,newdata=signatures[,unique(list.cpgs[list.cpgs$predictor %in% "DLBCL","class"]),drop=F],probability=T)
		res=data.frame(attributes(svm.pred)$probabilities,check.names=F)
		suggested.classes=apply(res>0.35,1,function(x,y){paste(y[x],collapse="|")},y=colnames(res))
		suggested.classes[apply(res,1,max)<0.65]="unclassified"
		colnames(res)=paste0("P(",colnames(res),")")
		res$svm.prediction=as.character(svm.pred)
		res$suggested.classes=suggested.classes
		output.list$DLBCL.predictor=res
	}
	if (which.predictor %in% c("all","MCL")){
		svm.pred=predict(fit.MCL,newdata=signatures[,unique(list.cpgs[list.cpgs$predictor %in% "MCL","class"]),drop=F],probability=T)
		res=data.frame(attributes(svm.pred)$probabilities,check.names=F)
		colnames(res)=substr(colnames(res),start=1,stop=2)
		suggested.classes=apply(res>0.35,1,function(x,y){paste(y[x],collapse="|")},y=colnames(res))
		suggested.classes[apply(res,1,max)<0.65]="unclassified"
		colnames(res)=paste0("P(",colnames(res),")")
		res$svm.prediction=substr(as.character(svm.pred),start=1,stop=2)
		res$suggested.classes=suggested.classes
		output.list$MCL.predictor=res
	}

	return(output.list)

}


