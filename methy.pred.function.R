
require("e1071")

methy.prediction=function(data,path="https://github.com/gclot/methy_predict/blob/master/methy.pred.RData",which.predictor="all"){

	# data: cpgs in columns and samples in rows

	# check if which predictor is correctly specified
	if (!which.predictor %in% c("all","main","ALL","CLL","DLBCL","MCL")) stop("Specify a valid predictor")

	# aixó comprobar si realment s'ho carrega, crec que si
	# warn if some of the global enviroment objects is going to be overwritten
	glob.env.obj=ls(envir = .GlobalEnv)
	if (any(c("fit.main","fit.ALL","fit.CLL","fit.MCL","fit.DLBCL","list.cpgs") %in% glob.env.obj)) {
		w=glob.env.obj[glob.env.obj %in% c("fit.main","fit.ALL","fit.CLL","fit.MCL","fit.DLBCL","list.cpgs")]
		print("The following objects have been overwritten:",paste(w,collapse=","))
	}

	# load predictor and cpgs
	load(path)

	# subset the cpg list according to which.predictor
	if (!which.predictor %in% "all"){
		list.cpgs=list.cpgs[list.cpgs$predictor %in% which.predictor,]
	}


	# check if any cpgs are missing in data
	if (!all(list.cpgs$cpg %in% colnames(data))){
		w=colnames(data)[colnames(data) %in% list.cpgs$cpg]
		stop(paste0("The following cpgs are missing in data object: ",paste(w,collapse=",")))	
	}

	# subset data to the relevant cpgs

	data=as.matrix(data[,list.cpgs$cpg,drop=F])

	# impute missings with a warning (impute value = ¿mean(class means)? or mean(class mean and mean of other classes?))

	if (any(is.na(data))){

		# pending work in here ...

		print("Found missing values in data object. A value have been imputed to missing data")
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
	{

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

	remove(fit.main,fit.ALL,fit.CLL,fit.MCL,fit.DLBCL,list.cpgs)

	return(output.list)

}


