# load libraries
# .libPaths("C:/Users/wchiu/Documents/R/win-library/3.3") # For testing
library("randomForest")
library("rcdk") # Note this requres JAVA that is that same build (e.g., 64-bit) as R

## Set up -- only do this once
load("ToxValRFModels.Rdata")
tv.residuals <- read.csv("Prediction_Residuals.csv",row.names=1)
row.names(tv.residuals) <- ToxValuesNames
dc <- get.desc.categories()
dnames<-c(get.desc.names(dc[1]),get.desc.names(dc[2]),get.desc.names(dc[3]),get.desc.names(dc[4]),get.desc.names(dc[5]))
dnames<-unique(dnames)
## Get file with reg tox values
tvdat <- read.csv("CTV_data_2016-xls.csv",as.is=TRUE)

## Get file with SMILES and MW
dat <- read.csv("EXAMPLE_CAS_SMILES_MW.csv",as.is=TRUE)

infile<-"EXAMPLE"
## Get descriptors from smiles - do this for each chemical
y.pred.df<-data.frame()
for (j in 1:100) {
  print(paste(j,dat$Name.Check[j],dat$CAS[j]))
  if (intvdat<-dat$CAS[j] %in% tvdat$CASRN) { # Use canonical SMILES from training dataset if same CASRN
    jmatch<-match(dat$CAS[j],tvdat$CASRN)
    newchem.smiles<-tvdat$CanonicalSMILES[jmatch]
  } else {
    newchem.smiles<-dat$SMILES[j]
  }
  if (newchem.smiles !="") {
    mw <- as.numeric(dat$AVERAGE_MASS[j])
    mol <- parse.smiles(newchem.smiles)[[1]]
    atomicnumbers <- unique(unlist(lapply(get.atoms(mol),get.atomic.number)))
    
    nocarbon <- !(6 %in% atomicnumbers)
    metals_etc <- c(5,13:14,21:33,39:52,57:71,89:103,72:84,104:112)
    hasmetals <- sum(metals_etc %in% atomicnumbers) > 0
    
    mol.desc<-cbind(data.frame(smiles=newchem.smiles,stringsAsFactors = FALSE),eval.desc(mol,dnames))
    
    ## Clean up descriptors
    mol.desc.clean <- mol.desc[,col.keep]
    frac.imputed <- sum(is.na(mol.desc.clean))/length(mol.desc.clean)
    
    mol.desc.clean.imputed <- mol.desc.clean
    
    # Impute NA with median of training set (if necessary)
    for (i in 1:ncol(mol.desc.clean.imputed)) {
      mol.desc.clean.imputed[is.na(mol.desc.clean.imputed[,i]),i]<-col.med[i]
    }
    
    ## Make predictions
    x.new<-mol.desc.clean.imputed[,-1]
    y.pred<-data.frame(tv=sub(".train","",ToxValuesNames),row.names=ToxValuesNames)
    y.pred$Name <- dat$Name.Check[j]
    y.pred$CAS <- dat$CAS[j]
    y.pred$SMILES <- newchem.smiles
    y.pred$MW <- mw
    y.pred$log10prediction<-0
    y.pred$log10lower95<-0
    y.pred$log10upper95<-0
    y.pred$prediction<-0
    y.pred$lower95<-0
    y.pred$upper95<-0
    y.pred$appl.domain<-0
    y.pred$Source<-"CTV"
    y.pred$Warnings<-""
    if (nocarbon | hasmetals | (frac.imputed > 0.1)) {
      if (nocarbon) y.pred$Warnings <- paste(y.pred$Warnings,"Inorganic (no carbons).")
      if (hasmetals) y.pred$Warnings <- paste(y.pred$Warnings,"Contains metals/metalloids.")
      if (frac.imputed > 0.1) y.pred$Warnings <- paste(y.pred$Warnings,">10% of molecular descriptors imputed.")
    }
    if (!nocarbon & !hasmetals) {
      for (tval in ToxValuesNames) {
        y.pred[tval,"log10prediction"] <- predict(tv.model[[tval]],x.new)
        y.pred[tval,"log10lower95"] <- y.pred[tval,"log10prediction"] + tv.residuals[tval,"CDK.ci.lb"]
        y.pred[tval,"log10upper95"] <- y.pred[tval,"log10prediction"] + tv.residuals[tval,"CDK.ci.ub"]
        if (tval == "OSF.train" | tval == "CPV.train" | tval == "IUR.train") {
          y.pred[tval,"prediction"] <- 10^y.pred[tval,"log10prediction"] / (mw*1000)
          y.pred[tval,"lower95"] <- 10^y.pred[tval,"log10lower95"] / (mw*1000)
          y.pred[tval,"upper95"] <- 10^y.pred[tval,"log10upper95"] / (mw*1000)
        } else {
          y.pred[tval,"prediction"] <- 10^(-y.pred[tval,"log10prediction"]) * (mw*1000)
          y.pred[tval,"lower95"] <- 10^(-y.pred[tval,"log10upper95"]) * (mw*1000)
          y.pred[tval,"upper95"] <- 10^(-y.pred[tval,"log10lower95"]) * (mw*1000)
        }
      }
      
      ## Check applicability domain Z-score
      
      for (tval in ToxValuesNames) {
        print(tval)
        ## Combine new chemical and tox value descriptors, and scale
        tv.desc <- tv.train[[tval]][,-(1:2)]
        all.desc <- rbind(x.new,tv.desc)
        all.desc.scale <- as.data.frame(scale(all.desc,center=TRUE,scale=TRUE))
        ## Write ".x" file for new chemical
        newchem_xfile<-paste(infile,"_newchem.x",sep="")
        x.new.scale <- all.desc.scale[1,]
        write.table(as.data.frame(t(dim(x.new.scale))),file=newchem_xfile,
                    row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
        write.table(as.data.frame(t(names(x.new.scale))),file=newchem_xfile,
                    row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
        write.table(cbind(newchem.smiles,x.new.scale),file=newchem_xfile,
                    row.names=TRUE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
        ## Write ".x" file for tox values
        tv.desc.scale <- all.desc.scale[-1,]
        tv_xfile<-paste(infile,"_",tval,".x",sep="")
        write.table(as.data.frame(t(dim(tv.desc.scale))),file=tv_xfile,
                    row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
        write.table(as.data.frame(t(names(tv.desc.scale))),file=tv_xfile,
                    row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
        write.table(cbind(tv.train[[tval]][,1],tv.desc.scale),file=tv_xfile,
                    row.names=TRUE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
        ## Calculate z-score
        gad_outfile <- paste(infile,"_newchem_by_",tval,".gad",sep="");
        system(paste("get_ad.exe ",tv_xfile," -4PRED=",newchem_xfile," -Z=3 -OUT=",gad_outfile,sep="")); 
        report.dat<-read.table(gad_outfile,header=TRUE,skip=3,comment.char="");
        y.pred[tval,"appl.domain"]<-report.dat$Z.score
      }
    }
    ## Replace with tabulated values if available
    if (intvdat) {
      tvdat.now<-tvdat[jmatch,]
      if (tvdat.now$RfD_Activity != "##") {
        y.pred[ToxValuesNames[1],]$log10prediction <- as.numeric(tvdat.now$RfD_Activity_LogMole)
        y.pred[ToxValuesNames[1],]$log10lower95 <- NA
        y.pred[ToxValuesNames[1],]$log10upper95 <- NA
        y.pred[ToxValuesNames[1],]$prediction <- as.numeric(tvdat.now$RfD_Activity)
        y.pred[ToxValuesNames[1],]$lower95 <- NA
        y.pred[ToxValuesNames[1],]$upper95 <- NA
        y.pred[ToxValuesNames[1],]$appl.domain <- NA
        y.pred[ToxValuesNames[1],]$Source <- tvdat.now$RfD_Source
      }
      if (tvdat.now$RfD_EffectLevel == "NOAEL" | tvdat.now$RfD_EffectLevel == "NOEL" |
          tvdat.now$RfD_EffectLevel == "NOAEL(ADJ)") {
        y.pred[ToxValuesNames[2],]$log10prediction <- as.numeric(tvdat.now$RfD_POD_LogMole)
        y.pred[ToxValuesNames[2],]$log10lower95 <- NA
        y.pred[ToxValuesNames[2],]$log10upper95 <- NA
        y.pred[ToxValuesNames[2],]$prediction <- as.numeric(tvdat.now$RfD_POD)
        y.pred[ToxValuesNames[2],]$lower95 <- NA
        y.pred[ToxValuesNames[2],]$upper95 <- NA
        y.pred[ToxValuesNames[2],]$appl.domain <- NA
        y.pred[ToxValuesNames[2],]$Source <- tvdat.now$RfD_Source
      }
      if (length(grep("BMDL",tvdat.now$RfD_EffectLevel))==1) {
        y.pred[ToxValuesNames[3],]$log10prediction <- as.numeric(tvdat.now$RfD_POD_LogMole)
        y.pred[ToxValuesNames[3],]$log10lower95 <- NA
        y.pred[ToxValuesNames[3],]$log10upper95 <- NA
        y.pred[ToxValuesNames[3],]$prediction <- as.numeric(tvdat.now$RfD_POD)
        y.pred[ToxValuesNames[3],]$lower95 <- NA
        y.pred[ToxValuesNames[3],]$upper95 <- NA
        y.pred[ToxValuesNames[3],]$appl.domain <- NA
        y.pred[ToxValuesNames[3],]$Source <- tvdat.now$RfD_Source
      }
      if (tvdat.now$Oral_Noncancer_BMD != "##") {
        y.pred[ToxValuesNames[4],]$log10prediction <- as.numeric(tvdat.now$Oral_Noncancer_moles)
        y.pred[ToxValuesNames[4],]$log10lower95 <- NA
        y.pred[ToxValuesNames[4],]$log10upper95 <- NA
        y.pred[ToxValuesNames[4],]$prediction <- as.numeric(tvdat.now$Oral_Noncancer_BMD)
        y.pred[ToxValuesNames[4],]$lower95 <- NA
        y.pred[ToxValuesNames[4],]$upper95 <- NA
        y.pred[ToxValuesNames[4],]$appl.domain <- NA
        y.pred[ToxValuesNames[4],]$Source <- tvdat.now$RfD_Source
      }
      if (tvdat.now$OSF_Activity != "##") {
        y.pred[ToxValuesNames[5],]$log10prediction <- as.numeric(tvdat.now$OSF_Activity_LogMole)
        y.pred[ToxValuesNames[5],]$log10lower95 <- NA
        y.pred[ToxValuesNames[5],]$log10upper95 <- NA
        y.pred[ToxValuesNames[5],]$prediction <- as.numeric(tvdat.now$OSF_Activity)
        y.pred[ToxValuesNames[5],]$lower95 <- NA
        y.pred[ToxValuesNames[5],]$upper95 <- NA
        y.pred[ToxValuesNames[5],]$appl.domain <- NA
        y.pred[ToxValuesNames[5],]$Source <- tvdat.now$OSF_Source
      }
      if (tvdat.now$CPV_Activity != "##") {
        y.pred[ToxValuesNames[6],]$log10prediction <- as.numeric(tvdat.now$CPV_Activity_LogMole)
        y.pred[ToxValuesNames[6],]$log10lower95 <- NA
        y.pred[ToxValuesNames[6],]$log10upper95 <- NA
        y.pred[ToxValuesNames[6],]$prediction <- as.numeric(tvdat.now$CPV_Activity)
        y.pred[ToxValuesNames[6],]$lower95 <- NA
        y.pred[ToxValuesNames[6],]$upper95 <- NA
        y.pred[ToxValuesNames[6],]$appl.domain <- NA
        y.pred[ToxValuesNames[6],]$Source <- tvdat.now$CPV_Source
      }
      if (tvdat.now$RfC_Activity != "##") {
        y.pred[ToxValuesNames[7],]$log10prediction <- as.numeric(tvdat.now$RfC_Activity_LogMole)
        y.pred[ToxValuesNames[7],]$log10lower95 <- NA
        y.pred[ToxValuesNames[7],]$log10upper95 <- NA
        y.pred[ToxValuesNames[7],]$prediction <- as.numeric(tvdat.now$RfC_Activity)
        y.pred[ToxValuesNames[7],]$lower95 <- NA
        y.pred[ToxValuesNames[7],]$upper95 <- NA
        y.pred[ToxValuesNames[7],]$appl.domain <- NA
        y.pred[ToxValuesNames[7],]$Source <- tvdat.now$RfC_Source
      }
      if (tvdat.now$IUR_Activity != "##") {
        y.pred[ToxValuesNames[8],]$log10prediction <- as.numeric(tvdat.now$IUR_Activity_LogMole)
        y.pred[ToxValuesNames[8],]$log10lower95 <- NA
        y.pred[ToxValuesNames[8],]$log10upper95 <- NA
        y.pred[ToxValuesNames[8],]$prediction <- as.numeric(tvdat.now$IUR_Activity)
        y.pred[ToxValuesNames[8],]$lower95 <- NA
        y.pred[ToxValuesNames[8],]$upper95 <- NA
        y.pred[ToxValuesNames[8],]$appl.domain <- NA
        y.pred[ToxValuesNames[8],]$Source <- tvdat.now$IUR_Source
      }
    }    
    ##
    y.pred.df <- rbind(y.pred.df,y.pred)  
  }
}
write.csv(y.pred.df,file=paste(infile,"ToxValuePredictions_1-100.csv",sep="_"),row.names = FALSE)
