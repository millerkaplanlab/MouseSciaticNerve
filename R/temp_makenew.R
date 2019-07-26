library(scClustViz)
library(Seurat)

filemap <- data.frame(old=c("inst/Inj9dBeads/Inj9dBeads_old.RData",
                            "inst/Inj9dBeadsMesenchymal/Inj9dBeadsMesenchymal_old.RData",
                            "inst/Inj9dFACS/Inj9dFACS_old.RData",
                            "inst/InjCombined/InjCombined_old.RData",
                            "inst/InjUninjMesenchymalCombined/InjUninjMesenchymalCombined_old.RData",
                            "inst/UninjMesenchymal/UninjMesenchymal_old.RData"),
                      new=c("../MouseSciaticNerve_eb1S/9dInjuredBeads/9dInjuredBeads_eb1S.RData",
                            "../MouseSciaticNerve_eb1S/9dInjuredBeadsMesenchymal/9dInjuredBeadsMesenchymal_eb1S.RData",
                            "../MouseSciaticNerve_eb1S/9dInjuredFACS/9dInjuredFACS_eb1S.RData",
                            "../MouseSciaticNerve_eb1S/CombinedInjured/CombinedInjured_eb1S.RData",
                            "../MouseSciaticNerve_eb1S/InjUninjMesCombined/InjUninjMesCombined_eb1S.RData",
                            "../MouseSciaticNerve_eb1S/InjUninjMesCombinedBatchCorrected/InjUninjMesCombinedBatchCorrected_eb1S.RData"),
                      v2=c("../MouseSciaticNerve_eb1S/9dInjuredBeads/9dInjuredBeads_eb1S_v2.RData",
                           "../MouseSciaticNerve_eb1S/9dInjuredBeadsMesenchymal/9dInjuredBeadsMesenchymal_eb1S_v2.RData",
                           "../MouseSciaticNerve_eb1S/9dInjuredFACS/9dInjuredFACS_eb1S_v2.RData",
                           "../MouseSciaticNerve_eb1S/CombinedInjured/CombinedInjured_eb1S_v2.RData",
                           "../MouseSciaticNerve_eb1S/InjUninjMesCombined/InjUninjMesCombined_eb1S_v2.RData",
                           "../MouseSciaticNerve_eb1S/InjUninjMesCombinedBatchCorrected/InjUninjMesCombinedBatchCorrected_eb1S_v2.RData"),
                      stringsAsFactors=F)
# for (i in 1:6) {
#   load(filemap$new[i])
#   eb1S <- UpdateSeuratObject(eb1S)
#   save(eb1S,file=filemap$v2[i])
# }


i <- 6
temp_name <- strsplit(filemap$old[i],"/",fixed=T)[[1]][2]
load(filemap$new[i])
load(filemap$v2[i])
load(filemap$old[i])
load(sub("old.RData","savedRes.RData",filemap$old[i]))
temp_cells <- sapply(colnames(data_for_scClustViz$nge),function(X) grep(X,colnames(getExpr(eb1S)),value=T))
temp_cl <- data_for_scClustViz$cl[names(temp_cells),savedRes,drop=F]
rownames(temp_cl) <- temp_cells
temp_seur <- DietSeurat(UpdateSeuratObject(eb1S),counts=F)
temp_seur <- subset(temp_seur,cells=temp_cells)
temp_seur@meta.data <- pDat[temp_cells,]
temp_scv <- CalcAllSCV(inD=temp_seur,clusterDF=temp_cl,assayType="RNA")
assign(paste(temp_name,"seur",sep="_"),temp_seur)
assign(paste(temp_name,"sCVdL",sep="_"),temp_scv)
save(list=c(paste(temp_name,"seur",sep="_"),paste(temp_name,"sCVdL",sep="_")),
     file=paste0("inst/",temp_name,"/",temp_name,".RData"))


for (i in 1:5) {
  temp_name <- strsplit(filemap$old[i],"/",fixed=T)[[1]][2]
  load(filemap$new[i])
  load(filemap$v2[i])
  load(filemap$old[i])
  load(sub("old.RData","savedRes.RData",filemap$old[i]))
  temp_seur <- DietSeurat(UpdateSeuratObject(eb1S),counts=F)
  temp_cl <- getMD(temp_seur)[,savedRes,drop=F]
  if (all(rownames(pDat) == rownames(getMD(temp_seur)))) {
    temp_seur@meta.data <- pDat
  } else {
    stop("rownames don't line up")
  }
  temp_scv <- CalcAllSCV(inD=temp_seur,clusterDF=temp_cl,assayType="RNA")
  assign(paste(temp_name,"seur",sep="_"),temp_seur)
  assign(paste(temp_name,"sCVdL",sep="_"),temp_scv)
  save(list=c(paste(temp_name,"seur",sep="_"),paste(temp_name,"sCVdL",sep="_")),
       file=paste0("inst/",temp_name,"/",temp_name,".RData"))
}
