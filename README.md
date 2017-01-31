# KT003_explore_data

```r
rm(list=ls())
setwd("~/Desktop/Microarray Analyses/Kevin's array/151118 data/") #set working directory (folder where Apac array files are stored)
```

##########################################################################################
## Import data

```r
print(load("./Data/Processed Data/0002.KTarray.Apac.cleaned.RData"))
files <- sort(list.files("./Data/Processed Data")[grep("import_data.RData", list.files("./Data/Processed Data"))], decreasing=T)[1]
print(load(paste("./Data/Processed Data/", files, sep=""))) 
```

##########################################################################################
##get dataframes into the same order

```r
annAnti <- annAnti[order(annAnti$Unique.ID),]
pheSera <- pheSera[order(pheSera$Sample),]
MNA <- raw_mean[order(rownames(raw_mean)),order(colnames(raw_mean))]
MedA <- raw_median[order(rownames(raw_median)),order(colnames(raw_median))]

pheSera$SampleID <- paste(gsub("Apac_", "", pheSera$study), pheSera$Sample, sep="_")
```


##########################################################################################
##remove "bad" data

```r
bad <- d1[d1$Flag=="Bad",]
for(row in 1:nrow(bad)){
  i <- bad[row,"study"]
  j <- bad[row,"sampleID"]
  k <- bad[row,"spot"]
  print(paste(gsub("Apac_", "", i), j, sep="_"))
  print(k)
  MNA[rownames(MNA)==k, colnames(MNA)==paste(gsub("Apac_", "", i), j, sep="_")] <- NA
  MedA[rownames(MedA)==k, colnames(MedA)==paste(gsub("Apac_", "", i), j, sep="_")] <- NA
}
```

##########################################################################################
##spot blank_104 is likely the missing BKPAR spot --> rename

```r
annAnti$antigen[annAnti$Unique.ID=="blank_104"] <- annAnti$antigen[annAnti$Unique.ID=="BKPAR_195"]
annAnti$description[annAnti$Unique.ID=="blank_104"] <- annAnti$description[annAnti$Unique.ID=="BKPAR_195"]
annAnti$Unique.ID[annAnti$Unique.ID=="blank_104"] <- "BKPAR_104"
annAnti$ID[annAnti$ID=="blank_104"] <- "BKPAR_104"

rownames(MNA) <- gsub("blank_104", "BKPAR_104", rownames(MNA))
rownames(MedA) <- gsub("blank_104", "BKPAR_104", rownames(MedA))
```

##########################################################################################
##explore data (blanks)

```r
blanks <- annAnti$Unique.ID[annAnti$description=="Blank"]

for(i in unique(pheSera$SampleID)){
  pheSera[pheSera$SampleID==i, "mean_blanks.MNA"] <- mean(MNA[rownames(MNA) %in% blanks, colnames(MNA)==i], na.rm=T)
  pheSera[pheSera$SampleID==i, "sd_blanks.MNA"] <- sd(MNA[rownames(MNA) %in% blanks, colnames(MNA)==i], na.rm=T)
  pheSera[pheSera$SampleID==i, "mean_blanks.MedA"] <- mean(MedA[rownames(MedA) %in% blanks, colnames(MedA)==i], na.rm=T)
  pheSera[pheSera$SampleID==i, "sd_blanks.MedA"] <- sd(MedA[rownames(MedA) %in% blanks, colnames(MedA)==i], na.rm=T)
}

pheSera[pheSera$mean_blanks.MNA > 8000,]
pheSera[pheSera$sd_blanks.MNA > 8000,]
pheSera[pheSera$mean_blanks.MedA > 8000,]
pheSera[pheSera$sd_blanks.MedA > 8000,]
```

###there seems to be a problem with Apac X1 slide 28 and Apac X2 slide 19 --> Flag these
```r
pheSera$Flag[(pheSera$study=="Apac_X1" & pheSera$Slide==28)|(pheSera$study=="Apac_X2" & pheSera$Slide==19)] <- "High Blank"
```

##########################################################################################
##explore data (PBS)

```r
pbs <- annAnti$Unique.ID[annAnti$description=="PBS"]

for(i in unique(pheSera$SampleID)){
  pheSera[pheSera$SampleID==i, "mean_pbs.MNA"] <- mean(MNA[rownames(MNA) %in% pbs, colnames(MNA)==i], na.rm=T)
  pheSera[pheSera$SampleID==i, "sd_pbs.MNA"] <- sd(MNA[rownames(MNA) %in% pbs, colnames(MNA)==i], na.rm=T)
  pheSera[pheSera$SampleID==i, "mean_pbs.MedA"] <- mean(MedA[rownames(MedA) %in% pbs, colnames(MedA)==i], na.rm=T)
  pheSera[pheSera$SampleID==i, "sd_pbs.MedA"] <- sd(MedA[rownames(MedA) %in% pbs, colnames(MedA)==i], na.rm=T)
}

pheSera[pheSera$mean_pbs.MNA > 8000,]
pheSera[pheSera$sd_pbs.MNA > 8000,]
pheSera[pheSera$mean_pbs.MedA > 8000,]
pheSera[pheSera$sd_pbs.MedA > 8000,]

#nothing obvious going on with the few cases with high PBS values
```

##########################################################################################
##cleanup workspace and save
```r
ls()
save(annAnti, pheSera, MNA, MedA, d1, file="./Data/Processed Data/0003.KTarray.Apac.explored.RData")
```
