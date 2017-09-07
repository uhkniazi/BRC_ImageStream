#install.packages("PerformanceAnalytics")
#install.packages("ggplot2")
#install.packages("corrplot")

library(ggplot2) 
source("http://www.sthda.com/upload/rquery_cormat.r")
library("PerformanceAnalytics")

rm(list = ls())  

setwd("C://Users/derinald/Desktop/ImageStream Data Set - October 2016/data")
main_path = "C://Users/derinald/Desktop/ImageStream Data Set - October 2016/"
input_dir=  paste(main_path, "data", sep="/")
out_dir=paste(main_path, "R/results", sep="/")

file_list <- list.files()
file_list<-file_list[grepl("^[0-9]+.txt", file_list)]

if(exists("dataset")){rm(dataset)}
for (file in file_list){
print(file)
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t", comment.char = "", na.strings = c("#VALUE!", "NaN"))

  }
  # if the merged dataset does exist, append to it
  else{
    temp_dataset <-read.table(file, header=TRUE, sep="\t", comment.char = "", na.strings = c("#VALUE!", "NaN"),colClasses=c(rep("character", 6), rep("numeric", 2),
    "character" , rep("numeric", 2), rep("NULL", count.fields(file,sep="\t")[2]-dim(dataset)[2])))
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  print(sum(duplicated(dataset)))
}
pat.data <-read.table("PASI report ISX July 2017.txt", header=TRUE, sep="\t", comment.char = "", na.strings = c("#VALUE!", "NaN"))
colnames(pat.data)[1] <- colnames(dataset)[1]
dataset<-merge(dataset,pat.data,by="Patient.ID")

drug.data <-read.table("Drug Levels.July2017.txt", header=TRUE, sep="\t", comment.char = "", na.strings = c("#VALUE!", "NaN"))
colnames(drug.data) <-c("Patient.ID", "Visit..Week.", "ADL.conc","USK.conc")
dataset$tmp<- paste(dataset$Patient.ID,dataset$Visit..Week.,sep=".")
drug.data$tmp<- paste(drug.data$Patient.ID,drug.data$Visit..Week.,sep=".")
dataset<-merge(dataset, drug.data, by="tmp")
dataset<-dataset[,!(names(dataset) %in% c("tmp","Patient.ID.y", "Visit..Week..y"))]
colnames(dataset)<-sub("Visit..Week..x", "Visit..Week",colnames(dataset))
colnames(dataset)<-sub("Patient.ID.x", "Patient.ID",colnames(dataset))

## CORRECTING MISPELLED WORDS
dataset$Transcription.factor<-as.factor(gsub("\\?", "k",dataset$Transcription.factor))
dataset$Stimulation<-as.factor(gsub("TNF\\?", "TNF-alpha",dataset$Stimulation))
dataset$Stimulation<-as.factor(gsub("IFN\\?", "IFN-gamma",dataset$Stimulation))


### AD-HOC CHANGES BASED ON ROSA'S INSTRUCTIONS
### Monocyte measurements should exclude patients 006-006 to 006-019 due to an issue we had with the protocol
p=c("006-006","006-007","006-008","006-009","006-010","006-011","006-012","006-013","006-014","006-015","006-016","006-017","006-018","006-019")
ix=(dataset$Patient.ID %in% p) & (dataset$Cell.type=="Monocytes")
dataset=dataset[!ix,]

# Patients 006-001 to 006-005 have different T cell subsets measured (the name of the cell types is different for T cells),
# given our time frame, I think it is better to exclude them for T cell measurements for this quick analysis.
p=c("006-001","006-002","006-003","006-004","006-005","006-006")
ix=(dataset$Patient.ID %in% p) & (dataset$Cell.type %in% c("T cells", "CLA+ T cells"))
dataset=dataset[!ix,]


### Additional ad-hoc modifications and additions
dataset$Week.12.PASI<- as.numeric(as.character(dataset$Week.12.PASI))
dataset$Week.4.PASI <- as.numeric(as.character(dataset$Week.4.PASI))
dataset$Week.1.PASI <- as.numeric(as.character(dataset$Week.1.PASI))
dataset$relativePASI<-dataset$Week.12.PASI/dataset$Baseline.PASI


dataset$PASI90<-dataset$Week.12.PASI < 0.10 *dataset$Baseline.PASI
dataset$PASI75<-dataset$Week.12.PASI < 0.25 *dataset$Baseline.PASI
dataset$PASI50<-dataset$Week.12.PASI < 0.50 *dataset$Baseline.PASI
dataset$relativePASI<-dataset$Week.12.PASI/dataset$Baseline.PASI
dataset$Rd.score<-as.numeric(as.character(dataset$Rd.score))
#dataset<-dataset[dataset$Cell.count > 10, ]     ##GETS READ OF RECORDS WHERE NOT ENOUGH CELLS WERE CAPTURED

##EXCLUDE PARTS OF THE EXPERIMENTS BECAUSE OF EXPERIMENTAL STAINING ISSUES  or the cytokine stimulation

#Experiments to be excluded due to staining problems:
#Patient ID	Visit	Transcription factor to be excluded
#006020	Week 12	NFkB, STAT3, STAT4      
#006022	Week 12	NFkB, STAT3, STAT4
#006025	Week 4	NFkB, STAT3, STAT4
#006030	Week 12	NFkB, STAT3, STAT4
#006032	Baseline	NFkB
#006033	Week 12	NFkB, STAT3, STAT4
#006036	Week 4	STAT3
#060002 	Week 12	NFkB
#060003	Week 4	NFkB
#    060 005 04: all NFkB signals
#    060 006 00: all TNF and TNF + IL17 signals 

ix0<- which(dataset$Patient.ID=="006-020" & dataset$Visit..Week=="Week 12"  & dataset$Transcription.factor %in% c("NF-kB", "STAT3", "STAT4"))
ix1<- which(dataset$Patient.ID=="006-022" & dataset$Visit..Week=="Week 12"  & dataset$Transcription.factor %in% c("NF-kB", "STAT3", "STAT4"))
ix2<- which(dataset$Patient.ID=="006-025" & dataset$Visit..Week=="Week 4"  & dataset$Transcription.factor %in% c("NF-kB", "STAT3", "STAT4"))
ix3<- which(dataset$Patient.ID=="006-030" & dataset$Visit..Week=="Week 12"  & dataset$Transcription.factor %in% c("NF-kB", "STAT3", "STAT4"))
ix4<- which(dataset$Patient.ID=="006-032" & dataset$Visit..Week=="Baseline"  & dataset$Transcription.factor %in% c("NF-kB"))
ix5<- which(dataset$Patient.ID=="006-033" & dataset$Visit..Week=="Week 12"  & dataset$Transcription.factor %in% c("NF-kB", "STAT3", "STAT4"))
ix6<- which(dataset$Patient.ID=="006-036" & dataset$Visit..Week=="Week 4"  & dataset$Transcription.factor %in% c("STAT3"))
ix7<- which(dataset$Patient.ID=="006-002" & dataset$Visit..Week=="Week 12"  & dataset$Transcription.factor %in% c("NF-kB"))
ix8<- which(dataset$Patient.ID=="006-003" & dataset$Visit..Week=="Week 4"  & dataset$Transcription.factor %in% c("NF-kB"))
ix9<- which(dataset$Patient.ID=="060-005" & dataset$Transcription.factor %in% c("NF-kB") & (dataset$Visit..Week == "Week 4"))
ix10<-which(dataset$Patient.ID=="060-006" & dataset$Stimulation %in% c("TNF-alpha", "TNF-alpha + IL-17") & (dataset$Visit..Week == "Baseline"))


ix=sort(unique(c(ix0,ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8,ix9,ix10)))
dataset<-dataset[-ix,]
write.table(dataset,row.names = FALSE, file = "merged.files_and_annotations.txt", sep="\t")
