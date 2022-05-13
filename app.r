options(shiny.maxRequestSize = 50*1024^2)
############################################
############################################
### BAGS - Barcode, Audit & Grade system ###
############################################
############################################

##### INSTALL NECESSARY PACKAGES
#install.packages(c("seqRFLP","bold","data.table","worms","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase"))

#RUN FROM GITHUB:
#1
#RUN:
#library(shiny)
#2
#RUN: 
#runGitHub("BAGS", "tadeu95")

##### LOAD NECESSARY PACKAGES 
library(seqRFLP)
library(bold)
library(data.table)
library(worms)
library(stringr)
library(readr)
library(fingerprint)
library(dplyr)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(snakecase)

#GET FILES FROM GITHUB
file_spb<-"https://raw.githubusercontent.com/tadeu95/BAGS/master/species_per_bin.txt"
spb<-fread(file_spb)
file_bps<-"https://raw.githubusercontent.com/tadeu95/BAGS/master/bin_per_species.txt"
bps<-fread(file_bps)


#####GRADES CHECKLIST
grades_checklist<-function(checklist,inputz,coordz){
  species<-unique(checklist[,1])
  x<-length(species)
  taxon_total = data.frame()
  y=x%/%300+2
  i <- 1
  while (i < y){
    ini<-1+(300*(i-1))
    fin<-300*i
    tmp <- bold_seqspec(taxon=species[ini:fin], response = TRUE)
    tt <- paste0(rawToChar(tmp$content, multiple = TRUE), collapse = "")
    Encoding(tt) <- "UTF-8"
    taxa <- utils::read.delim(text = tt, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    taxon_total <- rbind(taxon_total,taxa)
    i = i+1
  }
  taxon<-taxon_total
  taxon2<-taxon[taxon$species_name!=""|is.na(taxon$species_name),]
  taxon2<-taxon2[!(taxon2$bin_uri == "" | is.na(taxon2$bin_uri)), ]
  taxon2<-left_join(taxon2,bps,by="species_name")
  #taxon2$bin_uri.y=NULL
  names(taxon2)[names(taxon2) == "bin_uri.x"] <- "bin_uri"
  taxon2<-left_join(taxon2,spb,by="bin_uri")
  taxon2$species_name.y=NULL
  names(taxon2)[names(taxon2) == "species_name.x"] <- "species_name"
  taxon3<-taxon2[taxon2$markercode=="COI-5P",]
  taxon3$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon3$nucleotides)
  taxon3$nucleotides=gsub("-","",taxon3$nucleotides)
  patterns_non<-c("[0-9]+.+","[a-z,A-Z]+[0-9]+.+","^[A-Z,a-z] "," [A-Z,a-z]$","[0-9]+.*","[a-z,A-Z]+[0-9]+.*","[0-9]+"," +$"," cmplx$"," [A-Z]+.+$"," cmplx$")
  for (p in patterns_non){
    taxon3$species_name<-gsub(p,"",taxon3$species_name)
  }
  patterns_fixed<-c("-"," sp."," sp. nov"," complex."," f."," nr."," s.l."," grp."," type"," group")
  for (p in patterns_fixed){
    taxon3$species_name<-gsub(p,"",taxon3$species_name,fixed=TRUE)
  }
  taxon8<-data.frame(taxon3$species_name,taxon3$bin_uri,taxon3$nucleotides,taxon3$country,taxon3$family_name,taxon3$order_name,taxon3$class_name,taxon3$sampleid,taxon3$processid,taxon3$species_per_bin,taxon3$bin_per_species, taxon3$lat,taxon3$lon)
  names(taxon8)<-c("species_name","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude")
  taxon8$grade=NA
  names(taxon8)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","grade")
  taxon9=data.frame(taxon8$species,taxon8$BIN,taxon8$sequence,taxon8$country,taxon8$grade,taxon8$family,taxon8$order,taxon8$class,taxon8$sampleid,taxon8$processid,taxon8$species_per_bin,taxon8$bin_per_species,taxon8$lat,taxon8$lon)
  names(taxon9)=c("species","BIN","sequence","country","grade","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude")
  taxon19<-taxon9 %>%
    mutate(grade = ifelse(species_per_bin>1,"E",
                          ifelse(bin_per_species>1 & species_per_bin==1,"C",
                                 ifelse(bin_per_species==1 & species_per_bin==1,"AB","needs_update"))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "AB"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  taxon19$contains_dominant=NULL
  if (coordz){
    taxon19<-taxon19[!(is.na(taxon19$lat)) | taxon19$country!="",]
  }
  taxon19$base_number=str_count(taxon19$sequence, pattern="[A-Z]")
  taxon19<-taxon19[(taxon19$base_number>=inputz),]
  np=(str_count(taxon19$sequence, "N")/str_count(taxon19$sequence, "[A-Z]"))*100
  taxon19$n_percent=np
  taxon19<-subset(taxon19,taxon19$n_percent<1)
  taxon19$n_percent=NULL
  taxon19=subset(taxon19, lengths(gregexpr("\\w+", taxon19$species)) > 1)
  num_species=table(taxon19$species)
  num_species=as.data.frame(num_species)
  names(num_species)=c("species","frequency_species")
  taxon19<-inner_join(taxon19,num_species)
  taxon19<-taxon19 %>%
    mutate(grade = ifelse(grade=="E","E",
                          ifelse(frequency_species<3,"D",
                                 ifelse(grade=="C","C",
                                        ifelse(grade=="AB" & frequency_species<11,"B",
                                               ifelse(grade=="AB" & frequency_species>10,"A","needs_update"))))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "D"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "B"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "A"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  taxon19$contains_dominant=NULL
  taxon19$BIN_per_species=NULL
  taxon19$species_per_bin=NULL
  taxon19$species_frequency=NULL
  taxon19$species_per_bin=NULL
  taxon19$bin_per_species=NULL
  taxon19$frequency_species=NULL
  taxon19<-taxon19[order(taxon19$species),]
  assign('taxon19',taxon19,envir=.GlobalEnv)
}
#####IMPLEMENT RANKING SYSTEM FOR ALL TAXA
grades2<-function(groups,inputz,coordz){
  tmp <- bold_seqspec(taxon=groups, response = TRUE)
  tt <- paste0(rawToChar(tmp$content, multiple = TRUE), collapse = "")
  Encoding(tt) <- "UTF-8"
  taxon <- utils::read.delim(text = tt, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  taxon2<-taxon[taxon$species_name!=""|is.na(taxon$species_name),]
  taxon2<-taxon2[!(taxon2$bin_uri == "" | is.na(taxon2$bin_uri)), ]
  taxon2<-left_join(taxon2,bps,by="species_name")
  taxon2$bin_uri.y=NULL
  names(taxon2)[names(taxon2) == "bin_uri.x"] <- "bin_uri"
  taxon2<-left_join(taxon2,spb,by="bin_uri")
  taxon2$species_name.y=NULL
  names(taxon2)[names(taxon2) == "species_name.x"] <- "species_name"
  taxon3<-taxon2[taxon2$markercode=="COI-5P",]
  taxon3$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon3$nucleotides)
  taxon3$nucleotides=gsub("-","",taxon3$nucleotides)
  patterns_non<-c("[0-9]+.+","[a-z,A-Z]+[0-9]+.+","^[A-Z,a-z] "," [A-Z,a-z]$","[0-9]+.*","[a-z,A-Z]+[0-9]+.*","[0-9]+"," +$"," cmplx$"," [A-Z]+.+$"," cmplx$")
  for (p in patterns_non){
    taxon3$species_name<-gsub(p,"",taxon3$species_name)
  }
  patterns_fixed<-c("-"," sp."," sp. nov"," complex."," f."," nr."," s.l."," grp."," type"," group")
  for (p in patterns_fixed){
    taxon3$species_name<-gsub(p,"",taxon3$species_name,fixed=TRUE)
  }
  taxon8<-data.frame(taxon3$species_name,taxon3$bin_uri,taxon3$nucleotides,taxon3$country,taxon3$family_name,taxon3$order_name,taxon3$class_name,taxon3$sampleid,taxon3$processid,taxon3$species_per_bin,taxon3$bin_per_species, taxon3$lat,taxon3$lon)
  names(taxon8)<-c("species_name","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude")
  taxon8$grade=NA
  names(taxon8)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","grade")
  taxon9=data.frame(taxon8$species,taxon8$BIN,taxon8$sequence,taxon8$country,taxon8$grade,taxon8$family,taxon8$order,taxon8$class,taxon8$sampleid,taxon8$processid,taxon8$species_per_bin,taxon8$bin_per_species,taxon8$lat,taxon8$lon)
  names(taxon9)=c("species","BIN","sequence","country","grade","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude")
  taxon19<-taxon9 %>%
    mutate(grade = ifelse(species_per_bin>1,"E",
                          ifelse(bin_per_species>1 & species_per_bin==1,"C",
                                 ifelse(bin_per_species==1 & species_per_bin==1,"AB","needs_update"))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "AB"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  taxon19$contains_dominant=NULL
  if (coordz){
    taxon19<-taxon19[!(is.na(taxon19$lat)) | taxon19$country!="",]
  }
  taxon19$base_number=str_count(taxon19$sequence, pattern="[A-Z]")
  taxon19<-taxon19[(taxon19$base_number>inputz),]
  np=(str_count(taxon19$sequence, "N")/str_count(taxon19$sequence, "[A-Z]"))*100
  taxon19$n_percent=np
  taxon19<-subset(taxon19,taxon19$n_percent<1)
  taxon19$n_percent=NULL
  taxon19=subset(taxon19, lengths(gregexpr("\\w+", taxon19$species)) > 1)
  num_species=table(taxon19$species)
  num_species=as.data.frame(num_species)
  names(num_species)=c("species","frequency_species")
  taxon19<-inner_join(taxon19,num_species)
  taxon19<-taxon19 %>%
    mutate(grade = ifelse(grade=="E","E",
                          ifelse(frequency_species<3,"D",
                                 ifelse(grade=="C","C",
                                        ifelse(grade=="AB" & frequency_species<11,"B",
                                               ifelse(grade=="AB" & frequency_species>10,"A","needs_update"))))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "D"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "B"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "A"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  taxon19$contains_dominant=NULL
  taxon19$BIN_per_species=NULL
  taxon19$species_per_bin=NULL
  taxon19$species_frequency=NULL
  taxon19$species_per_bin=NULL
  taxon19$bin_per_species=NULL
  taxon19$frequency_species=NULL
  taxon19<-taxon19[order(taxon19$species),]
  assign('taxon19',taxon19,envir=.GlobalEnv)
}
####IMPLEMENT RANKING SYSTEM MARINE TAXA
grades<-function(groups,inputz,coordz){
  tmp <- bold_seqspec(taxon=groups, response = TRUE)
  tt <- paste0(rawToChar(tmp$content, multiple = TRUE), collapse = "")
  Encoding(tt) <- "UTF-8"
  taxon <- utils::read.delim(text = tt, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  taxon2<-taxon[taxon$species_name!=""|is.na(taxon$species_name),]
  taxon2<-taxon2[!(taxon2$bin_uri == "" | is.na(taxon2$bin_uri)), ]
  taxon2<-left_join(taxon2,bps,by="species_name")
  taxon2$bin_uri.y=NULL
  names(taxon2)[names(taxon2) == "bin_uri.x"] <- "bin_uri"
  taxon2<-left_join(taxon2,spb,by="bin_uri")
  taxon2$species_name.y=NULL
  names(taxon2)[names(taxon2) == "species_name.x"] <- "species_name"
  taxon3<-taxon2[taxon2$markercode=="COI-5P",]
  taxon3$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon3$nucleotides)
  taxon3$nucleotides=gsub("-","",taxon3$nucleotides,fixed=TRUE)
  patterns_non<-c("[0-9]+.+","[a-z,A-Z]+[0-9]+.+","^[A-Z,a-z] "," [A-Z,a-z]$","[0-9]+.*","[a-z,A-Z]+[0-9]+.*","[0-9]+"," +$"," cmplx$"," [A-Z]+.+$"," cmplx$")
  for (p in patterns_non){
    taxon3$species_name<-gsub(p,"",taxon3$species_name)
  }
  patterns_fixed<-c("-"," sp."," sp. nov"," complex."," f."," nr."," s.l."," grp."," type"," group")
  for (p in patterns_fixed){
    taxon3$species_name<-gsub(p,"",taxon3$species_name,fixed=TRUE)
  }
  taxon4<-wormsbynames(as.character(unique(taxon3$species_name)), marine_only=FALSE, ids=TRUE)
  taxon4<-taxon4%>%
    mutate(new_name = ifelse(status=="unaccepted" & !is.na(status),valid_name,
                             ifelse(status=="accepted" & !is.na(status),name,name)))
  colnames(taxon4)[colnames(taxon4) == "name"] <- "species_name"
  taxon5<-subset(taxon4, taxon4$isMarine==1 | taxon4$isBrackish==1)
  taxon3<-left_join(taxon3,taxon5,by="species_name")
  taxon6<-subset(taxon3,as.character(taxon3$species_name)%in% as.character(taxon5$species_name))
  taxon7<-taxon6
  taxon8<-data.frame(taxon7$species_name,taxon7$bin_uri,taxon7$nucleotides,taxon7$country,taxon7$family_name,taxon7$order_name,taxon7$class_name,taxon7$sampleid,taxon7$processid,taxon7$species_per_bin,taxon7$bin_per_species, taxon7$lat,taxon7$lon,taxon7$new_name)
  names(taxon8)<-c("species_name","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","worms_accepted_name")
  taxon8$grade=NA
  names(taxon8)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","worms_accepted_name","grade")
  taxon9=data.frame(taxon8$species,taxon8$BIN,taxon8$sequence,taxon8$country,taxon8$grade,taxon8$family,taxon8$order,taxon8$class,taxon8$sampleid,taxon8$processid,taxon8$species_per_bin,taxon8$bin_per_species,taxon8$lat,taxon8$lon,taxon8$worms_accepted_name)
  names(taxon9)=c("species","BIN","sequence","country","grade","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","worms_accepted_name")
  taxon19<-taxon9 %>%
    mutate(grade = ifelse(species_per_bin>1,"E",
                          ifelse(bin_per_species>1 & species_per_bin==1,"C",
                                 ifelse(bin_per_species==1 & species_per_bin==1,"AB","needs_update"))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "AB"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  taxon19$contains_dominant=NULL
  if (coordz){
    taxon19<-taxon19[!(is.na(taxon19$lat)) | taxon19$country!="",]
  }
  taxon19$base_number=str_count(taxon19$sequence, pattern="[A-Z]")
  taxon19<-taxon19[(taxon19$base_number>inputz),]
  np=(str_count(taxon19$sequence, "N")/str_count(taxon19$sequence, "[A-Z]"))*100
  taxon19$n_percent=np
  taxon19<-subset(taxon19,taxon19$n_percent<1)
  taxon19$n_percent=NULL
  taxon19=subset(taxon19, lengths(gregexpr("\\w+", taxon19$species)) > 1)
  num_species=table(taxon19$species)
  num_species=as.data.frame(num_species)
  names(num_species)=c("species","frequency_species")
  taxon19<-inner_join(taxon19,num_species)
  taxon19<-taxon19 %>%
    mutate(grade = ifelse(grade=="E","E",
                          ifelse(frequency_species<3,"D",
                                 ifelse(grade=="C","C",
                                        ifelse(grade=="AB" & frequency_species<11,"B",
                                               ifelse(grade=="AB" & frequency_species>10,"A","needs_update"))))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "D"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "B"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "A"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  taxon19$contains_dominant=NULL
  taxon19$BIN_per_species=NULL
  taxon19$species_per_bin=NULL
  taxon19$species_frequency=NULL
  taxon19$species_per_bin=NULL
  taxon19$bin_per_species=NULL
  taxon19$frequency_species=NULL
  taxon19$frequency_species=NULL
  taxon19<-taxon19[order(taxon19$species),]
  assign('taxon19',taxon19,envir=.GlobalEnv)
}
####IMPLEMENT RANKING SYSTEM NONMARINE TAXA
grades_nonmarine<-function(groups,inputz,coordz){
  tmp <- bold_seqspec(taxon=groups, response = TRUE)
  tt <- paste0(rawToChar(tmp$content, multiple = TRUE), collapse = "")
  Encoding(tt) <- "UTF-8"
  taxon <- utils::read.delim(text = tt, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  taxon2<-taxon[taxon$species_name!=""|is.na(taxon$species_name),]
  taxon2<-taxon2[!(taxon2$bin_uri == "" | is.na(taxon2$bin_uri)), ]
  taxon2<-left_join(taxon2,bps,by="species_name")
  taxon2$bin_uri.y=NULL
  names(taxon2)[names(taxon2) == "bin_uri.x"] <- "bin_uri"
  taxon2<-left_join(taxon2,spb,by="bin_uri")
  taxon2$species_name.y=NULL
  names(taxon2)[names(taxon2) == "species_name.x"] <- "species_name"
  taxon3<-taxon2[taxon2$markercode=="COI-5P",]
  taxon3$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon3$nucleotides)
  taxon3$nucleotides=gsub("-","",taxon3$nucleotides)
  patterns_non<-c("[0-9]+.+","[a-z,A-Z]+[0-9]+.+","^[A-Z,a-z] "," [A-Z,a-z]$","[0-9]+.*","[a-z,A-Z]+[0-9]+.*","[0-9]+"," +$"," cmplx$"," [A-Z]+.+$"," cmplx$")
  for (p in patterns_non){
    taxon3$species_name<-gsub(p,"",taxon3$species_name)
  }
  patterns_fixed<-c("-"," sp."," sp. nov"," complex."," f."," nr."," s.l."," grp."," type"," group")
  for (p in patterns_fixed){
    taxon3$species_name<-gsub(p,"",taxon3$species_name,fixed=TRUE)
  }
  taxon4<-wormsbynames(as.character(unique(taxon3$species_name)), marine_only=FALSE, ids=TRUE)
  taxon4<-taxon4%>%
    mutate(new_name = ifelse(status=="unaccepted" & !is.na(status),valid_name,
                             ifelse(status=="accepted" & !is.na(status),name,name)))
  colnames(taxon4)[colnames(taxon4) == "name"] <- "species_name"
  taxon5<-subset(taxon4, taxon4$isMarine!=1)
  taxon5b<-subset(taxon5, taxon5$isBrackish!=1)
  taxon3<-left_join(taxon3,taxon5b,by="species_name")
  taxon6<-subset(taxon3,(as.character(taxon3$species_name)%in% as.character(taxon5b$species_name)))
  taxon7<-taxon6
  taxon8<-data.frame(taxon7$species_name,taxon7$bin_uri,taxon7$nucleotides,taxon7$country,taxon7$family_name,taxon7$order_name,taxon7$class_name,taxon7$sampleid,taxon7$processid,taxon7$species_per_bin,taxon7$bin_per_species, taxon7$lat,taxon7$lon,taxon7$new_name)
  names(taxon8)<-c("species_name","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","worms_accepted_name")
  taxon8$grade=NA
  names(taxon8)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","worms_accepted_name","grade")
  taxon9=data.frame(taxon8$species,taxon8$BIN,taxon8$sequence,taxon8$country,taxon8$grade,taxon8$family,taxon8$order,taxon8$class,taxon8$sampleid,taxon8$processid,taxon8$species_per_bin,taxon8$bin_per_species,taxon8$lat,taxon8$lon,taxon8$worms_accepted_name)
  names(taxon9)=c("species","BIN","sequence","country","grade","family","order","class","sampleid","processid","species_per_bin","bin_per_species","lattitude","longitude","worms_accepted_name")
  taxon19<-taxon9 %>%
    mutate(grade = ifelse(species_per_bin>1,"E",
                          ifelse(bin_per_species>1 & species_per_bin==1,"C",
                                 ifelse(bin_per_species==1 & species_per_bin==1,"AB","needs_update"))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "AB"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  taxon19$contains_dominant=NULL
  if (coordz){
    taxon19<-taxon19[!(is.na(taxon19$lat)) | taxon19$country!="",]
  }
  taxon19$base_number=str_count(taxon19$sequence, pattern="[A-Z]")
  taxon19<-taxon19[(taxon19$base_number>inputz),]
  np=(str_count(taxon19$sequence, "N")/str_count(taxon19$sequence, "[A-Z]"))*100
  taxon19$n_percent=np
  taxon19<-subset(taxon19,taxon19$n_percent<1)
  taxon19$n_percent=NULL
  taxon19=subset(taxon19, lengths(gregexpr("\\w+", taxon19$species)) > 1)
  num_species=table(taxon19$species)
  num_species=as.data.frame(num_species)
  names(num_species)=c("species","frequency_species")
  taxon19<-inner_join(taxon19,num_species)
  taxon19<-taxon19 %>%
    mutate(grade = ifelse(grade=="E","E",
                          ifelse(frequency_species<3,"D",
                                 ifelse(grade=="C","C",
                                        ifelse(grade=="AB" & frequency_species<11,"B",
                                               ifelse(grade=="AB" & frequency_species>10,"A","needs_update"))))))
  dominant_grade <- "E"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "D"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "C"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "B"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  dominant_grade <- "A"
  dt <- as.data.table(taxon19)
  dt[, contains_dominant := any(grade == dominant_grade), by=species]
  dt[contains_dominant == TRUE, grade := dominant_grade]
  taxon19 <- setDF(dt)
  taxon19$contains_dominant=NULL
  taxon19$BIN_per_species=NULL
  taxon19$species_per_bin=NULL
  taxon19$species_frequency=NULL
  taxon19$species_per_bin=NULL
  taxon19$bin_per_species=NULL
  taxon19$frequency_species=NULL
  taxon19<-taxon19[order(taxon19$species),]
  assign('taxon19',taxon19,envir=.GlobalEnv)
}
###Barplot of the frequency of each grade assigned.
#Specimens
plot_summary_specimens=function(taxon19){
  frequency_grades=c(length(which(taxon19$grade=="A")),length(which(taxon19$grade=="B")),length(which(taxon19$grade=="C")),length(which(taxon19$grade=="D")),length(which(taxon19$grade=="E")))
  names_grades=c("A","B","C","D","E")
  dataframe_grades=data.frame(names_grades,frequency_grades)
  ggplot(dataframe_grades, aes(main="Number of specimens per grade",x=names_grades, y=frequency_grades, fill=names_grades))+
    geom_bar(stat="identity", width=0.7)+
    labs(x = "Grades", y = "Number of specimens", fill="Grades")+
    scale_fill_manual(values=c("deepskyblue2", "darkolivegreen3", "darkgoldenrod1","darkorange","brown1"))+
    #ggtitle("Specimens per grade")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,vjust=-5),text=element_text(size=12,face="bold"),axis.text=element_text(size=10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA))
}
#Species
plot_summary_species=function(taxon19){
  frequency_grades=c(as.numeric(length(unique(taxon19$species[taxon19$grade=="A"]))),as.numeric(length(unique(taxon19$species[taxon19$grade=="B"]))),as.numeric(length(unique(taxon19$species[taxon19$grade=="C"]))),as.numeric(length(unique(taxon19$species[taxon19$grade=="D"]))),as.numeric(length(unique(taxon19$species[taxon19$grade=="E"]))))
  names_grades=c("A","B","C","D","E")
  dataframe_grades=data.frame(names_grades,frequency_grades)
  ggplot(dataframe_grades, aes(main="Number of species per grade",x=names_grades, y=frequency_grades, fill=names_grades))+
    geom_bar(stat="identity", width=0.7)+
    labs(x = "Grades", y = "Number of species", fill="Grades")+
    scale_fill_manual(values=c("deepskyblue3", "darkolivegreen4", "darkgoldenrod2","darkorange2","brown3"))+
    #ggtitle("Species per grade")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5,vjust=-5),text=element_text(size=12,face="bold"),axis.text=element_text(size=10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA))
}
####Create fasta files individual
#grade A
create_A=function(taxon19){
  taxon_a=taxon19[taxon19$grade=="A",]
  taxon_a=taxon_a[rowSums(is.na(taxon_a)) != ncol(taxon_a), ]
  taxon_a$species=paste(taxon_a$species,taxon_a$BIN,taxon_a$processid,sep=" | ")
  taxon_a2<-data.frame(taxon_a$species,taxon_a$sequence)
  names(taxon_a2)<-c("name","sequence")
  return(taxon_a2)
}
#grade B
create_B=function(taxon19){
  taxon_b=taxon19[taxon19$grade=="B",]
  taxon_b=taxon_b[rowSums(is.na(taxon_b)) != ncol(taxon_b), ]
  taxon_b$species=paste(taxon_b$species,taxon_b$BIN,taxon_b$processid,sep=" | ")
  taxon_b2<-data.frame(taxon_b$species,taxon_b$sequence)
  names(taxon_b2)<-c("name","sequence")
  return(taxon_b2)
}
#grade C
create_C=function(taxon19){
  taxon_c=taxon19[taxon19$grade=="C",]
  taxon_c=taxon_c[rowSums(is.na(taxon_c)) != ncol(taxon_c), ]
  taxon_c$species=paste(taxon_c$species,taxon_c$BIN,taxon_c$processid,sep=" | ")
  taxon_c2<-data.frame(taxon_c$species,taxon_c$sequence)
  names(taxon_c2)<-c("name","sequence")
  return(taxon_c2)
}
#grade D
create_D=function(taxon19){
  taxon_d=taxon19[taxon19$grade=="D",]
  taxon_d=taxon_d[rowSums(is.na(taxon_d)) != ncol(taxon_d), ]
  taxon_d$species=paste(taxon_d$species,taxon_d$BIN,taxon_d$processid,sep=" | ")
  taxon_d2<-data.frame(taxon_d$species,taxon_d$sequence)
  names(taxon_d2)<-c("name","sequence")
  return(taxon_d2)
}
#grade E
create_E=function(taxon19){
  taxon_e=taxon19[taxon19$grade=="E",]
  taxon_e=taxon_e[rowSums(is.na(taxon_e)) != ncol(taxon_e), ]
  taxon_e$species=paste(taxon_e$species,taxon_e$BIN,taxon_e$processid,sep=" | ")
  taxon_e2<-data.frame(taxon_e$species,taxon_e$sequence)
  names(taxon_e2)<-c("name","sequence")
  return(taxon_e2)
}
####Create fasta files grouped
#Only grades A and B
create_AB=function(taxon19){
  taxon_ab=taxon19[taxon19$grade=="A" | taxon19$grade=="B",]
  taxon_ab=taxon_ab[rowSums(is.na(taxon_ab)) != ncol(taxon_ab), ]
  taxon_ab$species=paste(taxon_ab$species,taxon_ab$BIN,taxon_ab$processid,sep=" | ")
  taxon_ab2<-data.frame(taxon_ab$species,taxon_ab$sequence)
  names(taxon_ab2)<-c("name","sequence")
  return(taxon_ab2)
}
#Only grades A , B and C
create_ABC=function(taxon19){
  taxon_abc=taxon19[taxon19$grade=="A" | taxon19$grade=="B" | taxon19$grade=="C",]
  taxon_abc=taxon_abc[rowSums(is.na(taxon_abc)) != ncol(taxon_abc), ]
  taxon_abc$species=paste(taxon_abc$species,taxon_abc$BIN,taxon_abc$sprocessid,sep=" | ")
  taxon_abc2<-data.frame(taxon_abc$species,taxon_abc$sequence)
  names(taxon_abc2)<-c("name","sequence")
  return(taxon_abc2)
}
#Only grades A , B , C and D
create_ABCD=function(taxon19){
  taxon_abcd=taxon19[taxon19$grade=="A" | taxon19$grade=="B" | taxon19$grade=="C" | taxon19$grade=="D",]
  taxon_abcd=taxon_abcd[rowSums(is.na(taxon_abcd)) != ncol(taxon_abcd), ]
  taxon_abcd$species=paste(taxon_abcd$species,taxon_abcd$BIN,taxon_abcd$processid,sep=" | ")
  taxon_abcd2<-data.frame(taxon_abcd$species,taxon_abcd$sequence)
  names(taxon_abcd2)<-c("name","sequence")
  return(taxon_abcd2)
}
#All grades
create_fasta=function(taxon19){
  taxon19=taxon19[!is.na(taxon19$grade), ]
  taxon19$species=paste(taxon19$species, taxon19$BIN,taxon19$processid, sep=" | ")
  taxon20<-data.frame(taxon19$species,taxon19$sequence)
  names(taxon20)<-c("name","sequence")
  return(taxon20)
}

#########################
##################
############ APP BEGINS
##################
#########################

##################
############
#########
#####USER INTERFACE

ui <- navbarPage(title=tags$em(tags$b("BAGS: Barcode, Audit & Grade System v1.0.2")),inverse=TRUE,windowTitle="BAGS: Barcode, Audit & Grade System",
                 
                 ####HOME TAB
                 tabPanel(title="HOME", setBackgroundColor(
                   color = c("#e6f9ff", "#7aa8b8"),
                   gradient = "linear",
                   direction = "top"
                 ),
                 
                 fluidRow(column(5,align="left",tags$span(style="color:#803300",(tags$h1(tags$em(tags$strong("BAGS: Barcode, Audit & Grade System")))))),
                          column(2,div(style="display:inline-block",tags$a(href="https://www.uminho.pt/PT",tags$img(src='https://i.ibb.co/RccHm87/uminhoooo.png', width = "141px", height ="64px"),target="_blank"))),column(2,div(style="display:inline-block",tags$a(href="https://www.ntnu.edu/",tags$img(src='https://i.ibb.co/W0DMn9B/Logo-Ntnu-svg.png', width = "70px", height ="70px"),target="_blank"))),column(3,div(style="display:inline-block",tags$a(href="https://cbma.uminho.pt/",tags$img(src='https://i.ibb.co/1nnPWrX/Imagem2.png', width = "197px", height ="70"),target="_blank")))),tags$br(),
                 tabsetPanel(tabPanel(tags$span(style="color:#19194d",tags$h5(tags$b("OVERVIEW"))),tags$br(),fluidRow(
                   column(7,
                          tags$div(style="text-align:justify",column(12,align="center",tags$h3(tags$strong(tags$u("OVERVIEW")))),
                                   tags$h4(tags$p("BAGS emerged as a response to the growing awareness of the susceptibility of DNA barcode reference libraries to several types of errors and inconsistencies. 
                 These can arise at various stages of the barcoding pipeline, from the collection and identification of the specimen, through the DNA sequencing and subsequent 
                 uploading of the data to DNA sequence repositories, thus becoming potential liabilities for scientific studies which use DNA barcodes as their basis, such as metabarcoding.",tags$br(),tags$br(),
                                                  "BAGS enables the user to generate reference libraries which point out incongruencies between the species names and the sequences clustered in BINs, optimizing the process of selecting the most reliable specimen records and species to work with, according to the available data. This application is also meant to facilitate revision and curation of the reference libraries. 
                 Indeed, we encourage users of BAGS and of publicly available reference libraries, to retribute to the communitty by either contributing DNA barcodes to further expand the libraries or to review and curate data.",
                                                  tags$br(),tags$br(),"Given one or more taxonomic groups present in the",tags$a(href="http://www.boldsystems.org/","BOLD database,",target="_blank"),"or a user-provided species list (in the form of a tsv file), BAGS mines and subsequently performs post-barcoding auditing and annotation 
                 of a DNA barcode library of COI-5P sequences in an automated way.  BAGS features the following tools and options:"))),
                          tags$div(style="text-align:justify",tags$br(),tags$h4(tags$ol(
                            tags$li("User library selection - taxa search or user-provided species list."),tags$br(),
                            tags$li("Library compilation - application of quality filters to the sequences and specimen data"),tags$br(),
                            tags$li("Optional marine taxa selection/exclusion filter through the ",tags$a(href="http://www.marinespecies.org/","WoRMS database.",target="_blank")),tags$br(),
                            tags$li("Auditing and annotation - implementation of the grade ranking system."),tags$br(),
                            tags$li("Output and annotation-based file sorting - fasta compilation according to grades and auditing report."))),tags$br(),tags$br(),tags$br(),tags$br()
                          )),column(5,align="center",tags$br(),
                                    tags$img(src='https://i.ibb.co/6W5hN4s/global-scheme.png',width="433px",height="881"))),tags$br(),tags$br()),
                   tabPanel(tags$span(style="color:#19194d",tags$h5(tags$b("WORKFLOW"))), fluidRow(column(1,align="left"),column(10,align="center",tags$u(tags$h3(tags$b("WORKFLOW"))),tags$br(),
                                                                                                                                 tags$div(style="text-align:justify",tags$h4(tags$p("Firstly, to use the application you should make sure the taxonomic group or groups you want to annotate are present at the",
                                                                                                                                                                                    tags$a(href="http://www.boldsystems.org/","BOLD Systems database,",target="_blank"),
                                                                                                                                                                                    "considering that intermediate taxa are usually the most likely to be absent. Additionally, the spelling of the taxa should be identical to the spelling according to the information at BOLD.",
                                                                                                                                                                                    tags$br(), tags$br(),
                                                                                                                                                                                    "The app has three main options: download a tsv library for every species belonging to the taxa,
                                                    download a tsv library for only marine species belonging to the taxa or download a tsv library for only non-marine species belonging to the taxa.
                                                    ",tags$br(),tags$br(),"Then, after you enter the name of the taxonomic group/groups or enter a file with a list of species names, a data set
                                                     will be created and curated following these steps:", tags$br(),tags$br(),tags$div(style="text-align:justify",tags$ol(
                                                       tags$li("Downloading the data set in tsv file format, consisting of specimen data and its respective COI-5P sequence belonging to the chosen taxa, from the",
                                                               tags$a(href="http://boldsystems.org/index.php/resources/api", "BOLD Public Data Portal.", target="_blank")),tags$br(),
                                                       tags$li("Filtering", tags$strong("out"), "the following from the data set:",tags$br(),tags$br(),tags$ul(tags$li("Specimens with sequences of length below the threshold chosen by the user"), tags$br(),
                                                                                                                                                               tags$li("Specimens without data on species name, BIN, lattitude or country of origin"),tags$br(),
                                                                                                                                                               tags$li("Ambiguous characters occasionally present in the species name and COI-5P sequences"),tags$br(),
                                                                                                                                                               tags$li("Specimens with sequences consisting of > 1% Ns, which are usually the most commonambiguous character"))),tags$br(), 
                                                       tags$li("In the case of the marine or non-marine taxa options, the data set is filtered once again, retaining or excluding only the species known to be from marine or brackish habitats, using their species name as reference.
                           This is achieved using the", tags$a(href="http://www.marinespecies.org/","WoRMS database,",target="_blank"),"therefore, the download and annotation will take longer."), tags$br(),
                                                       tags$li("Lastly, according to the quality and availability of the data of each specimen, qualitative grades from A-E are assigned to each species present in the data set.
                           Then, several reference libraries in fasta format are created, which can be downloaded individually."))))))),column(1,align="right"))),
                   tabPanel(tags$span(style="color:#19194d",tags$h5(tags$b("GRADES"))),fluidRow(column(12,align="center",tags$u(tags$h3(tags$b("GRADES"))))),tags$br(),
                            
                            fluidRow(column(1,align="left"),column(10,align="center",tags$div(style="text-align:justify",tags$h4(tags$p("The assignment of each grade is based on the quality, availability and replicability of the data and metadata for each species, as well as
                                                     the quality and congruence of the COI-5P sequences, evaluated in accordance to their ",
                                                                                                                                        tags$a(href="http://www.boldsystems.org/index.php/Public_BarcodeIndexNumber_Home", "Barcode Index Number (BIN).", target="_blank"),
                                                                                                                                        column(1,align="right"),tags$br(),
                                                                                                                                        "The BIN System is an online framework at",tags$a(href="http://www.boldsystems.org/", "BOLD", target="_blank"),"that generates Operational Taxonomic Units (OTUS)
                    by clustering barcode sequences algorithmically, grouping them in a manner that ideally, mirrors their respective specimen morphological identification."),tags$br(),
                                                                                                                                 tags$p("The grades are attributed to each species according to the following criteria:"),tags$br(),
                                                                                                                                 tags$ul(tags$li(tags$strong(tags$span(style="color:#1a1aff","Grade A "),"Consolidated concordance"),"The morphospecies is assigned a unique BIN, which is also assigned
uniquely to that species, plus the species has more than 10 specimens present in the library"), tags$br(),
                                                                                                                                         tags$li(tags$strong(tags$span(style="color:#00b33c","Grade B "),"Basal concordance"),"The morphospecies is assigned a unique BIN, which is also assigned
uniquely to that species, plus the species has 10 or less specimens present in the reference library"),tags$br(),
                                                                                                                                         tags$li(tags$strong(tags$span(style="color:#cca300","Grade C "),"Multiple BINs"),"The morphospecies is assigned more than one different
BINs, but each of those BINs are assigned exclusively to that species"),tags$br(),
                                                                                                                                         tags$li(tags$strong(tags$span(style="color:#ff6600","Grade D "),"Insufficient data"),"Species is not assigned discordantly, but it has less than 3 specimens available in the reference library "),tags$br(),
                                                                                                                                         tags$li(tags$strong(tags$span(style="color:#cc0000","Grade E "),"Discordant species assignment"),"Species assigned to a BIN that is
assigned to more than one different species. The specimen may match with a different
species or display paraphyly or polyphyly"),tags$br(),tags$br()))),div(style="display:inline-block",tags$img(src='https://i.ibb.co/WpxmPd4/scheme-grades3a2aece1f92bb98a.png', width = "808px", height ="499px")),tags$br(),
                                                                   tags$br(),tags$br(),tags$p(tags$strong(("The grades were adapted from the following studies:"))),
                                                                   tags$div(style="text-align:justify",tags$ul(tags$li("Costa, Filipe O., Landi, M., Martins, R., Costa, M. H., Costa, M. E., Carneiro, M., . Carvalho, G. R. (2012). A ranking system for 
                          reference libraries of DNA barcodes: application to marine fish species from Portugal. PloS One, 7(4), 1-9. doi: 10.1371/journal.pone.0035858"),
                                                                                                               tags$br(),
                                                                                                               tags$li("Oliveira, L. M., Knebelsberger, T., Landi, M., Soares, P., Raupach, M. J., & Costa, F. O. (2016). Assembling and auditing a 
                          comprehensive DNA barcode reference library for European marine fishes. Journal of Fish Biology, 89(6), 2741-2754. doi: 10.1111/jfb.13169")))))))),
                 
                 ########### DOWNLOAD AND AUDIT DATA SETS TAB
                 
                 tabPanel(title="TAXA FOR AUDITING",
                          #NORMAL
                          tabsetPanel(type="tabs",tabPanel(tags$span(style="color:#19194d",tags$h4(tags$b("TAXA SELECTION"))),tags$br(),tabsetPanel(type="pills", 
                                                                                                                                                    tabPanel(tags$span(style="color:#262626",tags$h5(tags$b("Any biome"))),tags$br(),column(12,align="center",tags$span(style="color:#000000", tags$h3(align="center",tags$em(tags$strong(tags$u("Download, audit and annotate library for all species"))))),tags$br(),
                                                                                                                                                                                                                                            textInputIcon(inputId="taxa2",icon=icon("search"),width="500px",label=tags$h5(tags$strong("Enter the name of the taxonomic group or groups separated by commas, without spaces:")),placeholder="Example: Carnivora,Ursidae,Artiodactyla,Soricomorpha"),
                                                                                                                                                                                                                                            downloadButton("downloadData_2","Download"))),
                                                                                                                                                    #MARINE  
                                                                                                                                                    tabPanel(tags$span(style="color:#262626",tags$h5(tags$b("Marine Taxa Only"))),tags$br(),column(12,align="center",tags$span(style="color:#000000", tags$h3(align="center",tags$em(tags$strong(tags$u("Download, audit and annotate library for marine species"))))),tags$br(),
                                                                                                                                                                                                                                                   textInputIcon(inputId="taxa",icon=icon("search"),width="500px",
                                                                                                                                                                                                                                                                  label=tags$h5(tags$strong("Enter the name of the taxonomic group or groups separated by commas, without spaces:")),placeholder="Example: Cetacea,Hippocampus,Octopoda"),tags$span(style="color:#b94646",tags$h6(tags$b("NOTE: This option selects only species that are considered as marine and/or brackish at",tags$a(href="http://www.marinespecies.org/","WoRMS.",target="_blank")))),
                                                                                                                                                                                                                                                   downloadButton("downloadData","Download"))),
                                                                                                                                                    
                                                                                                                                                    #NON-MARINE  
                                                                                                                                                    tabPanel(tags$span(style="color:#262626",tags$h5(tags$b("Excluding Marine Taxa"))),tags$br(),column(12,align="center",tags$span(style="color:#000000", tags$h3(align="center",tags$em(tags$strong(tags$u("Download, audit and annotate library for non-marine species"))))),tags$br(),
                                                                                                                                                                                                                                                        textInputIcon(inputId="taxa_non",icon=icon("search"),width="500px",
                                                                                                                                                                                                                                                                       label=tags$h5(tags$strong("Enter the name of the taxonomic group or groups separated by commas, without spaces:")),placeholder="Example: Palaemonidae,Salmoniformes"),tags$span(style="color:#b94646",tags$h6(tags$b("NOTE: This option excludes all species which are assigned exclusively to marine and/or brackish at",tags$a(href="http://www.marinespecies.org/","WoRMS.",target="_blank")))),
                                                                                                                                                                                                                                                        downloadButton("downloadData_non","Download"))))),
                                      #CHECKLIST
                                      tabPanel(tags$span(style="color:#19194d",tags$h4(tags$b("UPLOAD SPECIES LIST"))),
                                               column(12,align="center", tags$br(),
                                                      tags$h3(align="center",tags$em(tags$strong(tags$u(tags$span(style="color:#000000","Download, audit and annotate library for species list"))))), tags$br(),
                                                      fileInput("file1", "Upload a txt or tsv file comprising a list of species:",buttonLabel = "Browse for a file",
                                                                multiple = FALSE, width="350px",
                                                                accept = c("text/tsv",
                                                                           "text/tab-separated-values,text/plain",
                                                                           ".tsv",".txt")),
                                                      tags$span(style="color:#b94646",tags$h6(tags$b("NOTE: Untick the header checkbox if your file does not have a header for the species column."))),
                                                      checkboxInput("header", "Header", TRUE),
                                                      radioButtons("disp", "Display",
                                                                   choices = c(Head = "head",
                                                                               All = "all"),
                                                                   selected = "head"),tableOutput("contents"),
                                                      downloadButton("download_1","Download"))),
                                      fluidRow(column(12, align="center",tags$br(),
                                                      sliderInput("seqsize", "Minimum sequence size in base pairs:",min = 300, max = 650, value = 500),
                                                      checkboxInput("rmv", "Remove records without data on country of origin or latitude", TRUE))),
                                      fluidRow(column(12,align="center",tags$br(),tags$br(),
                                                      tags$span(style="color:#990000", tags$h6(tags$b("NOTE: Since the download process includes the auditing and annotation of the library, the report is ready once the download is concluded."))),
                                                      tags$span(style="color:#990000", tags$h6(tags$b("Make sure to refresh the page every time you are about to download a new library."))))))),
                 
                 
                 ############ DOWNLOAD REFERENCE LIBRARIES IN FASTA FORMAT TAB
                 tabPanel(title="DOWNLOAD GRADED LIBRARIES",
                          
                          fluidRow(column(12,align="center",
                                          #INDIVIDUAL
                                          tabsetPanel(type="tabs",tabPanel(tags$span(style="color:#19194d",tags$h4(tags$b("INDIVIDUAL"))), tags$br(),
                                                                           tags$span(style="color:#000000", tags$h3(align="center",
                                                                                                                    tags$u(tags$em(tags$strong("Download graded libraries in fasta format"))))), tags$br(),
                                                                           tags$h4(tags$strong("Choose which grades to include in your library:")),
                                                                           tags$span(style="color:#b94646", tags$h6(tags$b("NOTE: Make sure the data set download is already completed"))),
                                                                           tags$br(),
                                                                           
                                                                           tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#468bb9","grade A"))),
                                                                           tags$style(type="text/css", "#downloadData2 {background-color:white;color:#0066ff}"),
                                                                           downloadButton('downloadData2', 'Download A library',class = "butt1"), tags$br(),tags$br(),
                                                                           
                                                                           tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#4d804d","grade B"))),
                                                                           tags$style(type="text/css", "#downloadData21 {background-color:white;color:#009900}"),
                                                                           downloadButton('downloadData21','Download B library',class = "butt2"),tags$br(),tags$br(),
                                                                           
                                                                           tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#ccad33","grade C"))),
                                                                           tags$style(type="text/css", "#downloadData22 {background-color:white;color:#ccad33}"),
                                                                           downloadButton('downloadData22','Download C library',class="butt3"),tags$br(),tags$br(),
                                                                           
                                                                           tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#ff6600","grade D"))),
                                                                           tags$style(type="text/css", "#downloadData23 {background-color:white;color:#ff6600}"),
                                                                           downloadButton('downloadData23','Download D library',class="butt4"),tags$br(),tags$br(),
                                                                           
                                                                           tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#cc0000","grade E"))),
                                                                           tags$style(type="text/css", "#downloadData24 {background-color:white;color:#ff0000}"),
                                                                           downloadButton('downloadData24','Download E library',class="butt5"),tags$br(),tags$br(),tags$br(),tags$br()),
                                                      #GROUPED
                                                      tabPanel(tags$span(style="color:#19194d",tags$h4(tags$b("GROUPED"))), tags$br(),
                                                               tags$span(style="color:#000000", tags$h3(align="center",
                                                                                                        tags$u(tags$em(tags$strong("Download graded libraries in fasta format"))))), tags$br(),
                                                               tags$h4(tags$strong("Choose which grades to include in your library:")),
                                                               tags$span(style="color:#b94646", tags$h6(tags$b("NOTE: Make sure the data set download is already completed"))),
                                                               tags$br(),                                      tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#4d804d","grades A and B"))),
                                                               tags$style(type="text/css", "#downloadData3 {background-color:white;color:#4d804d}"),
                                                               downloadButton('downloadData3','Download AB library',class="butt6"),
                                                               tags$br(),
                                                               tags$br(),
                                                               tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#ccad33","grades A, B and C"))),
                                                               tags$style(type="text/css", "#downloadData4 {background-color:white;color:#ccad33}"),
                                                               downloadButton('downloadData4','Download ABC library',class="butt7"),
                                                               tags$br(),
                                                               tags$br(),
                                                               
                                                               tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#ff6600","grades A, B, C and D"))),
                                                               tags$style(type="text/css", "#downloadData6 {background-color:white;color:#ff6600}"),
                                                               downloadButton('downloadData6','Download ABCD library',class="butt8"),
                                                               tags$br(),
                                                               tags$br(),
                                                               
                                                               tags$h5(tags$strong("Graded library including species with",tags$span(style="color:#cc0000","all grades assigned"))),
                                                               tags$style(type="text/css", "#downloadData5 {background-color:white;color:#cc0000}"),
                                                               downloadButton('downloadData5','Download ABCDE library',class="butt9"),
                                                               tags$br(),
                                                               tags$br()))))),
                 
                 ########## SUMMARY UI 
                 tabPanel(title="AUDITING REPORT",
                          fluidRow(tags$span(style="color:#000000", tags$h3(align="center",
                                                                            tags$em(tags$u(tags$strong("Library auditing report")))))), tags$br(),
                          fluidRow(column(12,align="center",actionBttn("clicks3",size="sm",icon =icon("arrow-circle-right"),label="REPORT"))),tags$br(), fluidRow(column(12,align="center",tags$h4(uiOutput("summary")))),
                          fluidRow(tags$span(style="color:#000000", tags$h3(align="center",
                                                                            tags$em(tags$u(tags$strong("Barplots display"))))),tags$br(),tags$br()),
                          fluidRow(column(6,align="center",
                                          actionBttn("clicks",size="sm",icon=icon("arrow-circle-right"),label="NUMBER OF SPECIMENS PER GRADE")),
                                   column(6,align="center",actionBttn("clicks2",size="sm",icon=icon("arrow-circle-right"),label="NUMBER OF SPECIES PER GRADE"))),
                          fluidRow(column(6,plotOutput("bar1")),
                                   column(6,plotOutput("bar2")))),
                 
                 #CONTACTS/RESOURCES
                 tabPanel(title="CONTACT AND RESOURCES", fluidRow(column(12,align="center",
                                                                         tags$h3(tags$strong("Citing:")),tags$h4("JoÃ£o Tadeu Fontes, Pedro Vieira, TorbjÃ¸rn Ekrem, Pedro Soares, Filipe O Costa"),
                                                                         tags$h4(tags$a(href="https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13262","BAGS: An automated Barcode, Audit & Grade System for DNA barcode reference libraries",target="_blank")),tags$br(),
                                                                         tags$h3(tags$strong("Useful links:")),
                                                                         tags$h4(tags$a(href="http://www.boldsystems.org/", "BOLD", target="_blank")),
                                                                         tags$h4(tags$a(href="http://marinespecies.org/", "WoRMS", target="_blank")),
                                                                         tags$h4(tags$a(href="https://ibol.org/", "iBOL", target="_blank")),
                                                                         tags$h4(tags$a(href="https://dnaqua.net/", "DNAqua-Net", target="_blank")),
                                                                         tags$h4(tags$a(href="https://cbma.uminho.pt/", "CBMA", target="_blank")),
                                                                         tags$h4(tags$a(href="http://ib-s.uminho.pt/", "IB-S", target="_blank")),
                                                                         tags$h4(tags$a(href="https://www.researchgate.net/lab/ME-Barcode-Molecular-Ecology-Biodiversity-and-DNA-barcoding-Filipe-O-Costa", "ME-Barcode", target="_blank")),tags$br(),tags$br(),tags$br(),
                                                                         fluidRow(column(12,align="center",icon("github","fa-5x"),tags$strong(tags$h3(tags$a(href="https://github.com/tadeu95/BAGs","GitHub repository",target="_blank"))))))),tags$br(),tags$br(),tags$br(),
                          fluidRow(column(1,align="left"),column(10,align="center",  tags$div(style="text-align:justify",tags$h5(tags$strong("Disclaimer")),
                                                                                              tags$h5(tags$p("Despite the fact that utmost care has been taken by us to guarantee the effectivness and reliability of the web application,
                          the use of the application is without any kind of warranty, expressed or implied. In no event shall the authors be liable for any damages of any type.")))),column(1,align="right"))))


##########################
##################
#########
##### APP SERVER 

server <- function(input, output){
  output$contents <- renderTable({
    
    
    req(input$file1)
    
    tryCatch(
      {
        checklist <<- read.delim(input$file1$datapath,
                                 header = input$header)
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
    
    if(input$disp == "head") {
      return(head(checklist))
    }
    else {
      return(checklist)
    }
    
  })
  taxaInput_10 <- reactive({grades_checklist(checklist,as.integer(input$seqsize),input$rmv)})
  output$download_1<- downloadHandler(
    filename = function() {
      paste("Library_checklist",".tsv")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading and annotating library for species list"), detail='This may take several minutes',
        value=10,
        {
          shiny::incProgress(10/10)
          write_tsv(taxaInput_10(), file)
        }
      )
    }
  )
  ##################### DOWNLOAD GENERAL TSV
  taxaInput_2 <- reactive({grades2(unlist(strsplit(input$taxa2, ",")),as.integer(input$seqsize),input$rmv)})
  output$downloadData_2 <- downloadHandler(
    filename = function() {
      paste(to_upper_camel_case(input$taxa2,sep_out=","), ".tsv")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading and annotating library for ",to_upper_camel_case(input$taxa2,sep_out=",")), detail='This may take several minutes',
        value=10,
        {
          shiny::incProgress(10/10)
          write_tsv(taxaInput_2(), file)
        }
      )
    }
  )
  
  #################### DOWNLOAD MARINE TSV   
  taxaInput <- reactive({grades(unlist(strsplit(input$taxa, ",")),as.integer(input$seqsize),input$rmv)})
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(to_upper_camel_case(input$taxa,sep_out=","), ".tsv")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading and annotating library for marine ",to_upper_camel_case(input$taxa,sep_out=",")), detail='This may take several minutes',
        value=10,
        {
          shiny::incProgress(10/10)
          write_tsv(taxaInput(), file)
        }
      )
    }
  )
  #################### DOWNLOAD NON-MARINE TSV   
  taxaInput_non <- reactive({grades_nonmarine(unlist(strsplit(input$taxa_non, ",")),as.integer(input$seqsize),input$rmv)})
  output$downloadData_non <- downloadHandler(
    filename = function() {
      paste(to_upper_camel_case(input$taxa_non,sep_out=","), ".tsv")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading and annotating library for non-marine ",to_upper_camel_case(input$taxa_non,sep_out=",")), detail='This may take several minutes',
        value=10,
        {
          shiny::incProgress(10/10)
          write_tsv(taxaInput_non(), file)
        }
      )
    }
  )
  ####### REPORT EVENTS
  summary_reac<-eventReactive(input$clicks3,{
    
    tagList(  
      
      tags$span(style="color:#2e2e1f",tags$h4(tags$strong(paste0("Taxa name: ",to_upper_camel_case(input$taxa2,sep_out=",")))),
                tags$h4(tags$strong(paste0("Marine taxa name: ",to_upper_camel_case(input$taxa,sep_out=",")))),
                tags$h4(tags$strong(paste0("Non-marine taxa name: ",to_upper_camel_case(input$taxa_non,sep_out=","))))),tags$br(),
      tags$p(tags$strong("Number of species:"),length(unique(taxon19$species))),
      tags$p(tags$strong("Number of BINs:"),length(unique(taxon19$BIN))),
      tags$p(tags$strong("Total Number of specimens in reference library:"),length(taxon19$species)),
      tags$br(),
      tags$h4(tags$strong("Specimens")),
      tags$p(tags$strong("Number of specimens with", tags$span(style="color:#1a1aff","A") ,"grade:"),as.numeric(length(taxon19$species[taxon19$grade=="A"]))," 
               ",tags$strong("| Percentage:"),signif((as.numeric(length(taxon19$species[taxon19$grade=="A"])*100)/length(taxon19$species)),digits=3),"%"),
      tags$p(tags$strong("Number of specimens with",  tags$span(style="color:#00b33c","B"), "grade:"),as.numeric(length(taxon19$species[taxon19$grade=="B"]))," ",tags$strong("| Percentage:")
             ,signif((as.numeric(length(taxon19$species[taxon19$grade=="B"])*100)/length(taxon19$species)),digits=3),"%"),
      tags$p(tags$strong("Number of specimens with" ,tags$span(style="color:#cca300","C"), "grade:"),as.numeric(length(taxon19$species[taxon19$grade=="C"]))," ",tags$strong("| Percentage:"),
             signif((as.numeric(length(taxon19$species[taxon19$grade=="C"])*100)/length(taxon19$species)),digits=3),"%"),
      tags$p(tags$strong("Number of specimens with", tags$span(style="color:#ff6600","D"), "grade:"),as.numeric(length(taxon19$species[taxon19$grade=="D"]))," ",tags$strong("| Percentage:"),
             signif((as.numeric(length(taxon19$species[taxon19$grade=="D"])*100)/length(taxon19$species)),digits=3),"%"),
      tags$p(tags$strong("Number of specimens with" ,tags$span(style="color:#cc0000","E") ,"grade:"),as.numeric(length(taxon19$species[taxon19$grade=="E"]))," ",tags$strong("| Percentage:"),
             signif((as.numeric(length(taxon19$species[taxon19$grade=="E"])*100)/length(taxon19$species)),digits=3),"%"),
      
      tags$br(),
      tags$br(),
      tags$h4(tags$strong("Species")),
      tags$p(tags$strong("Number of species with", tags$span(style="color:#1a1aff","A") ,"grade:"),as.numeric(length(unique(taxon19$species[taxon19$grade=="A"])))," 
               ",tags$strong("| Percentage:"),signif((as.numeric(length(unique(taxon19$species[taxon19$grade=="A"]))*100)/length(unique(taxon19$species))),digits=3),"%"),
      tags$p(tags$strong("Number of species with",  tags$span(style="color:#00b33c","B"), "grade:"),as.numeric(length(unique(taxon19$species[taxon19$grade=="B"])))," ",tags$strong("| Percentage:")
             ,signif((as.numeric(length(unique(taxon19$species[taxon19$grade=="B"]))*100)/length(unique(taxon19$species))),digits=3),"%"),
      tags$p(tags$strong("Number of species with" ,tags$span(style="color:#cca300","C"), "grade:"),as.numeric(length(unique(taxon19$species[taxon19$grade=="C"])))," ",tags$strong("| Percentage:"),
             signif((as.numeric(length(unique(taxon19$species[taxon19$grade=="C"]))*100)/length(unique(taxon19$species))),digits=3),"%"),
      tags$p(tags$strong("Number of species with", tags$span(style="color:#ff6600","D"), "grade:"),as.numeric(length(unique(taxon19$species[taxon19$grade=="D"])))," ",tags$strong("| Percentage:"),
             signif((as.numeric(length(unique(taxon19$species[taxon19$grade=="D"]))*100)/length(unique(taxon19$species))),digits=3),"%"),
      tags$p(tags$strong("Number of species with" ,tags$span(style="color:#cc0000","E") ,"grade:"),as.numeric(length(unique(taxon19$species[taxon19$grade=="E"])))," ",tags$strong("| Percentage:"),
             signif((as.numeric(length(unique(taxon19$species[taxon19$grade=="E"]))*100)/length(unique(taxon19$species))),digits=3),"%")
      
      
      
    )
  })
  output$summary <- renderUI({ summary_reac() })
  
  barplot1_reac<-eventReactive(input$clicks,{plot_summary_specimens(taxon19)})
  output$bar1<-renderPlot(barplot1_reac(),bg="transparent", execOnResize = TRUE)
  barplot2_reac<-eventReactive(input$clicks2,{plot_summary_species(taxon19)})
  output$bar2<-renderPlot(barplot2_reac(),bg="transparent", execOnResize = TRUE)
  
  ####### DOWNLOAD REFERENCE LIBRARIES IN FASTA FORMAT
  #INDIVIDUAL
  fastaInput <- reactive({create_A(taxon19)})
  fastaInput11 <- reactive({create_B(taxon19)})
  fastaInput12 <- reactive({create_C(taxon19)})
  fastaInput13 <- reactive({create_D(taxon19)})
  fastaInput14 <- reactive({create_E(taxon19)})
  output$downloadData2 <- downloadHandler(
    
    filename = function() {
      paste("grade_A_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grade A "),
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput(),file)
        })
    })
  
  output$downloadData21 <- downloadHandler(
    
    filename = function() {
      paste("grade_B_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grade B "),
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput11(),file)
        })
    })
  
  output$downloadData22 <- downloadHandler(
    
    filename = function() {
      paste("grade_C_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grade C "),
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput12(),file)
        })
    })
  
  output$downloadData23 <- downloadHandler(
    
    filename = function() {
      paste("grade_D_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grade D "),
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput13(),file)
        })
    })
  
  output$downloadData24 <- downloadHandler(
    
    filename = function() {
      paste("grade_E_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grade E "),
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput14(),file)
        })
    })
  #GROUPED
  fastaInput2 <- reactive({create_AB(taxon19)})
  fastaInput3 <- reactive({create_ABC(taxon19)})
  fastaInput5 <- reactive({create_ABCD(taxon19)})
  fastaInput4 <- reactive({create_fasta(taxon19)})
  
  output$downloadData3 <- downloadHandler(
    filename = function() {
      paste("grade_AB_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grades A and B "),
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput2(),file)
        })
    })
  output$downloadData4 <- downloadHandler(
    filename = function() {
      paste("grade_ABC_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grades A, B and C "),
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput3(),file)
        })
    })
  output$downloadData5 <- downloadHandler(
    filename = function() {
      paste("grade_whole_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grades A, B, C, D and E "),
        value=0,
        {
          
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput4(),file)
        })
    })
  output$downloadData6 <- downloadHandler(
    filename = function() {
      paste("grade_ABCD_library", ".fasta")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading fasta reference library with grades A, B, C and D"),
        value=0,
        {
          
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(5/10)
          dataframe2fas(fastaInput5(),file)
        })
    })
  
}

############## RUNNING THE APP
shinyApp(ui=ui,server=server)


