################
#####################
#####################################
#BAGs - Barcode, Audit & Grade system
#####################################
#####################
################

##### INSTALL NECESSARY PACKAGES
#install.packages(c("seqRFLP","bold","data.table","worms","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase"))

#RUN FROM GITHUB:
#1
#RUN:
#library(shiny)
#2
#RUN: 
#runGitHub("BAGs", "tadeu95")

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


#####IMPLEMENT RANKING SYSTEM FOR ALL TAXA
grades2<-function(groups){
  taxon<-bold_seqspec(taxon=groups, format = "tsv", marker="COI-5P")
  taxon2<-taxon[taxon$markercode=="COI-5P",]
  taxon2<-taxon2[(!(is.na(taxon2$lat)) | taxon2$country!="") & (taxon2$species_name!=""),]
  taxon2<-taxon2[!(taxon2$bin_uri == "" | is.na(taxon2$bin_uri)), ]
  taxon2$number<-str_count(taxon2$nucleotides, pattern="[A-Z]")
  taxon3<-taxon2[(taxon2$number>499),]
  np=(str_count(taxon3$nucleotides, "N")/str_count(taxon3$nucleotides, "[A-Z]"))*100
  taxon3$n_percent=np
  taxon3<-subset(taxon3,taxon3$n_percent<1)
  taxon3$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon3$nucleotides)
  taxon3$nucleotides=gsub("-","",taxon3$nucleotides)
  bins_list<-as.list(as.character(unique(taxon3$bin_uri)))
  a<-lapply(bins_list, function (x) bold_stats(bin = x))
  c<-lapply(a, '[[', c('species','count'))
  d<-unlist(c)
  d=data.frame(d)
  b<-cbind(unique(taxon3$bin_uri),d$d)
  b=data.frame(b)
  names(b)<-c("bin_uri","bold_species_per_bin")
  taxon3<-inner_join(taxon3,b)
  taxon3$bold_species_per_bin=as.numeric(taxon3$bold_species_per_bin)
  taxon3$species_name<-gsub("-", "", taxon3$species_name)
  taxon3$species_name<-gsub("sp.", "", taxon3$species_name)
  taxon3$species_name<-gsub("sp. nov", "", taxon3$species_name)
  taxon3$species_name<-gsub("cf.", "", taxon3$species_name)
  taxon3$species_name<-gsub("complex.", "", taxon3$species_name)
  taxon3$species_name<-gsub("cmplx.", "", taxon3$species_name)
  taxon3$species_name<-gsub("pr.", "", taxon3$species_name)
  taxon3$species_name<-gsub("f.", "", taxon3$species_name)
  taxon3$species_name<-gsub("nr.", "", taxon3$species_name)
  taxon3$species_name<-gsub("s.l.", "", taxon3$species_name)
  taxon3$species_name<-gsub("grp.", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z]+.+$", "", taxon3$species_name)
  taxon3$species_name<-gsub(" type", "", taxon3$species_name)
  taxon3$species_name<-gsub(" group", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+.+", "", taxon3$species_name)
  taxon3$species_name<-gsub("[a-z,A-Z]+[0-9]+.+", "", taxon3$species_name)
  taxon3$species_name<-sub("^[A-Z,a-z] ", "", taxon3$species_name)
  taxon3$species_name<-gsub("^[A-Z,a-z][A-Z,a-z] ", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z,a-z]$", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z,a-z][A-Z,a-z]$", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+.*", "", taxon3$species_name)
  taxon3$species_name<-gsub("[a-z,A-Z]+[0-9]+.*", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+", "", taxon3$species_name)
  taxon3$species_name<-gsub(" +$", "", taxon3$species_name)
  taxon6=subset(taxon3, lengths(gregexpr("\\w+", taxon3$species_name)) > 1)
  taxon7=taxon6
  taxon8<-data.frame(taxon7$species_name,taxon7$bin_uri,taxon7$nucleotides,taxon7$country,taxon7$family_name,taxon7$order_name,taxon7$class_name,taxon7$sampleid,taxon7$processid,taxon7$bold_species_per_bin)
  names(taxon8)<-c("species_name","BIN","nucleotides","country","family","order","class","sampleid","processid","bold_species_per_bin")
  num_species=table(taxon8$species_name)
  num_species=as.data.frame(num_species)
  names(num_species)=c("species","frequency_species")
  taxon8$grade=NA
  names(taxon8)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","bold_species_per_bin","grade")
  taxon9<-inner_join(taxon8,num_species)
  taxon9$base_number=str_count(taxon9$sequence, pattern="[A-Z]")
  taxon9=data.frame(taxon9$species,taxon9$BIN,taxon9$sequence,taxon9$country,taxon9$grade,taxon9$frequency_species,taxon9$base_number,taxon9$family,taxon9$order,taxon9$class,taxon9$sampleid,taxon9$processid,taxon9$bold_species_per_bin)
  names(taxon9)=c("species","BIN","sequence","country","grade","species_frequency","base_number","family","order","class","sampleid","processid","bold_species_per_bin")
  taxon10<-taxon9%>% 
    group_by(species) %>%
    summarise(occurrence = n_distinct(BIN),
              BIN = str_c(unique(BIN), collapse = ","))
  names(taxon10)<-c("species","bin_per_species","BIN")
  taxon11<-taxon9%>% 
    group_by(BIN) %>%
    summarise(occurrence = n_distinct(species),
              species = str_c(unique(species), collapse = ","))
  names(taxon11)<-c("BIN","species_per_bin","species")
  taxon16<-full_join(taxon9,taxon10,by = "species")
  taxon17<-data.frame(taxon16$species,taxon16$BIN.x,taxon16$sequence,taxon16$country,taxon16$grade,taxon16$species_frequency,taxon16$base_number,taxon16$bin_per_species,taxon16$family,taxon16$order,taxon16$class,taxon16$sampleid,taxon16$processid,taxon16$bold_species_per_bin)
  names(taxon17)<-c("species","BIN","COI_sequence","country","grade","species_frequency","base_number","BIN_per_species","family","order","class","sampleid","processid","bold_species_per_bin")
  taxon18<-full_join(taxon17,taxon11,by="BIN")
  colnames(taxon18)[colnames(taxon18)=="species.x"]<- "species"
  taxon18$species.y<-NULL
  taxon19<-taxon18 %>%
    mutate(grade = ifelse(species_per_bin>1 | bold_species_per_bin>1,"E",
                          ifelse(species_frequency<3,"D",
                                 ifelse(BIN_per_species>1 & (species_per_bin==1 | bold_species_per_bin==1),"C",
                                        ifelse(BIN_per_species==1 & (species_per_bin==1 | bold_species_per_bin==1) & species_frequency<11,"B",
                                               ifelse(BIN_per_species==1 & (species_per_bin==1 | bold_species_per_bin==1) & species_frequency>10,"A",NA)) ))))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_e = sum(grade=='E')),by='species') %>%
    mutate(grade = ifelse(sum_e>0,"E",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_b=sum(grade=='B')),by='species') %>%
    mutate(grade = ifelse(sum_b>0,"B",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_a=sum(grade=='A')),by='species') %>%
    mutate(grade = ifelse(sum_a>0,"A",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_c=sum(grade=='C')),by='species') %>%
    mutate(grade = ifelse(sum_c>0,"C",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_d=sum(grade=='D')),by='species') %>%
    mutate(grade = ifelse(sum_d>0,"D",grade))
  taxon19$sum_a=NULL
  taxon19$sum_d=NULL
  taxon19$sum_c=NULL
  taxon19$sum_b=NULL
  taxon19$sum_e=NULL
  taxon19$BIN_per_species=NULL
  taxon19$species_per_bin=NULL
  taxon19$species_frequency=NULL
  taxon19$bold_species_per_bin=NULL
  assign('taxon19',taxon19,envir=.GlobalEnv)
}


####IMPLEMENT RANKING SYSTEM FOR NON-MARINE TAXA
grades3<-function(groups){
  taxon<-bold_seqspec(taxon=groups, format = "tsv", marker="COI-5P")
  taxon2<-taxon[taxon$markercode=="COI-5P",]
  taxon2<-taxon2[(!(is.na(taxon2$lat)) | taxon2$country!="") & (taxon2$species_name!=""),]
  taxon2<-taxon2[!(taxon2$bin_uri == "" | is.na(taxon2$bin_uri)), ]
  taxon2$number<-str_count(taxon2$nucleotides, pattern="[A-Z]")
  taxon3<-taxon2[(taxon2$number>499),]
  np=(str_count(taxon3$nucleotides, "N")/str_count(taxon3$nucleotides, "[A-Z]"))*100
  taxon3$n_percent=np
  taxon3<-subset(taxon3,taxon3$n_percent<1)
  taxon3$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon3$nucleotides)
  taxon3$nucleotides=gsub("-","",taxon3$nucleotides)
  bins_list<-as.list(as.character(unique(taxon3$bin_uri)))
  a<-lapply(bins_list, function (x) bold_stats(bin = x))
  c<-lapply(a, '[[', c('species','count'))
  d<-unlist(c)
  d=data.frame(d)
  b<-cbind(unique(taxon3$bin_uri),d$d)
  b=data.frame(b)
  names(b)<-c("bin_uri","bold_species_per_bin")
  b$bin_uri=as.character(b$bin_uri)
  taxon3<-inner_join(taxon3,b)
  taxon3$bold_species_per_bin=as.numeric(taxon3$bold_species_per_bin)
  taxon3$species_name<-gsub("-", "", taxon3$species_name)
  taxon3$species_name<-gsub("sp.", "", taxon3$species_name)
  taxon3$species_name<-gsub("sp. nov", "", taxon3$species_name)
  taxon3$species_name<-gsub("cf.", "", taxon3$species_name)
  taxon3$species_name<-gsub("complex.", "", taxon3$species_name)
  taxon3$species_name<-gsub("cmplx.", "", taxon3$species_name)
  taxon3$species_name<-gsub("pr.", "", taxon3$species_name)
  taxon3$species_name<-gsub("f.", "", taxon3$species_name)
  taxon3$species_name<-gsub("nr.", "", taxon3$species_name)
  taxon3$species_name<-gsub("s.l.", "", taxon3$species_name)
  taxon3$species_name<-gsub("grp.", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z]+.+$", "", taxon3$species_name)
  taxon3$species_name<-gsub(" type", "", taxon3$species_name)
  taxon3$species_name<-gsub(" group", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+.+", "", taxon3$species_name)
  taxon3$species_name<-gsub("[a-z,A-Z]+[0-9]+.+", "", taxon3$species_name)
  taxon3$species_name<-sub("^[A-Z,a-z] ", "", taxon3$species_name)
  taxon3$species_name<-gsub("^[A-Z,a-z][A-Z,a-z] ", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z,a-z]$", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z,a-z][A-Z,a-z]$", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+.*", "", taxon3$species_name)
  taxon3$species_name<-gsub("[a-z,A-Z]+[0-9]+.*", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+", "", taxon3$species_name)
  taxon3$species_name<-gsub(" +$", "", taxon3$species_name)
  taxon3=subset(taxon3, lengths(gregexpr("\\w+", taxon3$species_name)) > 1)
  taxon4<-wormsbynames(as.character(taxon3$species_name), marine_only=FALSE, ids=TRUE)
  taxon4<-taxon4%>%
    mutate(species_name = ifelse(status=="unaccepted" & !is.na(status),valid_name,
                                 ifelse(status=="accepted" & !is.na(status),name,name)))
  taxon5<-subset(taxon4, taxon4$isMarine!="1")
  taxon6<-subset(taxon3,as.character(taxon3$species_name)%in% as.character(taxon5$species_name))
  taxon7=taxon6
  taxon8<-data.frame(taxon7$species_name,taxon7$bin_uri,taxon7$nucleotides,taxon7$country,taxon7$family_name,taxon7$order_name,taxon7$class_name,taxon7$sampleid,taxon7$processid,taxon7$bold_species_per_bin)
  names(taxon8)<-c("species_name","BIN","nucleotides","country","family","order","class","sampleid","processid","bold_species_per_bin")
  num_species=table(taxon8$species_name)
  num_species=as.data.frame(num_species)
  names(num_species)=c("species","frequency_species")
  taxon8$grade=NA
  names(taxon8)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","bold_species_per_bin","grade")
  taxon9<-inner_join(taxon8,num_species)
  taxon9$base_number=str_count(taxon9$sequence, pattern="[A-Z]")
  taxon9=data.frame(taxon9$species,taxon9$BIN,taxon9$sequence,taxon9$country,taxon9$grade,taxon9$frequency_species,taxon9$base_number,taxon9$family,taxon9$order,taxon9$class,taxon9$sampleid,taxon9$processid,taxon9$bold_species_per_bin)
  names(taxon9)=c("species","BIN","sequence","country","grade","species_frequency","base_number","family","order","class","sampleid","processid","bold_species_per_bin")
  taxon10<-taxon9%>% 
    group_by(species) %>%
    summarise(occurrence = n_distinct(BIN),
              BIN = str_c(unique(BIN), collapse = ","))
  names(taxon10)<-c("species","bin_per_species","BIN")
  taxon11<-taxon9%>% 
    group_by(BIN) %>%
    summarise(occurrence = n_distinct(species),
              species = str_c(unique(species), collapse = ","))
  names(taxon11)<-c("BIN","species_per_bin","species")
  taxon16<-full_join(taxon9,taxon10,by = "species")
  taxon17<-data.frame(taxon16$species,taxon16$BIN.x,taxon16$sequence,taxon16$country,taxon16$grade,taxon16$species_frequency,taxon16$base_number,taxon16$bin_per_species,taxon16$family,taxon16$order,taxon16$class,taxon16$sampleid,taxon16$processid,taxon16$bold_species_per_bin)
  names(taxon17)<-c("species","BIN","COI_sequence","country","grade","species_frequency","base_number","BIN_per_species","family","order","class","sampleid","processid","bold_species_per_bin")
  taxon18<-full_join(taxon17,taxon11,by="BIN")
  colnames(taxon18)[colnames(taxon18)=="species.x"]<- "species"
  taxon18$species.y<-NULL
  taxon19<-taxon18 %>%
    mutate(grade = ifelse(species_per_bin>1 | bold_species_per_bin>1,"E",
                          ifelse(species_frequency<3,"D",
                                 ifelse(BIN_per_species>1 & (species_per_bin==1 | bold_species_per_bin==1),"C",
                                        ifelse(BIN_per_species==1 & (species_per_bin==1 | bold_species_per_bin==1) & species_frequency<11,"B",
                                               ifelse(BIN_per_species==1 & (species_per_bin==1 | bold_species_per_bin==1) & species_frequency>10,"A",NA)) ))))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_e = sum(grade=='E')),by='species') %>%
    mutate(grade = ifelse(sum_e>0,"E",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_b=sum(grade=='B')),by='species') %>%
    mutate(grade = ifelse(sum_b>0,"B",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_a=sum(grade=='A')),by='species') %>%
    mutate(grade = ifelse(sum_a>0,"A",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_c=sum(grade=='C')),by='species') %>%
    mutate(grade = ifelse(sum_c>0,"C",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_d=sum(grade=='D')),by='species') %>%
    mutate(grade = ifelse(sum_d>0,"D",grade))
  taxon19$sum_a=NULL
  taxon19$sum_d=NULL
  taxon19$sum_c=NULL
  taxon19$sum_b=NULL
  taxon19$sum_e=NULL
  taxon19$BIN_per_species=NULL
  taxon19$species_per_bin=NULL
  taxon19$species_frequency=NULL
  taxon19$bold_species_per_bin=NULL
  assign('taxon19',taxon19,envir=.GlobalEnv)
}
####IMPLEMENT RANKING SYSTEM MARINE TAXA
grades<-function(groups){
  taxon<-bold_seqspec(taxon=groups, format = "tsv", marker="COI-5P")
  taxon2<-taxon[taxon$markercode=="COI-5P",]
  taxon2<-taxon2[(!(is.na(taxon2$lat)) | taxon2$country!="") & (taxon2$species_name!=""),]
  taxon2<-taxon2[!(taxon2$bin_uri == "" | is.na(taxon2$bin_uri)), ]
  taxon2$number<-str_count(taxon2$nucleotides, pattern="[A-Z]")
  taxon3<-taxon2[(taxon2$number>499),]
  np=(str_count(taxon3$nucleotides, "N")/str_count(taxon3$nucleotides, "[A-Z]"))*100
  taxon3$n_percent=np
  taxon3<-subset(taxon3,taxon3$n_percent<1)
  taxon3$nucleotides=gsub("[^ATGCNRYSWKMBDHV]+", "", taxon3$nucleotides)
  taxon3$nucleotides=gsub("-","",taxon3$nucleotides)
  bins_list<-as.list(as.character(unique(taxon3$bin_uri)))
  a<-lapply(bins_list, function (x) bold_stats(bin = x))
  c<-lapply(a, '[[', c('species','count'))
  d<-unlist(c)
  d=data.frame(d)
  b<-cbind(unique(taxon3$bin_uri),d$d)
  b=data.frame(b)
  names(b)<-c("bin_uri","bold_species_per_bin")
  b$bin_uri=as.character(b$bin_uri)
  taxon3<-inner_join(taxon3,b)
  taxon3$bold_species_per_bin=as.numeric(taxon3$bold_species_per_bin)
  taxon3$species_name<-gsub("-", "", taxon3$species_name)
  taxon3$species_name<-gsub("sp.", "", taxon3$species_name)
  taxon3$species_name<-gsub("sp. nov", "", taxon3$species_name)
  taxon3$species_name<-gsub("cf.", "", taxon3$species_name)
  taxon3$species_name<-gsub("complex.", "", taxon3$species_name)
  taxon3$species_name<-gsub("cmplx.", "", taxon3$species_name)
  taxon3$species_name<-gsub("pr.", "", taxon3$species_name)
  taxon3$species_name<-gsub("f.", "", taxon3$species_name)
  taxon3$species_name<-gsub("nr.", "", taxon3$species_name)
  taxon3$species_name<-gsub("s.l.", "", taxon3$species_name)
  taxon3$species_name<-gsub("grp.", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z]+.+$", "", taxon3$species_name)
  taxon3$species_name<-gsub(" type", "", taxon3$species_name)
  taxon3$species_name<-gsub(" group", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+.+", "", taxon3$species_name)
  taxon3$species_name<-gsub("[a-z,A-Z]+[0-9]+.+", "", taxon3$species_name)
  taxon3$species_name<-sub("^[A-Z,a-z] ", "", taxon3$species_name)
  taxon3$species_name<-gsub("^[A-Z,a-z][A-Z,a-z] ", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z,a-z]$", "", taxon3$species_name)
  taxon3$species_name<-gsub(" [A-Z,a-z][A-Z,a-z]$", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+.*", "", taxon3$species_name)
  taxon3$species_name<-gsub("[a-z,A-Z]+[0-9]+.*", "", taxon3$species_name)
  taxon3$species_name<-gsub("[0-9]+", "", taxon3$species_name)
  taxon3$species_name<-gsub(" +$", "", taxon3$species_name)
  taxon3=subset(taxon3, lengths(gregexpr("\\w+", taxon3$species_name)) > 1)
  taxon4<-wormsbynames(as.character(taxon3$species_name), marine_only=FALSE, ids=TRUE)
  taxon4<-taxon4%>%
    mutate(species_name = ifelse(status=="unaccepted" & !is.na(status),valid_name,
                                 ifelse(status=="accepted" & !is.na(status),name,name)))
  taxon5<-subset(taxon4, taxon4$isMarine=="1" | taxon4$isBrackish=="1")
  taxon6<-subset(taxon3,as.character(taxon3$species_name)%in% as.character(taxon5$species_name))
  taxon7=taxon6
  taxon8<-data.frame(taxon7$species_name,taxon7$bin_uri,taxon7$nucleotides,taxon7$country,taxon7$family_name,taxon7$order_name,taxon7$class_name,taxon7$sampleid,taxon7$processid,taxon7$bold_species_per_bin)
  names(taxon8)<-c("species_name","BIN","nucleotides","country","family","order","class","sampleid","processid","bold_species_per_bin")
  num_species=table(taxon8$species_name)
  num_species=as.data.frame(num_species)
  names(num_species)=c("species","frequency_species")
  taxon8$grade=NA
  names(taxon8)=c("species","BIN","sequence","country","family","order","class","sampleid","processid","bold_species_per_bin","grade")
  taxon9<-inner_join(taxon8,num_species)
  taxon9$base_number=str_count(taxon9$sequence, pattern="[A-Z]")
  taxon9=data.frame(taxon9$species,taxon9$BIN,taxon9$sequence,taxon9$country,taxon9$grade,taxon9$frequency_species,taxon9$base_number,taxon9$family,taxon9$order,taxon9$class,taxon9$sampleid,taxon9$processid,taxon9$bold_species_per_bin)
  names(taxon9)=c("species","BIN","sequence","country","grade","species_frequency","base_number","family","order","class","sampleid","processid","bold_species_per_bin")
  taxon10<-taxon9%>% 
    group_by(species) %>%
    summarise(occurrence = n_distinct(BIN),
              BIN = str_c(unique(BIN), collapse = ","))
  names(taxon10)<-c("species","bin_per_species","BIN")
  taxon11<-taxon9%>% 
    group_by(BIN) %>%
    summarise(occurrence = n_distinct(species),
              species = str_c(unique(species), collapse = ","))
  names(taxon11)<-c("BIN","species_per_bin","species")
  taxon16<-full_join(taxon9,taxon10,by = "species")
  taxon17<-data.frame(taxon16$species,taxon16$BIN.x,taxon16$sequence,taxon16$country,taxon16$grade,taxon16$species_frequency,taxon16$base_number,taxon16$bin_per_species,taxon16$family,taxon16$order,taxon16$class,taxon16$sampleid,taxon16$processid,taxon16$bold_species_per_bin)
  names(taxon17)<-c("species","BIN","COI_sequence","country","grade","species_frequency","base_number","BIN_per_species","family","order","class","sampleid","processid","bold_species_per_bin")
  taxon18<-full_join(taxon17,taxon11,by="BIN")
  colnames(taxon18)[colnames(taxon18)=="species.x"]<- "species"
  taxon18$species.y<-NULL
  taxon19<-taxon18 %>%
    mutate(grade = ifelse(species_per_bin>1 | bold_species_per_bin>1,"E",
                          ifelse(species_frequency<3,"D",
                                 ifelse(BIN_per_species>1 & (species_per_bin==1 | bold_species_per_bin==1),"C",
                                        ifelse(BIN_per_species==1 & (species_per_bin==1 | bold_species_per_bin==1) & species_frequency<11,"B",
                                               ifelse(BIN_per_species==1 & (species_per_bin==1 | bold_species_per_bin==1) & species_frequency>10,"A",NA)) ))))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_e = sum(grade=='E')),by='species') %>%
    mutate(grade = ifelse(sum_e>0,"E",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_b=sum(grade=='B')),by='species') %>%
    mutate(grade = ifelse(sum_b>0,"B",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_a=sum(grade=='A')),by='species') %>%
    mutate(grade = ifelse(sum_a>0,"A",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_c=sum(grade=='C')),by='species') %>%
    mutate(grade = ifelse(sum_c>0,"C",grade))
  taxon19<-taxon19 %>% left_join(taxon19 %>% 
                                   group_by(species) %>% 
                                   summarize(sum_d=sum(grade=='D')),by='species') %>%
    mutate(grade = ifelse(sum_d>0,"D",grade))
  taxon19$sum_a=NULL
  taxon19$sum_d=NULL
  taxon19$sum_c=NULL
  taxon19$sum_b=NULL
  taxon19$sum_e=NULL
  taxon19$BIN_per_species=NULL
  taxon19$species_per_bin=NULL
  taxon19$species_frequency=NULL
  taxon19$bold_species_per_bin=NULL
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
    ggtitle("Specimens per grade")+
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
    ggtitle("Species per grade")+
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
  taxon_a$name=paste(taxon_a$species,taxon_a$BIN,taxon_a$processid,sep=" | ")
  taxon_a2<-data.frame(taxon_a$name,taxon_a$COI_sequence)
  names(taxon_a2)<-c("name","sequence")
  assign('taxon_a2',taxon_a2,envir=.GlobalEnv)
}
#grade B
create_B=function(taxon19){
  taxon_b=taxon19[taxon19$grade=="B",]
  taxon_b$name=paste(taxon_b$species,taxon_b$BIN,taxon_b$processid,sep=" | ")
  taxon_b2<-data.frame(taxon_b$name,taxon_b$COI_sequence)
  names(taxon_b2)<-c("name","sequence")
  assign('taxon_b2',taxon_b2,envir=.GlobalEnv)
}
#grade C
create_C=function(taxon19){
  taxon_c=taxon19[taxon19$grade=="C",]
  taxon_c$name=paste(taxon_c$species,taxon_c$BIN,taxon_c$processid,sep=" | ")
  taxon_c2<-data.frame(taxon_c$name,taxon_c$COI_sequence)
  names(taxon_c2)<-c("name","sequence")
  assign('taxon_c2',taxon_c2,envir=.GlobalEnv)
}
#grade D
create_D=function(taxon19){
  taxon_d=taxon19[taxon19$grade=="D",]
  taxon_d$name=paste(taxon_d$species,taxon_d$BIN,taxon_d$processid,sep=" | ")
  taxon_d2<-data.frame(taxon_d$name,taxon_d$COI_sequence)
  names(taxon_d2)<-c("name","sequence")
  assign('taxon_d2',taxon_d2,envir=.GlobalEnv)
}
#grade E
create_E=function(taxon19){
  taxon_e=taxon19[taxon19$grade=="E",]
  taxon_e$name=paste(taxon_e$species,taxon_e$BIN,taxon_e$processid,sep=" | ")
  taxon_e2<-data.frame(taxon_e$name,taxon_e$COI_sequence)
  names(taxon_e2)<-c("name","sequence")
  assign('taxon_e2',taxon_e2,envir=.GlobalEnv)
}

####Create fasta files grouped

#Only grades A and B
create_AB=function(taxon19){
  taxon_ab=taxon19[taxon19$grade=="A" | taxon19$grade=="B",]
  taxon_ab$name=paste(taxon_ab$species,taxon_ab$BIN,taxon_ab$processid,sep=" | ")
  taxon_ab2<-data.frame(taxon_ab$name,taxon_ab$COI_sequence)
  names(taxon_ab2)<-c("name","sequence")
  assign('taxon_ab2',taxon_ab2,envir=.GlobalEnv)
}

#Only grades A , B and C
create_ABC=function(taxon19){
  taxon_abc=taxon19[taxon19$grade=="A" | taxon19$grade=="B" | taxon19$grade=="C",]
  taxon_abc$name=paste(taxon_abc$species,taxon_abc$BIN,taxon_abc$sprocessid,sep=" | ")
  taxon_abc2<-data.frame(taxon_abc$name,taxon_abc$COI_sequence)
  names(taxon_abc2)<-c("name","sequence")
  assign('taxon_abc2',taxon_abc2,envir=.GlobalEnv)
}
#Only grades A , B , C and D
create_ABCD=function(taxon19){
  taxon_abcd=taxon19[taxon19$grade=="A" | taxon19$grade=="B" | taxon19$grade=="C" | taxon19$grade=="D",]
  taxon_abcd$name=paste(taxon_abcd$species,taxon_abcd$BIN,taxon_abcd$processid,sep=" | ")
  taxon_abcd2<-data.frame(taxon_abcd$name,taxon_abcd$COI_sequence)
  names(taxon_abcd2)<-c("name","sequence")
  assign('taxon_abcd2',taxon_abcd2,envir=.GlobalEnv)
}

#All grades
create_fasta=function(taxon19){
  taxon19$name=paste(taxon19$species, taxon19$BIN,taxon19$processid, sep=" | ")
  taxon20<-data.frame(taxon19$name,taxon19$COI_sequence)
  names(taxon20)<-c("name","sequence")
  assign('taxon20',taxon20,envir=.GlobalEnv)
}


#########################
##################
############ APP BEGINS
##################
#####################################################################

##################
############
#########
#####USER INTERFACE

ui <- navbarPage(title=tags$em(tags$b("BAGs: Barcode, Audit & Grade system v0.1")),inverse=TRUE,windowTitle="BAGs: Barcode, Audit & Grade system",
       
####HOME TAB
                 tabPanel(title="HOME", setBackgroundColor(
    color = c("#e6f9ff", "#7aa8b8"),
    gradient = "linear",
    direction = "top"
  ),
      
  tags$span(style="color:#803300",(tags$h1(tags$em(tags$strong("BAGs: Barcode, Audit & Grade system"))))),  tags$hr(), tags$br(),
  fluidRow(
  column(9,
  tags$div(style="text-align:justify",tags$h3(tags$strong("Motivation")),
  tags$h4(tags$p("The purpose of this web application is to, given one or more taxonomic groups present at the",tags$a(href="http://www.boldsystems.org/","BOLD Systems database,",target="_blank"),"perform post barcoding auditing and annotation of a 
                                                     DNA barcode reference library
                                                     of COI-5P sequences, in a semi-automated way. Subsequently, qualitative grades from A-E are assigned to each species present in the reference library, according to the quality and availability of the data.
                                                     ",tags$br(),tags$br(),
                                                     "This comes as a response to the fact that the data present in DNA barcode reference libraries is prone to several types of errors and inconsistencies. These can arise at any
                                                    moment, from the collection and identification of the specimen, through the DNA sequencing and subsequent update of the data to biological databases, thus becoming potential liabilities for scientific
                                                    studies which use DNA barcodes as their basis. Therefore, this system allows the user to obtain the most useful and congruent reference library
                                                    for the taxa in question, optimizing the process of choosing the best specimens and species to work with.
                                                     "
                                                     )), tags$br(),
           tags$h3(tags$strong("Workflow")),
                                                    tags$h4(tags$p("Firstly, to use this app you should make sure the taxonomic group or groups you want to annotate are present at the",
                                                                   tags$a(href="http://www.boldsystems.org/","BOLD Systems database,",target="_blank"),
                                                                   "considering that intermediate taxa are usually the most likely to be absent.",
                                                      tags$br(), tags$br(),
                                                    
                                          
                                                      "The app has three main options: download a tsv library for every species belonging to the taxa;
                                                    download a tsv library for only marine species belonging to the taxa; and
                                                    download a tsv library for only non-marine species belonging to the taxa;
                                                    ",tags$br(),tags$br(),
                                                    "Then, after you enter the name of the taxonomic group in one of the four download fields, a data set
                                                     will be created and curated following these steps:", tags$br(),
                 tags$div(style="text-align:justify",tags$ol(
                   tags$li("Downloading the data set in tsv file format, consisting of specimen data and its respective COI-5P sequence belonging to the taxonomic group chosen, from the",
                           tags$a(href="http://boldsystems.org/index.php/resources/api", "BOLD Public Data Portal.", target="_blank")),tags$br(),
                   tags$li("Filtering", tags$strong("out"), "the following from the data set:",tags$br(),tags$br(),tags$ul(tags$li("Specimens with sequences of length < 500bp"), tags$br(),
                                                                                  tags$li("Specimens without data on species name, BIN, lattitude or country of origin"),tags$br(),
                                                                                  tags$li("Ambiguous characters occasionally present in the species name and COI-5P sequences"),tags$br(),
                                                                                  tags$li("Specimens with sequences consisting of > 1% Ns, which are usually the most common
                                                                                          ambiguous character"))),tags$br(), 
                   tags$li("In the case of the marine and non-marine taxa options, the data set is filtered once again, retaining only the species known to be from those habitats, using their species name as reference.
                           This is achieved using the", tags$a(href="http://www.marinespecies.org/","WoRMS database,",target="_blank"),"therefore, the download and annotation will take longer."), tags$br(),
                   tags$li("Lastly, according to the quality and availability of the data of each specimen, qualitative grades from A-E are assigned to each species present in the data set.
                           Then, several reference libraries in fasta format are created, which can be downloaded individually."))))),
  tags$div(style="text-align:justify",tags$h3(tags$strong("Disclaimer")),
           tags$h4(tags$p("Despite the fact that utmost care has been taken by us to guarantee the effectivness and reliability of the web application,
                          the use of the application is without any kind of warranty. In not event shall the authors be liable for any damages(...)"))))),column(3,align="center",tags$br(),tags$br(),
                div(style="display:inline-block",tags$img(src='http://biodiversitygenomics.net/site/wp-content/uploads/2016/01/logo_bold.png', width = "250px", height ="165px"),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  div(style="display:inline-block",tags$img(src='http://marinespecies.org/images/layout/worms_logo.png', width = "325px", height ="76px")))))),

###DEFINITION OF THE GRADES TAB
  tabPanel(title="GRADES",   setBackgroundColor(color = c("#e6f9ff", "#7aa8b8"),gradient = "linear",direction = "top"),
  fluidRow(column(12,align="center",tags$u(tags$h2(tags$em(tags$b("Definition of the grades")))))),tags$br(),
    
  fluidRow(column(2,align="left"),column(8,align="center",tags$div(style="text-align:justify",tags$h4(tags$p("The assignment of each grade is based on the quality, availability and replicability of the data and metadata for each species, as well as
                                                     the quality and congruence of the COI-5P sequences, evaluated in accordance to their ",
                                                                                                             tags$a(href="http://www.boldsystems.org/index.php/Public_BarcodeIndexNumber_Home", "Barcode Index Number (BIN).", target="_blank"),
          column(2,align="right"),tags$br(),
          "The BIN System is an online framework at",tags$a(href="http://www.boldsystems.org/", "BOLD", target="_blank"),"that generates Operational Taxonomic Units (OTUS)
                    by clustering barcode sequences algorithmically, grouping them in a manner that ideally, mirrors their respective specimen morphological identification."),tags$br(),
          tags$p("The grades are attributed to each species according to the following criteria:"),tags$br(),
             tags$ul(tags$li(tags$strong(tags$span(style="color:#1a1aff","Grade A "),"External concordance"),"Species is assigned a unique BIN, which is also assigned
uniquely to that species, plus the species has more than 10 specimens present in the reference library"), tags$br(),
                     tags$li(tags$strong(tags$span(style="color:#00b33c","Grade B "),"Internal concordance"),"Species is assigned a unique BIN, which is also assigned
uniquely to that species, plus the species has 10 or less specimens present in the reference library"),tags$br(),
                     tags$li(tags$strong(tags$span(style="color:#cca300","Grade C "),"Sub-optimal concordance"),"Species is assigned more than one different
BINs, but each of those BINs are assigned exclusively to that species"),tags$br(),
                     tags$li(tags$strong(tags$span(style="color:#ff6600","Grade D "),"Insufficient data"),"Species is not assigned discordantly, but it has less than 3 specimens available in the reference library "),tags$br(),
                     tags$li(tags$strong(tags$span(style="color:#cc0000","Grade E "),"Discordant species assignment"),"Species assigned to a BIN that is
assigned to more than one different species. The specimen may match with a different
species or display paraphyly or polyphyly"),tags$br(),tags$br()))),
          tags$p(tags$strong(("The grades were adapted from the following studies:"))),
          tags$div(style="text-align:justify",tags$ul(tags$li("Costa, Filipe O., Landi, M., Martins, R., Costa, M. H., Costa, M. E., Carneiro, M., . Carvalho, G. R. (2012). A ranking system for 
                          reference libraries of DNA barcodes: application to marine fish species from Portugal. PloS One, 7(4), 1-9. doi: 10.1371/journal.pone.0035858"),
                                                      tags$br(),
                  tags$li("Oliveira, L. M., Knebelsberger, T., Landi, M., Soares, P., Raupach, M. J., & Costa, F. O. (2016). Assembling and auditing a 
                          comprehensive DNA barcode reference library for European marine fishes. Journal of Fish Biology, 89(6), 2741-2754. doi: 10.1111/jfb.13169")))))),

########### DOWNLOAD AND AUDIT DATA SETS TAB

  tabPanel(title="SELECT TAXA FOR AUDITING",
           #NORMAL
           tabsetPanel(type="tabs",tabPanel(tags$span(style="color:#19194d",tags$h4(tags$b("ALL TAXA"))),tags$br(),column(12,align="center",tags$span(icon("globe-africa","fa-3x"),style="color:#000000", tags$h3(align="center",tags$em(tags$strong(tags$u("Download and audit data set for all species"))))),tags$br(),
               textInputAddon(inputId="taxa2",addon=icon("search"),width="500px",label=tags$h5(tags$strong("Enter the name of the taxonomic group or groups separated by commas, without spaces:")),placeholder="Example: Carnivora,Ursidae,Mustelidae"),
               downloadBttn("downloadData_2",size="sm","Download"))),
               #MARINE  
               tabPanel(tags$span(style="color:#19194d",tags$h4(tags$b("MARINE TAXA"))),tags$br(),column(12,align="center",tags$span(icon("fish","fa-3x"),style="color:#000000", tags$h3(align="center",tags$em(tags$strong(tags$u("Download and audit data set for marine species"))))),tags$br(),
                                                                                                         textInputAddon(inputId="taxa",addon=icon("search"),width="500px",
                                                                                                                        label=tags$h5(tags$strong("Enter the name of the taxonomic group or groups separated by commas, without spaces:")),placeholder="Example: Cetacea,Hippocampus,Octopoda"),tags$span(style="color:#b94646",tags$h6(tags$b("NOTE: The filtering of the marine species is fully done according to the data present at",tags$a(href="http://www.marinespecies.org/","WoRMS.",target="_blank")))),
                                                                                                         downloadBttn("downloadData",size="sm","Download"))),
            #NON MARINE
            tabPanel(tags$span(style="color:#19194d",tags$h4(tags$b("NON-MARINE TAXA"))),tags$br(),column(12,align="center",tags$span(icon("fish","fa-3x"),style="color:#000000", tags$h3(align="center",tags$em(tags$strong(tags$u("Download and audit data set for non-marine species"))))),tags$br(),textInputAddon(inputId="taxa3",addon=icon("search"),width="500px",
                label=tags$h5(tags$strong("Enter the name of the taxonomic group or groups separated by commas, without spaces:")),placeholder="Example: Palaemonidae,Salmoniformes,Nemacheilidae,Balitoridae"),tags$span(style="color:#b94646",tags$h6(tags$b("NOTE: The filtering of the non-marine aquatic species is fully done according to the data present at",tags$a(href="http://www.marinespecies.org/","WoRMS.",target="_blank")))),
                downloadBttn("downloadData_3",size="sm","Download"))),


            fluidRow(column(12,align="center",tags$br(),tags$br(),
                                            tags$span(style="color:#990000", tags$h5(tags$b("NOTE: Make sure to refresh the page every time you are about to download
                                                                                            a new data set"))))))),

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
                                                     downloadBttn('downloadData2',size="sm", 'Download A library'), tags$br(),tags$br(),
                                                     
                                                     tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#4d804d","grade B"))),
                                                     downloadBttn('downloadData21',size="sm", color="success",'Download B library'),tags$br(),tags$br(),
                                                     
                                                     tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#ccad33","grade C"))),
                                                     downloadBttn('downloadData22',size="sm", color="warning",'Download C library'),tags$br(),tags$br(),
                                                     
                                                     tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#ff6600","grade D"))),
                                                     downloadBttn('downloadData23',size="sm", color="warning",'Download D library'),tags$br(),tags$br(),
                                                     
                                                     tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#cc0000","grade E"))),
                                                     downloadBttn('downloadData24',size="sm", color="danger",'Download E library')),
                    #GROUPED
                                tabPanel(tags$span(style="color:#19194d",tags$h4(tags$b("GROUPED"))), tags$br(),
                                         tags$span(style="color:#000000", tags$h3(align="center",
                                                                                  tags$u(tags$em(tags$strong("Download graded libraries in fasta format"))))), tags$br(),
                                         tags$h4(tags$strong("Choose which grades to include in your library:")),
                                         tags$span(style="color:#b94646", tags$h6(tags$b("NOTE: Make sure the data set download is already completed"))),
                                         tags$br(),                                      tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#4d804d","grades A and B"))),
                                         downloadBttn('downloadData3', size="sm",color="success",'Download AB library'),
                                         tags$br(),
                                         tags$br(),
                                         
                                         
                                         tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#ccad33","grades A, B and C"))),
                                         downloadBttn('downloadData4', size="sm",color="warning",'Download ABC library'),
                                         tags$br(),
                                         tags$br(),
                                         
                                         tags$h5(tags$strong("Graded library including only species with",tags$span(style="color:#ff6600","grades A, B, C and D"))),
                                         downloadBttn('downloadData6', size="sm",color="warning",'Download ABCD library'),
                                         tags$br(),
                                         tags$br(),
                                         
                                         tags$h5(tags$strong("Graded library including species with",tags$span(style="color:#cc0000","all grades assigned"))),
                                         downloadBttn('downloadData5',size="sm",color="danger", 'Download ABCDE library'),
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
                  column(6,align="center",actionBttn("clicks2",size="sm",icon=icon("arrow-circle-right"),label="NUMBER OF SPECIES PER GRADE"))), tags$br(),tags$br(),
  fluidRow(column(6,plotOutput("bar1")),
           column(6,plotOutput("bar2")))),tags$br(),tags$br())



##########################
##################
#########
##### APP SERVER 

server <- function(input, output){
##################### DOWNLOAD GENERAL TSV
  taxaInput_2 <- reactive({grades2(unlist(strsplit(input$taxa2, ",")))})
  output$downloadData_2 <- downloadHandler(
    filename = function() {
      paste(to_upper_camel_case(input$taxa2,sep_out=","), ".tsv")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading and annotating dataset for ",to_upper_camel_case(input$taxa2,sep_out=",")), detail='This may take several minutes',
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(4/10)
          write_tsv(taxaInput_2(), file)
        }
      )
    }
  )
#################### DOWNLOAD NON-MARINE TSV 
  taxaInput_3 <- reactive({grades3(unlist(strsplit(input$taxa3, ",")))})
  output$downloadData_3 <- downloadHandler(
    filename = function() {
      paste(to_upper_camel_case(input$taxa3,sep_out=","), ".tsv")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading and annotating dataset for non-marine ",to_upper_camel_case(input$taxa3,sep_out=",")), detail='This may take several minutes',
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(4/10)
          write_tsv(taxaInput_3(), file)
        }
      )
    }
  )  

#################### DOWNLOAD MARINE TSV   
  taxaInput <- reactive({grades(unlist(strsplit(input$taxa, ",")))})
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(to_upper_camel_case(input$taxa,sep_out=","), ".tsv")
    },
    content = function(file) {
      shiny::withProgress(
        message=paste0("Downloading and annotating dataset for marine ",to_upper_camel_case(input$taxa,sep_out=",")), detail='This may take several minutes',
        value=0,
        {
          shiny::incProgress(1/10)
          Sys.sleep(1)
          shiny::incProgress(4/10)
          write_tsv(taxaInput(), file)
        }
      )
    }
  )

####### REPORT EVENTS
  summary_reac<-eventReactive(input$clicks3,{

    tagList(  
      
      tags$span(style="color:#2e2e1f",tags$h4(tags$strong(paste0("Taxa name: ",to_upper_camel_case(input$taxa2,sep_out=",")))),
                tags$h4(tags$strong(paste0("Non-marine taxa name: ",to_upper_camel_case(input$taxa3,sep_out=",")))),
                tags$h4(tags$strong(paste0("Marine taxa name: ",to_upper_camel_case(input$taxa,sep_out=","))))
                
                
                ),tags$br(),
      tags$p(tags$strong("Number of species:"),length(unique(taxon19$species))),
                        tags$p(tags$strong("Number of BINs:"),length(unique(taxon19$BIN))),
      tags$p(tags$strong("Total Number of specimens in reference library:"),length(taxon19$species)),
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
                   signif((as.numeric(length(unique(taxon19$species[taxon19$grade=="E"]))*100)/length(unique(taxon19$species))),digits=3),"%"),
      tags$br(),
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
             signif((as.numeric(length(taxon19$species[taxon19$grade=="E"])*100)/length(taxon19$species)),digits=3),"%")

            
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

