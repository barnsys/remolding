
library(seqinr)
library(stringr)
library(phangorn)

setwd("E:\\sequence_translation")

#трансляция стандартным кодом
  
rem<-c("Acanthogammarus_victorii", "Eulimnogammarus_vittatus", "Gmelinoides_fasciatus")

pas_dna<-read.fasta("Total_amph_COX1-1.fasta", as.string=T)
pas_names<-names(pas_dna)

data_cod_total<-read.csv("МT_genetic_code.tsv",sep="\t",header=T, stringsAsFactors=F)

n_rundom<-10

#получение набора данных без ремолдинга
wcon<-file("Total_amph_COX1-1_prot.fasta", open ="w")
for(j in 1:length(pas_names))
 {
  pas<-pas_dna[j][1]
  pas<-toupper(pas)
  n_codon<-floor(nchar(pas)/3)
  #трансляция без ремолдинга
  codons<-rep("0", n_codon)
  for(i in 1:n_codon) codons[i]<-substring(pas, 3*i-2, 3*i)
  prot_pas<-data.frame(codons=codons)
  amino_acids<-rep("0", n_codon)
  for(i in 1:n_codon)
   {
    n<-which(data_cod$codons==prot_pas$codons[i])
    amino_acids[i]<-data_cod$Standard_code[n]
   }
  writeLines(paste(">", pas_names[j], "\n", sep=""), wcon, sep="")
  writeLines(paste(amino_acids, collapse=""), wcon, sep="")
  writeLines("\n", wcon, sep="")
 }
close(wcon)
#/////////////////////////////////////////////////////


#получение набора данных с полным ремолдингом
wcon<-file("Total_amph_COX1-1_prot_fuul_remolding.fasta", open ="w")
for(j in 1:length(pas_names))
 {
  pas<-pas_dna[j][1]
  pas<-toupper(pas)
  n_codon<-floor(nchar(pas)/3)
  #трансляция без ремолдинга
  codons<-rep("0", n_codon)
  for(i in 1:n_codon) codons[i]<-substring(pas, 3*i-2, 3*i)
  prot_pas<-data.frame(codons=codons)
  amino_acids<-rep("0", n_codon)
  for(i in 1:n_codon)
   {
    n<-which(data_cod$codons==prot_pas$codons[i])
    amino_acids[i]<-data_cod$Standard_code[n]
   }
  writeLines(paste(">", pas_names[j], "\n", sep=""), wcon, sep="")
  writeLines(paste(amino_acids, collapse=""), wcon, sep="")
  writeLines("\n", wcon, sep="")
  #трансляция с полным ремолдингом
  if(length(which(rem %in% pas_names[j]))>0)
   {
    data_cod<-data_cod_total[, c("codons", "Standard_code", pas_names[j])]
	codons<-rep("0", n_codon)
    for(i in 1:n_codon) codons[i]<-substring(pas, 3*i-2, 3*i)
    prot_pas<-data.frame(codons=codons)
    amino_acids<-rep("0", n_codon)
    for(i in 1:n_codon)
     {
      n<-which(data_cod$codons==prot_pas$codons[i])
      amino_acids[i]<-data_cod[,3][n]
     }
    writeLines(paste(">", pas_names[j], "_full_remolding\n", sep=""), wcon, sep="")
    writeLines(paste(amino_acids, collapse=""), wcon, sep="")
    writeLines("\n", wcon, sep="")
   }
   
 }
close(wcon)
#/////////////////////////////////////////////////////////////////////////


#получение набора данных со случаным ремолдингом
wcon<-file("Total_amph_COX1-1_prot_rundom_remolding.fasta", open ="w")
for(j in 1:length(pas_names))
 {
  pas<-pas_dna[j][1]
  pas<-toupper(pas)
  n_codon<-floor(nchar(pas)/3)
  #трансляция без ремолдинга
  codons<-rep("0", n_codon)
  for(i in 1:n_codon) codons[i]<-substring(pas, 3*i-2, 3*i)
  prot_pas<-data.frame(codons=codons)
  amino_acids<-rep("0", n_codon)
  for(i in 1:n_codon)
   {
    n<-which(data_cod$codons==prot_pas$codons[i])
    amino_acids[i]<-data_cod$Standard_code[n]
   }
  writeLines(paste(">", pas_names[j], "\n", sep=""), wcon, sep="")
  writeLines(paste(amino_acids, collapse=""), wcon, sep="")
  writeLines("\n", wcon, sep="")
  #трансляция с полным ремолдингом
  if(length(which(rem %in% pas_names[j]))>0)
   {
    data_cod<-data_cod_total[, c("codons", "Standard_code", pas_names[j])]
	codons<-rep("0", n_codon)
    for(i in 1:n_codon) codons[i]<-substring(pas, 3*i-2, 3*i)
    prot_pas<-data.frame(codons=codons)
    amino_acids<-rep("0", n_codon)
    for(i in 1:n_codon)
     {
      n<-which(data_cod$codons==prot_pas$codons[i])
      amino_acids[i]<-data_cod[,3][n]
     }
    writeLines(paste(">", pas_names[j], "_full_remolding\n", sep=""), wcon, sep="")
    writeLines(paste(amino_acids, collapse=""), wcon, sep="")
    writeLines("\n", wcon, sep="")
	
	for(k in 1:n_rundom)
     {
      amino_acids<-rep("0", n_codon)
      for(i in 1:n_codon)
       {
        n<-which(data_cod$codons==prot_pas$codons[i])
	    p<-sample(c(2,3), 1, replace=FALSE)
        amino_acids[i]<-data_cod[n, p]
       }
	  writeLines(paste(">", pas_names[j], "_var", k, "_remolding\n", sep=""), wcon, sep="")
      writeLines(paste(amino_acids, collapse=""), wcon, sep="")
      writeLines("\n", wcon, sep="") 
	 } 
   }
 }
close(wcon)
#/////////////////////////////////////////////////////////////////////////



#оценка количества вариантов ремолдига

sp<-"Acanthogammarus_victorii"

pas<-pas_dna[[which(pas_names==sp)]][1]

pas<-toupper(pas)
n_codon<-floor(nchar(pas)/3)
n_codon

codons<-rep("0", n_codon)
for(i in 1:n_codon) codons[i]<-substring(pas, 3*i-2, 3*i)
prot_pas<-data.frame(codons=codons)

data_cod<-read.csv("МХ_геном_амфипод.txt",sep="\t",header=T, stringsAsFactors=F)

data_cod<-data_cod[, c("codons", "Standard_code", sp)]

nrep=100000
rep_amino_acids<-rep("0", nrep)

pb <- txtProgressBar(min=1, max=nrep, style = 3)
for(j in 1:nrep)
 {
  amino_acids<-rep("0", n_codon)
  for(i in 1:n_codon)
   {
    n<-which(data_cod$codons==prot_pas$codons[i])
	p<-sample(c(2,3), 1, replace=FALSE)
    amino_acids[i]<-data_cod[n, p]
   }
  rep_amino_acids[j]<-paste(amino_acids, collapse="")
  setTxtProgressBar(pb, j)
 }  
close(pb)

rep_amino_acids<-(unique(rep_amino_acids))
 
wcon<-file(paste(sp, "_CO1_prot.fasta", sep=""), open ="w")
for(i in 1:length(rep_amino_acids)) 
 {
  writeLines(paste(">var_", i, "_", sp, "\n", sep=""), wcon, sep="")
  writeLines(paste(rep_amino_acids[i], collapse=""), wcon, sep="")
  writeLines("\n", wcon, sep="")
 }
close(wcon)


