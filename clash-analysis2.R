#### COMPOSITION COLLAPSED READS
directory <- paste("~/Desktop/P3 miR-9/CLASH paper/CLASH datasets mouse brain/")
setwd("~/Desktop/P3 miR-9/datasets-CLASH/")
library(data.table)
library(stringr)

fasta <- read.table(file="sra_data_SRR2736023_CTTTGGTTATCTAGCTGTA.fasta.tsv", header = FALSE, fill = FALSE, sep = "\t")
fasta$noread <- regexpr(pattern = ">SRR",fasta$V1 )
fasta <- fasta[which(fasta$noread=="-1"),]

colnames(fasta) <- c("read", "V1")
collapsed <- fasta



collapsed$read <- as.character(collapsed$read)
collapsed$mir9can <- regexpr('TCTTTGGTTATCT', collapsed$read)
collapsed$mir9alt <- regexpr('CTTTGGTTATCT', collapsed$read)

mir9canonical <- collapsed[which(collapsed$mir9can==1),]
#mir9canonical$V2 <- substr(mir9canonical$read, 26, nchar(as.character(mir9canonical$read)))
mir9alternative <- collapsed[which(collapsed$mir9alt==1),]
#mir9alternative$V2 <- substr(mir9alternative$read, 26, nchar(as.character(mir9alternative$read)))
mir9alternative <- mir9alternative[order(mir9alternative$read),]

mir9alternative$target <- substr(mir9alternative$read, 19, nchar(as.character(mir9alternative$read)))
mir9canonical$target <- substr(mir9canonical$read, 19, nchar(as.character(mir9canonical$read)))

nrow(mir9alternative)/(nrow(mir9alternative)+nrow(mir9canonical))*100

#### Analysis alternative targets

#mir9alternative$alt8mer <- regexpr('.AACCAAAA|.GACCAAAA|.AGCCAAAA|.AACCGAAA|.AACCAGAA|.AACCAAGA', mir9alternative$target) #.NNNNNNNA allow 1 G:U pair
mir9alternative$alt8mer <- regexpr('.AACCAAAA', mir9alternative$target) #.NNNNNNNA
#mir9canonical$alt8mer <- regexpr('.AACCAAAA|.GACCAAAA|.AGCCAAAA|.AACCGAAA|.AACCAGAA|.AACCAAGA', mir9canonical$target) #.NNNNNNNA allow 1 G:U pair
mir9canonical$alt8mer <- regexpr('.AACCAAAA', mir9canonical$target) #.NNNNNNNA allow 1 G:U pair

nrow(mir9alternative[which(mir9alternative$alt8mer>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$alt8mer>-1),])/nrow(mir9canonical)*100

#mir9alternative$alt7merA1 <- regexpr('.[^A,G]ACCAAAA|.[^A]GCCAAAA|.[^A]ACCGAAA|.[^A]ACCAGAA|.[^A]ACCAAGA', mir9alternative$target) # .ENNNNNNA allow 1 G:U pair
mir9alternative$alt7merA1 <- regexpr('.[^A]ACCAAAA', mir9alternative$target) # .ENNNNNNA allow 1 G:U pair
#mir9canonical$alt7merA1 <- regexpr('.[^A,G]ACCAAAA|.[^A]GCCAAAA|.[^A]ACCGAAA|.[^A]ACCAGAA|.[^A]ACCAAGA', mir9canonical$target)# .ENNNNNNA allow 1 G:U pair
mir9canonical$alt7merA1 <- regexpr('.[^A]ACCAAAA', mir9canonical$target)# .ENNNNNNA allow 1 G:U pair

nrow(mir9alternative[which(mir9alternative$alt7merA1>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$alt7merA1>-1),])/nrow(mir9canonical)*100

mir9alternative$alt7merm8 <- regexpr('.AACCAAA[^A]', mir9alternative$target) #.NNNNNNNE 
mir9canonical$alt7merm8 <- regexpr('.AACCAAA[^A]', mir9canonical$target) #.NNNNNNNE 

nrow(mir9alternative[which(mir9alternative$alt7merm8>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$alt7merm8>-1),])/nrow(mir9canonical)*100

mir9alternative$alt6mer <- regexpr('ACCAAA', mir9alternative$target)
mir9canonical$alt6mer <- regexpr('ACCAAA', mir9canonical$target)

nrow(mir9alternative[which(mir9alternative$alt6mer>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$alt6mer>-1),])/nrow(mir9canonical)*100

#### Analysis canonical targets
mir9alternative$can8mer <- regexpr('.ACCAAAGA|.GCCAAAGA|.ACCGAAGA|.ACCAGAGA|.ACCAAGGA', mir9alternative$target) #.NNNNNNNA allow 1 G:U pair
mir9canonical$can8mer <- regexpr('.ACCAAAGA|.GCCAAAGA|.ACCGAAGA|.ACCAGAGA|.ACCAAGGA', mir9canonical$target) #.NNNNNNNA allow 1 G:U pair

nrow(mir9alternative[which(mir9alternative$can8mer>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$can8mer>-1),])/nrow(mir9canonical)*100

mir9alternative$can7merA1 <- regexpr(".[^A]CCAAAGA|.[^A]CCGAAGA|.[^A]CCAGAGA|.[^A]CCAAGGA", mir9alternative$target) # .ENNNNNNA allow 1 G:U pair
mir9canonical$can7merA1 <- regexpr(".[^A]CCAAAGA|.[^A]CCGAAGA|.[^A]CCAGAGA|.[^A]CCAAGGA", mir9canonical$target) # .ENNNNNNA allow 1 G:U pair

nrow(mir9alternative[which(mir9alternative$can7merA1>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$can7merA1>-1),])/nrow(mir9canonical)*100

mir9alternative$can7merm8 <- regexpr('.ACCAAAG[^A]|.GCCAAAG[^A]|.ACCGAAG[^A]|.ACCAGAG[^A]|.ACCAAGG[^A]', mir9alternative$target) #.NNNNNNNE allow 1 G:U pair
mir9canonical$can7merm8 <- regexpr('.ACCAAAG[^A]|.GCCAAAG[^A]|.ACCGAAG[^A]|.ACCAGAG[^A]|.ACCAAGG[^A]', mir9canonical$target) #.NNNNNNNE allow 1 G:U pair

nrow(mir9alternative[which(mir9alternative$can7merm8>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$can7merm8>-1),])/nrow(mir9canonical)*100

mir9alternative$can6mer <- regexpr('.[^A,G]CCAAAG[^A]|.[^A,G]CCGAAG[^A]|.[^A,G]CCAGAG[^A]|.[^A,G]CCAAGG[^A]', mir9alternative$target) #.ENNNNNNE allow 1 G:U pair
mir9canonical$can6mer <- regexpr('.[^A,G]CCAAAG[^A]|.[^A,G]CCGAAG[^A]|.[^A,G]CCAGAG[^A]|.[^A,G]CCAAGG[^A]', mir9canonical$target) #.ENNNNNNE allow 1 G:U pair

nrow(mir9alternative[which(mir9alternative$can6mer>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$can6mer>-1),])/nrow(mir9canonical)*100

mir9alternative$canoff6mer <- regexpr(".ACCAAA[^G][^A]|.GCCAAA[^G][^A]|.ACCGAA[^G][^A]|.ACCAGA[^G][^A]|.ACCAAG[^G][^A]", mir9alternative$target) #.NNNNNNEE allow 1 G:U pair
mir9canonical$canoff6mer <- regexpr(".ACCAAA[^G][^A]|.GCCAAA[^G][^A]|.ACCGAA[^G][^A]|.ACCAGA[^G][^A]|.ACCAAG[^G][^A]", mir9canonical$target) #.NNNNNNEE allow 1 G:U pair

nrow(mir9alternative[which(mir9alternative$canoff6mer>-1),])/nrow(mir9alternative)*100
nrow(mir9canonical[which(mir9canonical$canoff6mer>-1),])/nrow(mir9canonical)*100

