setwd("~/Downloads/")
require('stringdist')
#fclip <- read.table("fCLIP-hek293-1.isomir.sequence_info.tsv", header = TRUE, sep = "\t", fill = TRUE)
fclip <- read.table("fCLIP-hek293-2.isomir.sequence_info.tsv", header = TRUE, sep = "\t", fill = TRUE)
fclip <- fclip[,c(1,2,3,10)]

fclip$SEQUENCE <- substr(fclip$SEQUENCE, 5, nchar(as.character(fclip$SEQUENCE)))
fclip$first24 <-substr(fclip$SEQUENCE, 0, 24)
fclip$dist <- stringdist(fclip$first24, "TCTTTGGTTATCTAGCTGTATGAG", method="dl")

fclip$canonical <- as.numeric(gregexpr("TCTTTGGTTATCTAGCTGTATGA", fclip$SEQUENCE, perl=TRUE))
fclip$alternative <- as.numeric(gregexpr("CTTTGGTTATCTAGCTGTATGAG", fclip$SEQUENCE, perl=TRUE))
fclip$m1alternative <- as.numeric(gregexpr("ATCTTTGGTTATCTAGCTGTATG", fclip$SEQUENCE, perl=TRUE))
fclip$m2alternative <- as.numeric(gregexpr("TATCTTTGGTTATCTAGCTGTAT", fclip$SEQUENCE, perl=TRUE))
fclip$p1alternative <- as.numeric(gregexpr("TTTGGTTATCTAGCTGTATGAGT", fclip$SEQUENCE, perl=TRUE))
fclip$sumalt <-fclip$canonical+fclip$alternative+fclip$m1alternative+fclip$m2alternative+fclip$p1alternative

#fclip$middle <- as.numeric(gregexpr("GGTTATCTAGCT", fclip$SEQUENCE, perl=TRUE))
#fclip <- fclip[which(fclip$middle!="-1"),]

fclip <- fclip[which(fclip$sumalt!="-5"),]

fclip$canonical3p <- gregexpr("ATAAAGCTAGATAACCGAAAGT", fclip$SEQUENCE, perl=TRUE)
fclip$alternative3p <- gregexpr("ATAAAGCTAGATAACCGAAAGTA", fclip$SEQUENCE, perl=TRUE)
fclip$malternative3p <- gregexpr("ATAAAGCTAGATAACCGAAAG", fclip$SEQUENCE, perl=TRUE)



#mir91 <- fclip[which(fclip$MIRNA=="hsa-miR-9-1-LOOP" & fclip$MATCH==""),]
mir91 <- fclip[which(fclip$MIRNA=="hsa-miR-9-1-LOOP"),]
ratioaltmir91 <- nrow(mir91[which(mir91$alternative=="1" & mir91$canonical=="-1"),])/(nrow(mir91[which(mir91$alternative=="1" & mir91$canonical=="-1"),])+nrow(mir91[which(mir91$alternative=="2" & mir91$canonical=="1"),]))
alternative91 <- mir91[which(mir91$alternative=="1" & mir91$canonical=="-1"),]
canonical91 <- mir91[which(mir91$alternative=="2" & mir91$canonical=="1"),]


#mir92 <- fclip[which(fclip$MIRNA=="hsa-miR-9-2-LOOP" & fclip$MATCH==""),]
mir92 <- fclip[which(fclip$MIRNA=="hsa-miR-9-2-LOOP"),]
#ratioaltmir92 <- nrow(mir92[which(mir92$alternative!="-1"),])/(nrow(mir92[which(mir92$canonical!="-1"),])+nrow(mir92[which(mir92$alternative!="-1"),]))
ratioaltmir92 <- nrow(mir92[which(mir92$alternative=="1" & mir92$canonical=="-1"),])/(nrow(mir92[which(mir92$alternative=="1" & mir92$canonical=="-1"),])+nrow(mir92[which(mir92$alternative=="2" & mir92$canonical=="1"),]))
alternative92 <- mir92[which(mir92$alternative=="1" & mir92$canonical=="-1"),]
canonical92 <- mir92[which(mir92$alternative=="2" & mir92$canonical=="1"),]

#mir93 <- fclip[which(fclip$MIRNA=="hsa-miR-9-3-LOOP" & fclip$MATCH==""),]
mir93 <- fclip[which(fclip$MIRNA=="hsa-miR-9-3-LOOP"),]
ratioaltmir93 <- nrow(mir93[which(mir93$alternative!="-1"),])/(nrow(mir93[which(mir93$canonical!="-1"),])+nrow(mir93[which(mir93$alternative!="-1"),]))

hist(as.numeric(mir91$canonical))
hist(as.numeric(mir91$alternative))

hist(as.numeric(fclip$canonical))