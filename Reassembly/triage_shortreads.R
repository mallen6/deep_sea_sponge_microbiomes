#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# This script requires a file containing the coverage of each contig, determined by bbmap. 
# Instructions for preparing this file are in triage_shortreads.pbs

# Test if there is at least one argument: if not, return an error
if (length(args)==0) {
	stop("At least one argument must be supplied (input file) .n", call.=FALSE)
} else if (length(args)==1) {
	# default output file
	args[2] = "out.txt"
}

print(args)
print(getwd())

out_name <- paste0(getwd(), "/collected.", args[2])
print(out_name)

library(dplyr)
require(zoo)
library(reshape)
library(data.table)


contigs <- read.table(args[1],header=FALSE,sep="\t")
library(reshape)
contigs<-rename(contigs,c(V1="name", V2="cov", V3="node", V4="length", V5="kmerCov"))  #renames the header

df <- data.frame(contigs)

long <- df[df$length > 100000, ]
sd <- sd(long$cov)
mean <- mean(long$cov)

if (is.na(sd) == TRUE) {
	sd <- 1
	print("long sd was set to one")
} else if (is.na(sd) == FALSE) {
	print("sd is ok")
}

collected <- df[df$cov <= (mean + 2*(sd)) & df$cov >= (mean - 2*(sd)) & df$length > 4999, ]


if (nrow(long) > 0){
	write.table(collected, file=out_name, sep="\t", row.names=FALSE,col.names=FALSE)
	print(paste0("long contigs were present for ", args[2]))
	} else if (nrow(long) == 0) {
	medium <- df[df$length > 10000, ]
	sd <- sd(medium$cov)
print(sd)
	mean <- mean(medium$cov)
print(mean)
	if (is.na(sd) == TRUE) {
		sd <- 1
		print("sd was set to one")
		} else if (is.na(sd) == FALSE) {
			print("sd of medium contigs is ok")
		}

	collected <- df[df$cov <= (mean + 2.5*(sd)) & df$cov >= (mean - 2*(sd)) & df$length > 4999, ]
	write.table(collected, file=out_name, sep="\t", row.names=FALSE,col.names=FALSE)
	print(paste0("medium contigs were used for ", args[2]))
	}
