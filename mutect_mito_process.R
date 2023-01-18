#!/usr/bin/Rscript --no-save --no-restore
# mutect_mito_somatic.R
# workflow for processing IonTorrent mitochondrial sequencing data for somatic mutection detection


## FUNCTIONAL DEFINITIONS
is.between <- function(x, a, b) {
	# for allele frequency checks
	return((x <= b) & (x > a))
}
haplogroup.best <- function(r, id) {
	# select best hit by quality then FP and NFP counts
	x <- fread(r)
	x$SampleID <- id
	y <- x[,list(nNFP = 1/length(unlist(strsplit(Not_Found_Polys, " "))), nFP = length(unlist(strsplit(Found_Polys, " ")))), by = eval(colnames(x))]
	y$nNFP[!is.finite(y$nNFP)] <- 0
	x <- x[order(y$Quality, y$nFP, y$nNFP, decreasing = TRUE),]
	return(x[1,])
}
vcf.filter <- function(file, N = 16569) {
	y <- scan(file, comment.char = "", sep = "\n", what = character())
	pos <- sapply(y, function(x) as.numeric(unlist(strsplit(x, "\t"))[2]))
	pos <- unname(ifelse(is.na(pos), yes = 0, no = pos))
	z <- y[which(pos <= N)]
	return(z)
}

egl.retain <- function(x, n = 100) {
	# modified 09/01/22 to filter out reads with Alt < 100
	x$V6[which(is.na(x$V6))] <- 0
	x$V7[which(is.na(x$V7))] <- 0
	return(x[which(x$V7 >= n),])
}
na.set <- function(x, y) {
	x[is.na(x)] <- y
	return(x)
}
af.calc <- function(x) {
	y <- as.data.table(x)
	y <- y[,list(Pos = V2, Ref = sum(V6), Alt = sum(V7)), by = c("V2")]
	y <- y[,list(AF = (Alt + 0.01)/(Ref + Alt + 0.02)), by = c("Pos", "Ref", "Alt")]
	return(as.data.frame(y[,list(Pos, AF)]))
}


library(data.table)
## INPUT PARAMETERS

# initial filter ranges
inputFilter <- list(AF_max_run1 = 0.85, AF_min_run1 = 0.25, AF_max_run2 = 0.80, AF_min_run2 = 0.20, Z = 3, hits = 10)
inputSetup <- list(
	file_bamList = "/home/ikeda/221209cellline/221209cellline_list.txt",
	workDir = "/home/ikeda/221209cellline/",
	qsub_script = "/home/ikeda/221209cellline/mutect_mito_single.sh",  # modified to retain duplicated reads and improved EAGLE calls on 9/9/22
	ref_genome = "/home/ikeda/221209cellline/chrM_rCRS_16649.fa", # reference sequence, rCRS
	mito_symbol_annotation = "/home/ikeda/221209cellline/chrM_rCRS_16649_anno.bed")

inputSetup$output_timestamp <- "20221227"


## BATCH MUTECT & EAGLE CALLS
bamList <- read.table(inputSetup$file_bamList, header = FALSE, sep = "\t") # tab-delimited subject summary frile
bamList$file <- inputSetup$workDir
bamList$inputTumor  <- normalizePath(paste(bamList$V2, bamList$V1, sep = "/"))
bamList$outputVCF   <- paste(normalizePath(inputSetup$workDir), bamList$V3, sep = "/")
bamList$outputEagle <- gsub("\\.vcf$", ".egl.txt", bamList$outputVCF)
bamList$outputAA    <- gsub("\\.vcf$", ".aa", bamList$outputVCF)
bamList$outputID    <- gsub("\\.vcf$", "", bamList$V3)
bamList$outputSnpEff <- gsub("\\.vcf$", "_snpeff.vcf", bamList$outputVCF)
bamList$outputVcfMT <- gsub("\\.vcf$", "_mt.vcf", bamList$outputVCF)
bamList$haplogrepList <- paste(bamList$outputVCF, ".hpl", sep = "")
bamList$cmd <- paste("qsub", inputSetup$qsub_script, bamList$inputTumor, bamList$outputVCF, inputSetup$ref_genome, bamList$outputEagle) # writing qsub calls
bamList$vcfExists <- file.exists(bamList$outputVCF)
if (prod(file.exists(bamList$inputTumor)) < 1) {
	stop("the following files cannot be found:\n", paste(bamList$inputTumor[!file.exists(bamList$inputTumor)], collapse = "\n"))
	q()
}

setwd(inputSetup$workDir)
stats_mutect_calls <- sapply(bamList$cmd[!bamList$vcfExists], system)

## BATCH SNPEFF CALLS
cmd_snpEff <- paste("snpEff GRCh38.p14 ", bamList$outputVcfMT, " > ", bamList$outputSnpEff, sep = "")
vcf0 <- lapply(bamList$outputVCF, scan, what = character(), sep = "\n")
vcf_mt <- vcf0
for (i in 1:length(vcf_mt)) {
	if (!file.exists(bamList$outputSnpEff[i])) {
		message("Annotating ", bamList$outputVCF[i], " by snpEff...")
		vcf_mt[[i]] <- gsub("^chrM", "MT", vcf_mt[[i]])
		cat(vcf_mt[[i]], file = bamList$outputVcfMT[i], sep = "\n")
		if (system(cmd_snpEff[i]) == 0) {
			file.remove(bamList$outputVcfMT[i])
		}
	}
}
vcf_snpEff <- lapply(bamList$outputSnpEff, read.table, sep = "\t", header = FALSE)
names(vcf_snpEff) <- bamList$outputSnpEff
vcf_annot <- vector(mode = "list", length = length(vcf0))
vcf_aa <- vector(mode = "list", length = length(vcf0))
vcf_aa_mut <- vector(mode = "list", length = length(vcf0))
vcf_cdna <- vector(mode = "list", length = length(vcf0))
vcf_cdna_mut <- vector(mode = "list", length = length(vcf0))
for (i in 1:length(vcf_snpEff)) {
	vcf_annot[[i]] <- lapply(vcf_snpEff[[i]]$V8, function(x) unlist(strsplit(x, "\\|")))
	vcf_aa[[i]] <- lapply(vcf_annot[[i]], function(x) unique(x[grep("protein_coding", x) + 3]))
	vcf_cdna[[i]] <- lapply(vcf_annot[[i]], function(x) unique(x[grep("protein_coding", x) + 2]))
	vcf_aa_mut[[i]] <- vector(mode = "character", length = length(vcf_aa[[i]]))
	vcf_cdna_mut[[i]] <- vector(mode = "character", length = length(vcf_aa[[i]]))
	for (j in 1:length(vcf_aa[[i]])) {
		# take first one only
		vcf_aa_mut[[i]][j] <- vcf_aa[[i]][[j]][ vcf_aa[[i]][[j]] != "" ][1]
		if (is.na(vcf_aa_mut[[i]][j])) {
			vcf_aa_mut[[i]][j] <- ""
		}
		vcf_cdna_mut[[i]][j] <- vcf_cdna[[i]][[j]][ vcf_cdna[[i]][[j]] != "" ][1]
		if (is.na(vcf_cdna_mut[[i]][j])) {
			vcf_cdna_mut[[i]][j] <- ""
		}
	}
	vcf_snpEff[[i]]$AA <- vcf_aa_mut[[i]]
	vcf_snpEff[[i]]$Change <- vcf_cdna_mut[[i]]
}
aa <- vector(mode = "list", length = length(vcf_snpEff))
for (i in 1:length(vcf_snpEff)) {
	aa[[i]] <- data.frame(Chr = "chrM", Pos = vcf_snpEff[[i]]$V2, Ref = vcf_snpEff[[i]]$V4, Alt = vcf_snpEff[[i]]$V5, AA = vcf_snpEff[[i]]$AA, Change = vcf_snpEff[[i]]$Change)
	write.table(aa[[i]], file = bamList$outputAA[i], quote = FALSE, col.names = TRUE, sep = "\t", row.names = FALSE)
}


## HAPLOGREP, FIRST PASS BY VCF
vcfTemp <- sapply(1:length(bamList$outputVCF), function(x) tempfile())
haplogrepRun1_cmd <- paste("haplogrep classify --input=", vcfTemp, " --output=", bamList$haplogrepList, " --hits=", inputFilter$hits, " --format=vcf --extend-report --hetLevel=", inputFilter$AF_max_run1, sep = "") # allele frequency 85%
vcf1 <- suppressWarnings(lapply(bamList$outputVCF, vcf.filter))
for (i in 1:length(vcf1)) {
	cat(vcf1[[i]], sep = "\n", file = vcfTemp[i])
	system(haplogrepRun1_cmd[i])
}

haplogrepRun1 <- unique(rbindlist(lapply(1:length(bamList$haplogrepList), function(x) haplogroup.best(r = bamList$haplogrepList[x], id = bamList$outputID[x]))))

write.table(haplogrepRun1, file = paste("haplogrep_firstPass-", inputSetup$output_timestamp, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
haplogrepRun1_muts <- rbindlist(lapply(1:length(haplogrepRun1$SampleID), function(x) {
	y <- unlist(strsplit(haplogrepRun1$Remaining_Polys[x], "\\) "))
	if (length(y) < 1) {
		y <- NA
	}
	data.frame(ID = haplogrepRun1$SampleID[x], Mut = y)
}))
haplogrepRun1_muts$Pos <- gsub(" \\(.*$", "", haplogrepRun1_muts$Mut)
haplogrepRun1_muts$Mut <- gsub("\\)", "", gsub("^.*\\(", "", haplogrepRun1_muts$Mut))
haplogrepRun1_muts$Key <- paste(haplogrepRun1_muts$ID, haplogrepRun1_muts$Pos, sep = ":")
haplogrepRun1_aac <- rbindlist(lapply(1:length(haplogrepRun1$SampleID), function(x) {
	y <- unlist(strsplit(haplogrepRun1$AAC_In_Remainings[x], "\\] "))
	y <- ifelse(length(y) < 1, yes = NA, no = y)
	data.frame(ID = haplogrepRun1$SampleID[x], Mut = y)
}))
haplogrepRun1_aac$Pos <- gsub(" \\[.*$", "", haplogrepRun1_aac$Mut)
haplogrepRun1_aac$Mut <- gsub("\\]", "", gsub("^.*\\[", "", haplogrepRun1_aac$Mut))
haplogrepRun1_aac$Key <- paste(haplogrepRun1_aac$ID, haplogrepRun1_aac$Pos, sep = ":")


## EAGLE OUTPUT PROCESSING
eagleList <- bamList$outputEagle
eagleList.key <- gsub("\\.vcf$", "", bamList$V3)
names(eagleList.key) <- eagleList
eagleList <- eagleList[match(bamList$outputID, eagleList.key)]
eagleList <- eagleList[!is.na(eagleList)] # in case of no matches
egl <- lapply(eagleList, read.table, sep = "\t", header = FALSE)
egl <- lapply(egl, egl.retain)


## RECALCULATION OF AF BY POSITION (e.g. 1A>C and 1A>T are counted together)
transTables <- lapply(egl, af.calc)
transMergeKey <- unique(unlist(lapply(transTables, function(x) x$Pos)))
transAF <- lapply(transTables, function(x) x$AF[match(transMergeKey, x$Pos)])
transAF <- do.call(cbind, transAF) # replace_na(0) seems broken with lists
transAF <- apply(transAF, 2, na.set, y = 0)
colnames(transAF) <- eagleList
rownames(transAF) <- transMergeKey
write.csv(transAF, file = paste("testTransAF-", inputSetup$output_timestamp, ".csv", sep = ""), row.names = TRUE, quote = FALSE)


## IMPORTING AMINO ACID ANNOTATION
egl.key <- lapply(egl, function(x) paste(x$V2, x$V3, x$V4, sep = ":"))
aa.key <- lapply(aa, function(x) paste(x$Pos, x$Ref, x$Alt, sep = ":"))
for (i in 1:length(egl)) {
	egl[[i]]$Change <- aa[[i]]$Change[match(egl.key[[i]], aa.key[[i]])]
	egl[[i]]$AA <- aa[[i]]$AA[match(egl.key[[i]], aa.key[[i]])]
}


# ANNOTATION OF MITOCHONDRIA GENE SYMBOLS
mitoBed <- unique(data.table::rbindlist(lapply(egl, function(x) data.frame(x$V1, x$V2 -1, x$V2))))
mitoBed <- mitoBed[order(mitoBed[,2])]
mitoBed.f <- tempfile()
write.table(mitoBed, file = mitoBed.f, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
mitoBed.anno <- lapply(system(paste("bedtools intersect -a", mitoBed.f ,"-b", inputSetup$mito_symbol_annotation, "-wao"), intern = TRUE), function(x) unlist(strsplit(x, "\t")))
mitoBed.anno.pos <- sapply(mitoBed.anno, function(x) x[3])
mitoBed.anno.gene <- sapply(mitoBed.anno, function(x) x[7])
mitoBed.anno <- data.table(Pos = mitoBed.anno.pos, Gene = mitoBed.anno.gene)
mitoBed.anno <- mitoBed.anno[,list(Symbol = paste(unique(Gene), collapse = "\t")), c("Pos")]


## GENERATING COMBINED SAMPLE VARIANT LIST
for (i in 1:length(egl)) {
	egl[[i]]$AA[is.na(egl[[i]]$AA)] <- ""
	egl[[i]]$Symbol <- mitoBed.anno$Symbol[match(egl[[i]]$V2, mitoBed.anno$Pos)]
	egl[[i]]$SampleID <- eagleList.key[i]
}
eglS <- rbindlist(egl)
eglS <- eglS[eglS$V2 <= 16569,] # set maximum length for haploGrep
eglS$Symbol <- gsub("\t", "/", eglS$Symbol)
write.table(eglS, file = paste("eagleSummary_all-", inputSetup$output_timestamp, ".txt", sep = ""), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


## CALCULATING Z SCORE CUTOFF FOR STAT-BASED MUTECTION CALLS
AFbyPos <- as.data.frame(transAF)
AFbyPos$Average  <- apply(as.data.frame(transAF), 1, mean, na.rm = TRUE)
AFbyPos$SD <- apply(as.data.frame(transAF), 1, sd, na.rm = TRUE)
AFbyPos$OutlierZSD <- apply(AFbyPos, 1, function(y) {
	m <- sum(y[!(names(y) %in% c("Average", "SD", "OutlierZSD"))] > y["Average"]+ y["SD"] * inputFilter$Z)
	return(m)
})
zAFbyPos <- as.data.frame(transAF)
for (i in 1:dim(zAFbyPos)[1]) {
	zAFbyPos[i,] <- (zAFbyPos[i,] - AFbyPos$Average[i])/AFbyPos$SD[i]
}
colnames(zAFbyPos) <- eagleList.key[colnames(zAFbyPos)]
AFbyPos$Higher3SD <- AFbyPos$Average + (3 * AFbyPos$SD)
AFbyPos2 <- as.data.frame(transAF)
AFbyPos2[(AFbyPos2 < AFbyPos$Higher3SD)] <- NA


## INTEGRATION OF Z SCORE AND HAPLOGREP PASS 1 CALLS
colnames(eglS)[1:10] <- c("Contig", "Pos", "Ref", "Alt", "Reads", "refReads", "altReads", "p", "OR", "nbrs")
eglS$Haplogroup1 <- haplogrepRun1$Haplogroup[match(eglS$SampleID, haplogrepRun1$SampleID)]
eglS$AF <- (eglS$altReads + 0.01)/(eglS$refReads + eglS$altReads + 0.02)
eglS$Z <- NA
eglS$Recurrence <- NA
eglS_match <- data.frame(Pos = match(as.character(eglS$Pos), rownames(zAFbyPos)), ID = match(eglS$SampleID, colnames(zAFbyPos)))
eglS$Z <- as.vector(apply(eglS_match, 1, function(x) zAFbyPos[x[1], x[2]]))
eglS$Recurrence <- AFbyPos$OutlierZSD[eglS_match$Pos]
eglS_MutLength <- nchar(eglS$Alt) - nchar(eglS$Ref)
eglS_MutKey <- ifelse(eglS_MutLength == 0, yes = paste(eglS$Pos, eglS$Alt, sep = ""), no = "")
eglS_MutKey[eglS_MutLength < 0] <- paste(eglS$Pos[eglS_MutLength < 0], "del", sep = "") # del, length not counted
eglS_MutKey[eglS_MutLength > 0] <- 
	paste(eglS$Pos, ".", eglS_MutLength, substr(eglS$Alt, start = nchar(eglS$Ref) + 1, stop = nchar(eglS$Alt)), sep = "")[eglS_MutLength > 0] # insertion
eglS_sample_mutKey <- paste(eglS$SampleID, eglS_MutKey, sep = ":")


## HAPLOGREP RUN 1 RESULT ASSESSMENT
haplogrepRun1_snps <- rbindlist(lapply(1:length(haplogrepRun1$Found_Polys), function(x) {
	y <- unlist(strsplit(haplogrepRun1$Found_Polys[x], " "))
	data.frame(ID = haplogrepRun1$SampleID[x], SNP = y)
}))
haplogrepRun1_snps$Key <- paste(haplogrepRun1_snps$ID, haplogrepRun1_snps$SNP, sep = ":")
eglS$SNP1 <- eglS_sample_mutKey %in% haplogrepRun1_snps$Key
eglS$SNPpass1 <- ifelse(eglS$SNP1, yes = "Haplogrep", no = "")
eglS$SNPpass1[which(eglS$SNPpass1 == "" & (eglS$AF >= inputFilter$AF_max_run1))] <- "Probable SNP (> AF cutoff)"

eglS$MutPass1 <- haplogrepRun1_muts$Mut[match(eglS_sample_mutKey, haplogrepRun1_muts$Key)]
eglS$MutPass1[is.na(eglS$MutPass1)] <- ""

eglS$AACPass1 <- haplogrepRun1_aac$Mut[match(eglS_sample_mutKey, haplogrepRun1_aac$Key)]
eglS$AACPass1[is.na(eglS$AACPass1)] <- ""


## HAPLOGREP RUN 2, SNP: AF >= 80%, MUT: 20% < AF <80%, Z > 3
sampleHSD <- data.table(ID = eglS$SampleID, Mut = eglS_MutKey, SNP = eglS$SNP1, AF = eglS$AF, Z = eglS$Z, Range = "1-16569", Haplogroup = "?")
sampleHSD$Found_Polys <- haplogrepRun1$Found_Polys[match(sampleHSD$ID, haplogrepRun1$SampleID)]
sampleHSD_filter <- sort(unique(c(
	which(sampleHSD$SNP), # SNP from Run 1
	which((!sampleHSD$SNP) & (is.between(sampleHSD$AF, inputFilter$AF_min_run2, inputFilter$AF_max_run2)) & (sampleHSD$Z >= inputFilter$Z)), # Mutations with Z >= cutoff and AF between specified range
	which((!sampleHSD$SNP) & (sampleHSD$AF > inputFilter$AF_min_run2))))) # SNPs
sampleHSD <- sampleHSD[sampleHSD_filter,]
sampleHSD <- sampleHSD[,list(
	Polymorphisms = paste(sort(unique(
		c(Mut, unlist(strsplit(Found_Polys, " "))))), collapse = "\t")), c("ID", "Range", "Haplogroup")]
sampleHSD_nochange <- which(is.na(match(haplogrepRun1$SampleID, sampleHSD$ID)))
sampleHSD_update <- haplogrepRun1[sampleHSD_nochange, c("SampleID", "Range", "Haplogroup", "Found_Polys")]
colnames(sampleHSD_update) <- c("ID", "Range", "Haplogroup", "Polymorphisms")
sampleHSD_update$Polymorphisms <- gsub(" ", "\t", sampleHSD_update$Polymorphisms)
sampleHSD <- rbind(sampleHSD, sampleHSD_update)

## Manual fixes
#  implemented as command-line version of Haplogrep treats HSD and VCF inputs differently, where VCF will ignore entries with the following variants
#  while when HSD files are used as inputs the program will refuse to proceed. Manual inspections at problematic loci in IGV are done to determine
#  whether these putative mutections detected by Mutect are in fact present or otherwise artifacts
#sampleHSD$Polymorphisms <- gsub("10813CA", "10814A", sampleHSD$Polymorphisms) # acceptable as is by haplogrep
sampleHSD$Polymorphisms <- gsub("14346TG", "14346T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("704CC", "704C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("5601TT", "5601T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16092CC", "16092C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16092CT", "", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("10397GG", "10397G\t10398G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("8251AC", "8251A", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("4769GA", "4769G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("309TC", "310C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("12704AT", "12704A\t12705T", sampleHSD$Polymorphisms) 
sampleHSD$Polymorphisms <- gsub("4769GG", "4769G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16297CC", "16297C\t16298C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("8860GC", "8860G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("4833GC", "4833G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16399GT", "16399G\t16400T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16311CA", "16311C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("5601TC", "5601T", sampleHSD$Polymorphisms)
#  09/12/22
sampleHSD$Polymorphisms <- gsub("13753CCC", "13753C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("2768GAC", "2768G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16642GC", "16642G\t16643C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("3118.1C", "", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("6319.1G", "", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("9902.1A", "", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("13676GC", "13676G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("14668TA", "14668T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16126CA", "16126C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16518AC", "16518A\t16519C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("3394CA", "3394C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("4824GC", "4824G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("68AA", "68A\t69A", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("709AC", "709A\t710C", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("8412CAC", "8412C", sampleHSD$Polymorphisms)
# 10/13/22
sampleHSD$Polymorphisms <- gsub("16126CG", "16126C\t16127G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("3394CG", "3394C\t3395G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("68AG", "68A", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("709AT", "709A", sampleHSD$Polymorphisms)
# 11/07/22
sampleHSD$Polymorphisms <- gsub("10152AA", "10152A", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("10398GCC", "10398G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("16223TT", "16223T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("2706GA", "2706G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("3380AA", "3380A", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("3625AC", "3625A", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("3992TA", "3992T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("750GA", "750G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("751TT", "751T\t752T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("9494GT", "9494G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("2706GG", "2706G\t2707G", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("2706GT", "2706G\t2707T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("5809AA", "5809A", sampleHSD$Polymorphisms)
# 11/08/22
sampleHSD$Polymorphisms[match("L-C52", sampleHSD$ID)] <- gsub("8248.1G", "", sampleHSD$Polymorphisms[match("L-C52", sampleHSD$ID)])
sampleHSD$Polymorphisms[match("L-K27", sampleHSD$ID)] <- gsub("3912.1C", "", sampleHSD$Polymorphisms[match("L-K27", sampleHSD$ID)])
sampleHSD$Polymorphisms[match("L-K37", sampleHSD$ID)] <- gsub("3912.1C", "", sampleHSD$Polymorphisms[match("L-K37", sampleHSD$ID)])
sampleHSD$Polymorphisms[match("L-K9NT", sampleHSD$ID)] <- gsub("3912.1C", "", sampleHSD$Polymorphisms[match("L-K9NT", sampleHSD$ID)])
sampleHSD$Polymorphisms <- gsub("9308GT", "9308G", sampleHSD$Polymorphisms)
# 12/27/22
sampleHSD$Polymorphisms <- gsub("14667GT", "14667G\t14668T", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("6293CT", "", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("8989AC", "8989A", sampleHSD$Polymorphisms)
sampleHSD$Polymorphisms <- gsub("\t+", "\t", sampleHSD$Polymorphisms)

## Haplogrep, second pass
haplogrepRun2_input <- paste(sampleHSD$ID, "_haplogrep2_input.hsd", sep = "")
haplogrepRun2_output <- paste(sampleHSD$ID, "_haplogrep2_output.txt", sep = "")
haplogrepRun2_cmd  <- paste("haplogrep classify --input=", haplogrepRun2_input, " --output=", haplogrepRun2_output, " --hits=", inputFilter$hits, " --format=hsd --extend-report 2>&1", sep = "")
hr2_exception <- c()
for (i in 1:length(sampleHSD$ID)) {
	write.table(sampleHSD[i,], file = haplogrepRun2_input[i], col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
	hr2_res <- system(haplogrepRun2_cmd[i], intern = TRUE)
	hr2_cur_err <- grep("exceptions.parse.sample.HsdFileSampleParseException", hr2_res, value = TRUE)
	if (length(hr2_cur_err) > 0) {
		hr2_cur_err <- gsub("^exceptions\\.parse\\.sample\\.HsdFileSampleParseException.* polymorphism ", "", hr2_cur_err)
		hr2_cur_err <- gsub(" contains some.*$", "", hr2_cur_err)
		hr2_exception <- c(hr2_exception, paste("sampleHSD$Polymorphisms <- gsub(\"", hr2_cur_err, "\", \"", hr2_cur_err, "\", sampleHSD$Polymorphisms) # ", haplogrepRun2_input[i], sep = ""))
	}
	remove(hr2_res, hr2_cur_err)
}
cat(sort(hr2_exception), sep = "\n") # displays errors that require manual inspections during the run. repeat this section until errors disappear

## INTEGRATION OF HAPLOGREP PASS 2 CALLS
haplogrepRun2 <- unique(rbindlist(lapply(1:length(haplogrepRun2_output), function(x) haplogroup.best(r = haplogrepRun2_output[x], id = sampleHSD$ID[x]))))
haplogrepRun2_muts <- rbindlist(lapply(1:length(haplogrepRun2$SampleID), function(x) {
	y <- unlist(strsplit(haplogrepRun2$Remaining_Polys[x], "\\) "))
	if (length(y) < 1) {
		y <- NA
	}
	data.frame(ID = haplogrepRun2$SampleID[x], Mut = y)
}))
haplogrepRun2_muts$Pos <- gsub(" \\(.*$", "", haplogrepRun2_muts$Mut)
haplogrepRun2_muts$Mut <- gsub("\\)", "", gsub("^.*\\(", "", haplogrepRun2_muts$Mut))
haplogrepRun2_muts$Key <- paste(haplogrepRun2_muts$ID, haplogrepRun2_muts$Pos, sep = ":")
haplogrepRun2_aac <- rbindlist(lapply(1:length(haplogrepRun2$SampleID), function(x) {
	y <- unlist(strsplit(haplogrepRun2$AAC_In_Remainings[x], "\\] "))
	if (length(y) < 1) {
		y <- NA
	}
	data.frame(ID = haplogrepRun2$SampleID[x], Mut = y)
}))
haplogrepRun2_aac$Pos <- gsub(" \\[.*$", "", haplogrepRun2_aac$Mut)
haplogrepRun2_aac$Mut <- gsub("\\]", "", gsub("^.*\\[", "", haplogrepRun2_aac$Mut))
haplogrepRun2_aac$Key <- paste(haplogrepRun2_aac$ID, haplogrepRun2_aac$Pos, sep = ":")
eglS$Haplogroup2 <- haplogrepRun2$Haplogroup[match(eglS$SampleID, haplogrepRun2$SampleID)]
haplogrepRun2_snps <- rbindlist(lapply(1:length(haplogrepRun2$Found_Polys), function(x) {
	y <- unlist(strsplit(haplogrepRun2$Found_Polys[x], " "))
	data.frame(ID = haplogrepRun2$SampleID[x], SNP = y)
}))
haplogrepRun2_snps$Key <- paste(haplogrepRun2_snps$ID, haplogrepRun2_snps$SNP, sep = ":")
eglS$SNP2 <- eglS_sample_mutKey %in% haplogrepRun2_snps$Key
eglS$SNPpass2 <- eglS$SNP2 & (eglS$AF >= inputFilter$AF_min_run2) 
eglS$SNPpass2 <- ifelse(eglS$SNP2, yes = "Haplogrep", no = "")
eglS$SNPpass2[which(eglS$SNPpass2 == "" & (eglS$AF >= inputFilter$AF_max_run2))] <- "Probable SNP (> AF cutoff)"
eglS$MutPass2 <- haplogrepRun2_muts$Mut[match(eglS_sample_mutKey, haplogrepRun2_muts$Key)]
eglS$MutPass2[is.na(eglS$MutPass2)] <- ""
eglS$AACPass2 <- haplogrepRun2_aac$Mut[match(eglS_sample_mutKey, haplogrepRun2_aac$Key)]
eglS$AACPass2[is.na(eglS$AACPass2)] <- ""

# CONSENSUS CALLS AND SUMMARIZING OUTPUT TABLES
eglS1 <- eglS[,list(Classification1 = paste(SNPpass2, MutPass2, sep = "/"), Classification2 = paste(SNPpass2, MutPass2, sep = "/"))]
eglS$Classification1 <- eglS1$Classification1
eglS$Classification1 <- gsub("\\/$", "", gsub("^\\/", "", eglS$Classification1))
eglS$Classification1 <- gsub("hotspot", "hospot mutation", eglS$Classification1)
eglS$Classification2 <- eglS1$Classification2
eglS$Classification2 <- gsub("\\/$", "", gsub("^\\/", "", eglS$Classification2))
eglS$Classification2 <- gsub("hotspot", "hospot mutation", eglS$Classification2)

# Merging 2 runs
eglS$cHaplogroup <- ifelse(eglS$Haplogroup1 == eglS$Haplogroup2, yes = eglS$Haplogroup1, no = paste(eglS$Haplogroup1, eglS$Haplogroup2, sep = "/"))
eglS$cSNP <- NA
eglS$cSNP[sort(c(which(eglS$SNPpass1 != ""), which(eglS$SNPpass2 != "")))] <- TRUE
eglS$cMut <- NA
eglS$cMut[sort(c(which(eglS$MutPass1 != ""), which(eglS$MutPass2 != "")))] <- TRUE
eglS$cClass <- ifelse(eglS$Classification1 == eglS$Classification2, yes = eglS$Classification1, no = paste(eglS$Classification1, eglS$Classification2, sep = "/"))
eglS$AA[eglS$Symbol == "."] <- "" # remove AA cals from intergenic regions
write.table(eglS, file = paste("haplogrep_summary_2runs-", inputSetup$output_timestamp, ".txt", sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")