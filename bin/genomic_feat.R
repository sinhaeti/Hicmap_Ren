## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Mmusculus.UCSC.mm9")
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm9)
require(parallel) 

segment_chromsome <- function(l_seg, genome, sel.chr){
	seqs <- seq(0, round(seqlengths(genome)[[sel.chr]]), l_seg)
	segments <- data.frame(chr=sel.chr, start=seqs[1:length(seqs)-1], end=seqs[-1])
	return(IRanges(start=seqs[1:length(seqs)-1], end=seqs[-1]))
}

extract_feat_chromsome <- function(genome, sel.chr, bin_size, cs.pos, cs.neg){
	seg.chr.ir <- segment_chromsome(bin_size, genome, sel.chr)
	seg.chr.gr <- GRanges(sel.chr, seg.chr.ir, elen=rep(0, length(seg.chr.ir)), nGC=rep(0, length(seg.chr.ir)))	
	cs.pos.chr.ir <- reduce(IRanges(start=cs.pos$start[cs.pos$chr==sel.chr], end=cs.pos$end[cs.pos$chr==sel.chr]))
	cs.neg.chr.ir <- reduce(IRanges(start=cs.neg$start[cs.neg$chr==sel.chr], end=cs.neg$end[cs.neg$chr==sel.chr]))
	
	# Step 2. chose all cutter inverals
	ov.chr.pos <- findOverlaps(seg.chr.ir, cs.pos.chr.ir)
	ov.chr.neg <- findOverlaps(seg.chr.ir, cs.neg.chr.ir)
	
	# step 3. find all overlaps	
	ids.pos <- ov.chr.pos@queryHits
	ids.neg <- ov.chr.neg@queryHits

	wids.pos <- ranges(ov.chr.pos, seg.chr.ir, cs.pos.chr.ir)@width
	wids.neg <- ranges(ov.chr.neg, seg.chr.ir, cs.neg.chr.ir)@width
	tmp = rbind(data.frame(id=ids.pos, wid=wids.pos), data.frame(id=ids.neg, wid=wids.neg))
	
	tmp = tmp[order(tmp[,1]),]
	tmp.split <- split(tmp, tmp$id)
	res.ls <- lapply(tmp.split, function(x){ data.frame(seg_id=x$id[1], wid_sum=sum(x$wid))})
	res.df <- do.call(rbind, res.ls)
	seg.chr.gr$elen[res.df$seg_id] = res.df$wid_sum/(bin_size*2)

	# count GC content
	genome.chr <- genome[[sel.chr]]
	ref.C <- matchPattern("C", genome.chr, max.mismatch=0)
	ref.G <- matchPattern("G", genome.chr, max.mismatch=0)		
	all.genomic <- getSeq(genome, sel.chr, start(seg.chr.gr)+1, end(seg.chr.gr))									  
	nC = do.call(rbind, lapply(vmatchPattern("C", all.genomic, max.mismatch=0), length))
	nG = do.call(rbind, lapply(vmatchPattern("G", all.genomic, max.mismatch=0), length))
	seg.chr.gr$nGC = (nC + nG)/bin_size	
	return(seg.chr.gr)	
}

extract_feat_genome <- function(genome, chroms, num.cores, bin_size, cs.pos, cs.neg){ 
	res <- mclapply(chroms, function(x){ extract_feat_chromsome(genome, x, bin_size, cs.pos, cs.neg) }, mc.cores=num.cores)
	return(do.call(c, res))
}

bin_size = 5000
num.cores = 10
genome <- BSgenome.Mmusculus.UCSC.mm9

# step 2. read in effective cutter interval
cs.pos <- read.table("data/mm9.MboI.500bp.pos.merged.bed")
cs.neg <- read.table("data/mm9.MboI.500bp.neg.merged.bed")
colnames(cs.pos) = colnames(cs.neg) = c("chr", "start", "end")

chroms <- paste("chr", c(1:19, "X", "Y"), sep="")
res <- extract_feat_genome(genome, chroms, num.cores, bin_size, cs.pos, cs.neg)

write.table(as.data.frame(res)[,c(1,2,3,6,7)], file = "genomicFeatures.mm9.MboI.5k.bed", append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = FALSE, qmethod = c("escape", "double"),
                 fileEncoding = "")
				 