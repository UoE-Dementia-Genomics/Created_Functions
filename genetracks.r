genetracks <- function (data = NULL, chr = NULL, cpg = NULL, range = NULL, 
xlim = NULL, ylim = NULL, col = NULL, cexgene = NULL, result = NULL, pad = NULL, 
genelines = NULL, cpgcol = NULL){

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(Homo.sapiens)
require(AnnotationHub)

if(is.null(pad)) {
	pad = 30000
	}
		
if(missing(chr)) {
                        
if(missing(cpg))        {

stop("please provide either 'chr' and 'range' or 'cpg'")
                
                }	
		                }

if(! is.null(cpg)){

	data[cpg,"CHR"]->chr
	data[cpg,"MAPINFO"]->MI
	range = MI:MI
	
    }	
	
if(is.null(xlim)) {
      xlim <- c((min(range) - pad),(max(range) + pad))
    }

par(xpd=TRUE)

data[which(data$CHR == chr & data$MAPINFO >= min(range) - 
pad & data$MAPINFO <= max(range) + pad),]->y

if(is.null(chr)) {
	chr = y$CHR
	}
	
if(is.null(cexgene)) {
	cex = 0.5
	}
      
if(is.null(cpgcol)){
    cpgcol="forestgreen"
    } 
	
chr2=paste("chr", chr, sep="")

hub <- AnnotationHub()
query(hub, c("cpg","hg19"))
cpgs <- as.data.frame(hub[["AH5086"]])
cpgs[which(cpgs$seqnames == chr2 & cpgs$start >=min(range)-pad & 
cpgs$start <= max(range)+pad),]->cpg
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

gr <- GRanges(seqnames = chr2, ranges = IRanges(start = min(range)-pad, 
end = max(range)+pad))

subsetByOverlaps(transcripts(txdb), gr)->trans

res <- as.data.frame(transcripts(Homo.sapiens, columns=c("TXNAME","SYMBOL")))
as.data.frame(res[res$TXNAME %in% as.data.frame(trans)[,"tx_name"],"SYMBOL"])->m
colnames(m)<-"SYMBOL"
cbind(trans,m)->trans
trans[is.na(trans$SYMBOL) == "FALSE",]->trans

if(is.null(genelines)){
	genelines = nrow(trans)
	}
	
rep(1:max(genelines), times=ceiling(nrow(trans)/max(genelines)))->gl
gl[1:nrow(trans)]->gl
trans$gl <- gl
	
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)->exons
exons[rownames(trans)]->exons

par(mar=c(1,3,1,2)+ 0.1) 
plot(0,0,type="n", xlim=xlim, ylim=c(0,genelines +1), axes=FALSE, xlab="", 
ylab="", xpd=FALSE)

for(g in 1:nrow(cpg)){
polygon(c(cpg[g, "start"], cpg[g, "end"], cpg[g, "end"], cpg[g, "start"]), 
c((genelines+1)-0.2, (genelines+1)-0.2, (genelines+1)+0.2,(genelines+1)+0.2),
col=cpgcol, xpd=FALSE)
                     }

for(i in 1:nrow(trans)){

lines(c(as.data.frame(trans)[i,"start"], as.data.frame(trans)[i,"end"]), 
c(trans$gl[i],trans$gl[i]), xpd=FALSE)
text(as.data.frame(trans)[i,"start"]+((as.data.frame(trans)[i,"end"] - 
as.data.frame(trans)[i,"start"])/2) ,trans$gl[i] + 0.4 , trans$SYMBOL[i], 
cex=cexgene, xpd=FALSE)

for(x in 1:nrow(as.data.frame(exons[[i]]))){

polygon(c(as.data.frame(exons[[i]])[x,"start"], 
as.data.frame(exons[[i]])[x,"end"], as.data.frame(exons[[i]])[x,"end"], 
as.data.frame(exons[[i]])[x,"start"]), 
c(trans$gl[i]-0.2, trans$gl[i]-0.2, trans$gl[i]+0.2,trans$gl[i]+0.2),
col="black", xpd=FALSE)


                                           }
                      }
                      }

