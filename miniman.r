miniman <- function (data = NULL, chr = NULL, cpg = NULL, range = NULL, 
xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, pch = 1, main = NULL, 
col = "black", cex = 1, result = NULL, pad = 30000, multiply = NULL, 
nullcol = "black", negcol = "black", poscol = "black", ESlevel = 0, ESdat = NULL, 
cexgene = 0.5, geneplot = TRUE, genelines = NULL, chrcol = NULL, 
mapcol = NULL, cpgcol = "forestgreen", cex.axis = 1, cex.lab = 1, sigline = NULL, 
siglevel = 1e-7, siglinecol = NULL) {

# data should just be a dataframe
# chrcol and mapcol will automatically use 'CHR' and 'MAPINFO' but can be 
# replaced. chrcol column should just be numeric
# Can either provide a chr and range or cpg (must be rowname and have 
# mapinfo and chr in the file)
# result is the column who's p-values you want to plot
# pad is amount of padding to add ot the coordinates given
# cexgene is text size for labelling genes
# geneplot will be added by default but can be stopped if = "FALSE"
# genelines will autimatically be calculated by number of transcripts but 
# can be adjusted. If more transtcripts than genelines, several will be plotted 
# on the same line
# ESdat is the effect size column to colour points by
# ESlevel is target effect size to colour by
# nullcol is colour to plot if data doesn't meet ESlevel
# negcol is colour to plot if data meets -ESlevel 
# poscol is colour to plot if data meets +ESlevel 
# cpgcol is colour for cpg island track

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(Homo.sapiens)
require(AnnotationHub)

if(missing(data))        {
stop("please provide data")
                }	

if(is.null(xlab)) {
	xlab = "Genomic Position"
	}

if(is.null(multiply)){
	multiply=1
	}
	
if(is.null(ylab)) {
	ylab = "-log10(p)"
	}	

if(is.null(siglinecol)) {
siglinecol = "red"
	}

if(is.null(main)) {
	main = ""
	}

if(is.null(chrcol)){
    chrcol="CHR"
    }
    
if(is.null(mapcol)){
    mapcol="MAPINFO"
    }    
      
if(missing(chr)) {
                        
if(missing(cpg)){

stop("please provide either 'chr' and 'range' or 'cpg'")
                
                }	
		          }

if(! is.null(cpg)){

	data[cpg,chrcol]->chr
	data[cpg,mapcol]->MI
	range = MI:MI
	
                  }	
	
if(is.null(xlim)) {
      xlim <- c((min(range) - pad),(max(range) + pad))
                  }

if(missing(result)){

stop("please provide result you would like to plot")
                
                   }	

par(xpd=TRUE)

data[which(data[,chrcol] == chr & data[,mapcol] >= min(range) - 
pad & data[,mapcol] <= max(range) + pad),]->y

if(is.null(ylim)) {
	ylim = c(min(-log10(y[,result])),max(-log10(y[,result])))
	}
	
if(is.null(chr)) {
	chr = y[,chrcol]
	}
	
if(is.null(cexgene)) {
	cexgene = 0.5
	}
	
if(is.null(geneplot)) {
	geneplot = TRUE
	}

### Dataplot

ESlevel <- ESlevel/multiply

if(geneplot == "TRUE"){
layout(matrix(c(1,2), 2, 1, byrow = TRUE))
par(mar=c(3,3,1,2)+ 0.1, mgp=c(2,1,0)) 

if(!is.null(ESdat)){

plot(-log10(y[which(y[,ESdat] > -ESlevel & y[,ESdat] < ESlevel),result])~
y[which(y[,ESdat] > -ESlevel & y[,ESdat] < ESlevel),mapcol], xlim = xlim, 
ylim = ylim, cex = cex, pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, 
col=nullcol, cex.axis = cex.axis, cex.lab = cex.lab)
points(-log10(y[which(y[,ESdat] >= ESlevel),result])~
y[which(y[,ESdat] >= ESlevel),mapcol], col=poscol, cex = cex)
points(-log10(y[which(y[,ESdat] <= -ESlevel),result])~
y[which(y[,ESdat] <= -ESlevel),mapcol], col=negcol, cex = cex)
if(sigline == TRUE){
	siglevellog = -log10(siglevel)
	abline(h = siglevellog, col=siglinecol, xpd=FALSE)
} 
                   }

else{

plot(-log10(y[,result])~y[,mapcol], xlim = xlim, ylim = ylim, cex = cex, 
pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, col=nullcol)
if(sigline == TRUE){
	siglevellog = -log10(siglevel)
	abline(h = siglevellog, col=siglinecol, xpd=FALSE)
} 
     }

### Geneplots

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

subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)->exons
trans[rownames(trans) %in% names(exons),]->trans

if(is.null(genelines)){
	genelines = nrow(trans)
	                  }

rep(1:max(genelines), times=ceiling(nrow(trans)/max(genelines)))->gl
gl[1:nrow(trans)]->gl
trans$gl <- gl

exons[rownames(trans)]->exons

par(mar=c(1,3,1,2)+ 0.1) 
plot(0,0,type="n", xlim=xlim, ylim=c(0,genelines +1), axes=FALSE, xlab="", 
ylab="", xpd=FALSE)

if ((nrow(cpg) >= 1)){
for(g in 1:nrow(cpg)){
polygon(c(cpg[g, "start"], cpg[g, "end"], cpg[g, "end"], cpg[g, "start"]), 
c((genelines+1)-0.2, (genelines+1)-0.2, (genelines+1)+0.2,(genelines+1)+0.2),
col=cpgcol, xpd=FALSE, border=cpgcol)
                     }
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
 
if(geneplot == "FALSE"){
plot.new()
if(!is.null(ESdat)){
par(mar=c(3,3,1,2)+ 0.1, mgp=c(2,1,0)) 

plot(-log10(y[which(y[,ESdat] > -ESlevel & y[,ESdat] < ESlevel),result])~
y[which(y[,ESdat] > -ESlevel & y[,ESdat] < ESlevel),mapcol], xlim = xlim, 
ylim = ylim, cex = cex, pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, 
col=nullcol, cex.axis = cex.axis, cex.lab = cex.lab)
points(-log10(y[which(y[,ESdat] >= ESlevel),result])~
y[which(y[,ESdat] >= ESlevel),mapcol], col=poscol, cex = cex)
points(-log10(y[which(y[,ESdat] <= -ESlevel),result])~
y[which(y[,ESdat] <= -ESlevel),mapcol], col=negcol, cex = cex)

                   }

else{

plot(-log10(y[,result])~y[,mapcol], xlim = xlim, ylim = ylim, cex = cex, 
pch = pch, xlab = xlab, ylab = ylab , xpd=FALSE, col=nullcol, 
cex.axis = cex.axis, cex.lab = cex.lab)

     }

                        }
                                                                }
