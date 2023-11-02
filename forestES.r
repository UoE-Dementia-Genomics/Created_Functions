forestES<- function (x, data = NULL, xlim = NULL, ylim = NULL, xlab = NULL, 
ylab = NULL, pch = NULL, main = NULL, col = NULL, multiply = NULL, cex = NULL, 
names = NULL, poly = NULL, polycol = NULL, sepline = NULL, seplinecol = NULL,
septitle = NULL, septitlepos = NULL, septitleposx = NULL, namesx = NULL){

# 'data' must be list of 1) dataframe of effect sizes 2) dataframe of 
# standard errors with columns as analyses/cohorts and rows of cpgs or locations. 
# Columns should be in same order. Will put first column on bottom and last 
# column on top
# x is cpg site (must be rownames)
# multiply is multiplication factor added to ES and SE
# names are what lines of forest should be labelled
# poly are lines you want to be represented by diamonds rather than points 
# and lines
# polycol what colours do you want the poly to be (from bottom to top)
# sepline is where you want to put separation lines (0.5 between points)
# septitle is section labels for separation
# septitlepos where do you want section labels
# namesx position of names on x axis (is autimatically calculated but may 
# need adjustment)

Est<-as.numeric(data[[1]][x,])

SE<-as.numeric(data[[2]][x,])

if(is.null(ylim)) {
	ylim=c(0.6,length(Est))
	}
	
if(is.null(xlab)) {
	xlab="Effect Size"
	}

if(is.null(ylab)) {
	ylab=""
	}	
	
if(is.null(pch)) {
	pch = 19 
	}

if(is.null(main)) {
	main = x
	}

if(is.null(col)) {
	col = "black"
	}
	
if(is.null(multiply)) {
	multiply = 1
	}

if(is.null(seplinecol)) {
	seplinecol = "gray85"
	}

if(is.null(names)) {
	names=colnames(data[[1]])
	}

if(is.null(poly)) {
	poly=NULL
	}
	
if(is.null(polycol)) {
	polycol=replicate(length(poly),"black")
	}
		
Est*multiply->Est2
SE*multiply->SE2

if(is.null(xlim)) {
      xlim <- c(min(Est2-SE2)-0.05,max(Est2+SE2)+0.02)
    }

if(is.null(septitleposx)) {
	septitleposx = min(xlim)
	}
	
if(is.null(namesx)) {
	namesx<-min(xlim) + 0.01
	}

par(xpd=TRUE)

plot(1:length(Est2)~Est2, xlim=xlim, axes=F, ylim=ylim, pch=pch, main=main, 
col=col, cex=cex, ylab=ylab, xlab=xlab)
segments((Est2 - SE2), 1:length(Est2), (Est2 + SE2), 1:length(Est2), col=col)
abline(v=0, xpd=FALSE)
axis(1)

for(i in 1:length(poly)){
polygon(c(Est2[poly[i]] - SE2[poly[i]], Est2[poly[i]],Est2[poly[i]] + 
SE2[poly[i]],Est2[poly[i]]), c(poly[i], poly[i]-0.2, poly[i], poly[i]+0.2), 
col=polycol[i], border=polycol[i])
}

for (i in 1:length(names)){
text(namesx, i, names[i], pos = 4)
}

for (i in 1:length(sepline)){
abline(h=sepline[i], xpd=FALSE, col=seplinecol)
}

for (i in 1:length(septitle)){
text(septitleposx, septitlepos[i], septitle[i], srt=90, cex=1)
}

}
