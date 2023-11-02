# Example miniman

# pad is padding around site or range. Default is 30,000
# result is column with p values to plot
# cexgene is cex for text on gene names
# geneplot is 'TRUE' by default. Logical, should geneplot be drawn?
# genelines is number of rows for genes
# nullcol, negcol and poscol relate to effect size colours, no, negative and 
# positive 
# ESdat is column containing effect size
# ESlevel is effect size level to colour by
# multiply is multiplication factor for effect size

# Plot by cpg or row name. Object must contain chromosome and location 
miniman(data = data, cpg = "cg12307200", result = "meta.p", pad = 300000, 
ESlevel = 0.01, ESdat = "meta.es", negcol = "red", poscol = "forestgreen", 
cexgene = 0.6, multiply = 6, cpgcol ="forestgreen", sigline = TRUE, 
siglinecol = "red", siglevel = 1e-8)

# Plot by chr and genomic position. Object must contain chromosome and 
# location columns
miniman(data = data, chr = 7, range = 27153212:27154305, result = "meta.p", 
pad = 300000, ESlevel = 0.01, ESdat = "meta.es", negcol = "red", 
poscol = "forestgreen", cexgene = 0.6, multiply = 6, cpgcol ="forestgreen",
sigline = TRUE, siglinecol = "red", siglevel = 1e-8)

# Plot by chr and genomic position. Object must contain chromosome and 
# location columns. Includes genelines option
miniman(data = data, chr = 19, range = 45407868:45412647, result = "meta.p", 
pad = 300000, ESlevel = 0.01, ESdat = "meta.es", negcol = "red", 
poscol = "forestgreen", cexgene = 0.6, multiply = 6, cpgcol ="forestgreen"
genelines = 5, sigline = TRUE, siglinecol = "red", siglevel = 1e-8)
