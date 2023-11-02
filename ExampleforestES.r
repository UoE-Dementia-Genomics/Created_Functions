# Example forestES plot

# Create list object of effect sizes and standard errors
# It will plot first column at bottom and last column at top
# rownames should be what you call as x
ES <- data[,c("CB.Effect_Fixed", "AZ1.CB.Estimate", "LBB2.CB.Estimate", 
"LBB1.CB.Estimate", "meta.es", "EC.Effect_Fixed", "LBB2.E.Estimate", 
"LBB1.E.Estimate", "TG.Effect_Fixed", "AZ2.F.Estimate", "AZ1.F.Estimate", 
"MS.F.Estimate", "LBB1.F.Estimate", "PFC.Effect_Fixed", "ROSMAP.A.Estimate", 
"MS.A.Estimate", "LBB1.A.Estimate")]

SE <- data[,c("CB.SE_Fixed", "AZ1.CB.Std.Error", "LBB2.CB.Std.Error", 
"LBB1.CB.Std.Error", "meta.se", "EC.SE_Fixed", "LBB2.E.SE", 
"LBB1.E.SE", "TG.SE_Fixed", "AZ2.F.SE", "AZ1.F.SE", "MS.F.SE", "LBB1.F.SE", 
"PFC.SE_Fixed", "ROSMAP.A.SE", "MS.A.SE", "LBB1.A.SE")]

dat <- list(ES = ES, SE = SE)

# Which levels are polygons?
poly <- c(1, 5, 6, 9, 14)

# Colours of levels
col <- c(replicate(4, "#ff8e04"), "black", replicate(3, "blue"), 
replicate(5, "forestgreen"), replicate(4, "red"))

# Colours of polygons
polycol <- c("#ff8e04", "black", "blue", "forestgreen", "red")

# Where do you want separation lines
sepline <- c(4.5, 5.5, 8.5, 13.5)

# Where and what do you want section headers to be
septitle <- c("Cerebellum", "Entorhinal \nCortex", 
"Superior/Middle \nTemporal Gyrus", "Prefrontal \nCortex")
septitlepos <- c(2.5, 7, 11, 15.5)

# Labels for each level
names <- c("Pooled", "Arizona 1", "London 2", "London 1", "Cross-Cortex", 
"Pooled", "London 2", "London 1", "Pooled", "Arizona 2", "Arizona 1", 
"Mount Sinai", "London 1", "Pooled", "ROSMAP", "Mount Sinai", "London 1")

# Plot, multiply is multiplication factor for ES and SE
forestES("cg22962123", data = dat, multiply = 6, col = col, polycol = polycol, 
names = names, septitle = septitle, septitlepos = septitlepos, 
xlab = "Effect Size (beta)", sepline = sepline, poly = poly)
