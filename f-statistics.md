# f-statistics

## Getting Started

Reading materials: [Patterson et al 2012](https://www.genetics.org/content/192/3/1065) ; [Peter 2016](https://www.genetics.org/content/202/4/1485) ; [Harris and DeGiorgio](https://doi.org/10.13110/humanbiology.89.1.02)

## Prerequisites

```bash
#dataset
data=/mnt/NAS/workshop/dataset/ho

#software
admixtools=/usr/local/sw/AdmixTools/bin

#fstatfolder
fstatfolder=/mnt/NAS/workshop/f-stat
```

## 1. Outgroup f3

The outgroup f3 statistic measures **shared genetic drift** between two populations, given an **outgroup**, using allele frequency data across the genome.

[AdmixTools](https://github.com/DReichLab/AdmixTools) calculates the confidence in the estimates using a block jackknife approach.

### 1.1. Modern vs Ancient

We will calculate f3 between:
- an Anatolian Neolithic population, either Boncuklu or Tepecik
- a West Eurasian present-day population,
using Yoruba from west Africa as outgroup.

#### 1.1.1. Run f3

```bash
mkdir myfstats
cd myfstats
cat > parf3.modern <<EOF
genotypename:    ${data}.geno
snpname:         ${data}.snp
indivname:       ${data}.ind
popfilename:     ${fstatfolder}/f3out.modern.poplist
EOF
```

The ``f3out.modern.poplist`` contains 3 populations' names on each line:
- an Anatolian Neolithic population, either Boncuklu or Tepecik
- a West Eurasian present-day population
- Yoruba, as outgroup


```bash
head $fstatfolder/f3out.modern.poplist
```

```
Boncuklu Abkhasian Yoruba
Boncuklu Adygei Yoruba
Boncuklu Albanian Yoruba
Boncuklu Armenian Yoruba
Boncuklu Balkar Yoruba
Boncuklu Basque Yoruba
Boncuklu BedouinA Yoruba
Boncuklu BedouinB Yoruba
Boncuklu Belarusian Yoruba
Boncuklu Bulgarian Yoruba
```

Now run f3 on the data.

```bash
nohup $admixtools/qp3Pop -p parf3.modern > f3out.modern.log &
# because the .log file contains extra lines, we grep for lines starting with the standard "result" string
grep result f3out.modern.log > f3modern.out
```


#### 1.1.2. Visualization in R: Raw estimates

```R
tmp<-read.table("f3modern.out")
head(tmp)
f3 <- tmp[,c(2,3,5,6)]
colnames(f3) <- c("B", "A", "f3", "se")
dim(f3)
head(f3)
```

```
         B         A       f3       se
1 Boncuklu Abkhasian 0.159139 0.001397
2 Boncuklu    Adygei 0.157850 0.001360
3 Boncuklu  Albanian 0.165172 0.001455
4 Boncuklu  Armenian 0.159583 0.001376
5 Boncuklu    Balkar 0.157370 0.001398
6 Boncuklu    Basque 0.167222 0.001410
```

Now plot the f3 values between 'Boncuklu' or 'Tepecik' and each of the modern populations.

You may add error bars showing 2 times the standard error, a consensus approach to describe approximate **95% confidence intervals** for each estimate.

```R
# draw a histogram
hist(f3$f3)

# plot the graphs, for both Boncuklu and Tepecik, in 1 frame
pdf("f3modern_raw.pdf", width=10, height=8)
par(mfrow=c(2,1), mgp=c(1.5, 0.5, 0), mar=c(10,6,3,6), tck=-0.01)
for (myPop in c("Boncuklu", "Tepecik")) {
  x <- f3[f3$B == myPop,]
  # order the values to improve visualisation
  x <- x[order(x$f3),]
  plot(1:nrow(x), x$f3, xlab="", ylab="f3", pch=19, col="tomato", frame.plot = T, xaxt="n", main=myPop)
  # add population names across comparisons
  axis(side = 1, labels = as.character(x$A), at = 1:nrow(x), las=2)
  # add standard errors
  segments(1:nrow(x), (x[,"f3"] - 2*x[,"se"]), 1:nrow(x), (x[,"f3"] + 2*x[,"se"]), col="tomato")
}
dev.off()
```

#### 1.1.2. Visualization in R: Regression

You may also plot the two sets of f3 values against each other, to detect any populations showing higher affinity to one vs. the other:

```R
pdf("f3reg.pdf")
f3x <- f3[f3$B == "Boncuklu",][order(unique(f3$A)),]$f3
f3y <- f3[f3$B == "Tepecik",][order(unique(f3$A)),]$f3
plot(f3x, f3y, xlab = "Boncuklu", ylab="Tepecik", pch=20, col="tomato")
abline(lm(f3y~f3x), col="tomato")
text(f3x, f3y, unique(f3$A), pos = 4)
dev.off()
```

#### 1.1.2. Visualization in R: Map display

You may also plot the f3 results on a map to better reflect geographic clustering. 

We do this for Boncuklu. For this we'll use a file with the geographic location of west Eurasian individuals in the Human Origins dataset.

```R
install.packages(c("rworldxtra", "rworldmap"))
library(rworldxtra)
library(rworldmap)

# geographic location dataset
ho3 <- as.matrix(read.table("ho3.csv", head=T, sep=",", comment.char = "", quote="", fill=T))
dim(ho3)
head(ho3)
pops <- unique(ho3[,"PopulationID.In_our_dataset."])
intersect(pops, as.character(f3$A))
lat <- as.numeric(ho3[,"Latitude"])
lon <- as.numeric(ho3[,"Longitude"])
pop <- ho3[,"PopulationID.In_our_dataset."]
ho3_v2 <- unique(data.frame(cbind(pop, lat, lon)))
head(ho3_v2)
dim(ho3_v2)
# join the location and f3 datasets
ho3_f3 <- merge(ho3_v2, f3, by.x=1, by.y=2)
head(ho3_f3)
dim(ho3_f3)

# scale the f3 data for plotting
x <- ho3_f3[ho3_f3$B == "Boncuklu", ]
f3col <- x$f3/max(x$f3)
summary(f3col)
f3col <- 1 - f3col
summary(f3col)
f3col <- f3col/max(f3col)
summary(f3col)

a_map = getMap(resolution = "high")
pdf("f3map.pdf")
par(bg = "azure")
map1 = plot(a_map,
            xlim = c(-30, 60),
            ylim = c(20, 65),
            asp = 1.3, lwd = 0.2,
            col = "antiquewhite2" )
points(x = as.numeric(as.character(x$lon)), y = as.numeric(as.character(x$lat)), pch=19, col = rgb(1,f3col,f3col,1), cex=1.4 )
dev.off()
```


### 1.2. Ancient vs Ancient

Another approach is to calculate f3 between pairs of ancient populations.

#### 1.2.1. Run f3
```bash
cat > parf3.ancient <<EOF
genotypename:    ${data}.geno
snpname:         ${data}.snp
indivname:       ${data}.ind
popfilename:     ${fstatfolder}/f3out.ancient.poplist
EOF
```

``f3out.ancient.poplist`` contains 3 populations on each line:
- an ancient West Eurasian population,
- another ancient West Eurasian population,
- Yoruba, as outgroup

```bash
cat $fstatfolder/f3out.ancient.poplist
```

```
CHG EN Yoruba
CHG Iran_N Yoruba
CHG Levant_N Yoruba
CHG Natufian Yoruba
CHG Boncuklu Yoruba
CHG Tepecik Yoruba
EN Iran_N Yoruba
EN Levant_N Yoruba
EN Natufian Yoruba
EN Boncuklu Yoruba
EN Tepecik Yoruba
Iran_N Levant_N Yoruba
Iran_N Natufian Yoruba
Iran_N Boncuklu Yoruba
Iran_N Tepecik Yoruba
Levant_N Natufian Yoruba
Levant_N Boncuklu Yoruba
Levant_N Tepecik Yoruba
Natufian Boncuklu Yoruba
Natufian Tepecik Yoruba
Boncuklu Tepecik Yoruba
```
Here:
* CHG stands for Caucasus Hunter-Gatherer (Upper Paleolithic),
* Natufian for the Epipaleolithic Natufian culture from Levant,
* Levant_N for Neolithic Levant,
* Iran_N for Neolithic Iran,
* EN for European Neolithic.

```bash
$admixtools/qp3Pop -p parf3.ancient > f3out.ancient.log
grep result f3out.ancient.log > f3ancient.out
```


#### 1.2.2. Visualization in R: Raw estimates

We first plot the raw estimates for 4 populations: "Boncuklu", "Tepecik", "Levant_N", "Iran_N".

```R
tmp<-read.table("f3ancient.out")
head(tmp)
f3a <- tmp[,c(2,3,5,6)]
colnames(f3a) <- c("B", "A", "f3", "se")
dim(f3a)
head(f3a)
```

```
         B         A       f3       se
1 Boncuklu Abkhasian 0.159139 0.001397
2 Boncuklu    Adygei 0.157850 0.001360
3 Boncuklu  Albanian 0.165172 0.001455
4 Boncuklu  Armenian 0.159583 0.001376
5 Boncuklu    Balkar 0.157370 0.001398
6 Boncuklu    Basque 0.167222 0.001410
```

Plot these the same way, for Boncuklu, EN, Levant_N and Iran_N:

```R
pdf("f3ancient_raw.pdf", width=10, height=8)
par(mfrow=c(2,2), mgp=c(1.5, 0.5, 0), mar=c(10,6,3,6), tck=-0.01)
for (myPop in c("Boncuklu", "EN", "Levant_N", "Iran_N")) {
  x <- f3a[f3a$B == myPop | f3a$A == myPop,]  # because the order of the 2 pops in the columns A and B is arbitrary, we need to check both possibilities
  x <- x[order(x$f3),]
  otherPops <- apply(x[,1:2], 1, function(z) setdiff(z, myPop))  # at each line, choose the name of the population that is NOT myPop
  plot(1:nrow(x), x$f3, xlab="", ylab="f3", pch=19, col="tomato", frame.plot = T, xaxt="n", main=myPop)
  axis(side = 1, labels = otherPops, at = 1:nrow(x), las=2)
  segments(1:nrow(x), (x[,"f3"] - 3*x[,"se"]), 1:nrow(x), (x[,"f3"] + 3*x[,"se"]), col="tomato")
}
dev.off()
```

#### 1.2.3. Visualization in R: Multidimensional scaling

Next we perform a summary of genetic distances among populations using the pairwise f3 estimates.

For this we first need to fill in a **distance matrix** from the f3 matrix. We can convert the f3 values into a distance measure by subtracting it from 1.

We can then use **multidimensional scaling** to summarise the distance matrix into 2 dimensions.

```R
pops <- unique(c(as.matrix(f3a[,1:2])))
pops

# define the empty distance matrix
f3dist <- matrix(,length(pops),length(pops))
colnames(f3dist) <- pops
rownames(f3dist) <- pops

# fill in the distance matrix
for (i in pops) for (j in pops) {
  if (i == j) f3dist[i,j] <- 0
  else {
    f3val <- (f3a$f3[(f3a$B == i & f3a$A == j) | (f3a$A == i & f3a$B == j)])
    f3dist[i,j] <- (1 - f3val)
  }
}
round(f3dist,3)
```
```
           CHG    EN Iran_N Levant_N Natufian Boncuklu Tepecik
CHG      0.000 0.846  0.841    0.857    0.865    0.845   0.844
EN       0.846 0.000  0.858    0.834    0.847    0.824   0.826
Iran_N   0.841 0.858  0.000    0.862    0.872    0.856   0.853
Levant_N 0.857 0.834  0.862    0.000    0.828    0.839   0.833
Natufian 0.865 0.847  0.872    0.828    0.000    0.849   0.849
Boncuklu 0.845 0.824  0.856    0.839    0.849    0.000   0.823
Tepecik  0.844 0.826  0.853    0.833    0.849    0.823   0.000
```

Now summarise this using MDS:
```R
# calculate MDS
fit <- cmdscale(f3dist, eig=TRUE, k=2)

# plot
fit <- cmdscale(f3dist, eig=TRUE, k=2)
cols <- c("slateblue1", "chocolate1", "slateblue1", "springgreen4", "springgreen4", "chocolate1", "chocolate1")
pdf("f3ancient_mds.pdf", width=6, height=6)
plot(fit$points[,1:2],
    xlab="Coordinate 1", ylab="Coordinate 2",
    main="MDS on f3", xlim=c(min(fit$points[,1]), max(fit$points[,1])+0.1),
    pch=19, col=cols, cex=2)
text(fit$points[,1:2], labels = rownames(f3dist), cex=1, adj = -0.05)
dev.off()
```

Note that we lack any measure of confidence in this clustering. Hence, to claim that any set of populations form a cluster -are more related to each other than to others- we need to perform formal statistical tests.


## 2. f4 statistics

The f4 statistic is calculated for 4 populations:

**f4 (A, B; C, D)**,

and measures the treeness in the form of ((A,B),(C,D)). If A and B consistently cluster together to the exclusion of the C and D, the expectation of f4 will be 0.

If A or B are, for a significant part of their genome, closer to C or D than to each other, then f4 will be non-0. The ratio is calculated as:
https://bioone.org/ContentImages/Journals/hbio/89/1/humanbiology.89.1.02/graphic/e07_21.gif

Assume that D is an outgroup, while C is a population whose lineage that could have admixed with the lineages A or B. In this setting, an f4 value indistinguishable from 0 would indicate lack of evidence for gene flow. On the contrary, an f4 value significantly different from 0 would be compatible with gene flow from C to A or B.

The f4 is thus directly related to the D-statistic used for testing admixture.

The confidence is calculated using block jackknife:
https://i.stack.imgur.com/FfXIT.jpg

The variance estimated in the jackknife, just like bootstrap analysis, indicates how strongly do different loci across the genome (which in practice represent different sources of ancestry) support an estimate. We thus calculate a confidence.

#### 2.1. Run f4

```bash
cat > parf4stat <<EOF
genotypename:    ${data}.geno
snpname:         ${data}.snp
indivname:       ${data}.ind
popfilename:     ${fstatfolder}/f4.poplist
f4mode:          YES
printsd:         YES
EOF
```

``f4.poplist`` contains 4 populations on each line: pop1, pop2, pop3, pop4, where pop1 is outgroup.

We are testing if pop3 or pop4 may have received gene flow from pop2.

```bash
head $fstatfolder/f4.poplist
```

```
Yoruba Iran_N CHG EN
Yoruba Levant_N CHG EN
Yoruba Natufian CHG EN
Yoruba Boncuklu CHG EN
Yoruba Tepecik CHG EN
Yoruba EN CHG Iran_N
Yoruba Levant_N CHG Iran_N
Yoruba Natufian CHG Iran_N
Yoruba Boncuklu CHG Iran_N
Yoruba Tepecik CHG Iran_N
```

```bash
head $fstatfolder/f4.poplist
```

```bash
wc -l $fstatfolder/f4.poplist
```

```bash
$admixtools/qpDstat -p parf4stat > f4.log
grep result f4.log > f4.out
wc -l f4.out
```
```
105
```

#### 2.2. Hypothesis testing in R

We read in the table, and remove the outgroup and other columns (for simplicity).

```R
tmp<-read.table("f4.out")
head(tmp)
f4 <- tmp[,c(3,4,5,6,7,8,11)]
colnames(f4) <- c("X", "B", "A", "f4", "se", "Z", "snps")
dim(f4)
head(f4)
```
```
         X   B      A        f4       se       Z   snps
1   Iran_N CHG     EN -0.004192 0.000364 -11.513 562261
2 Levant_N CHG     EN  0.005724 0.000441  12.985 433391
3 Natufian CHG     EN  0.004451 0.000502   8.875 251187
4 Boncuklu CHG     EN  0.005221 0.000482  10.838 560372
5  Tepecik CHG     EN  0.004591 0.000409  11.238 396894
6       EN CHG Iran_N -0.002922 0.000385  -7.590 562261
```
Here, popX represents the same thing as popC in the above description. f4 values indicate:
* if positive, higher allele sharing between popX and popA,
* if negative, higher allele sharing between popX and popB.

Note that the order of the populations in the columns A and B is arbitrary. For example on row 1, we could have popB as EN and popA as CHG. In that case, the f4 value would be positive: 0.004192.

The Z value is the standard normal variate (mean/se), from which one may calculate a p-value based on the standard normal distribution. An absolute Z-value of 2 approximately corresponds to p<0.05, and an absolute Z-value of 3 to p<0.01.

Let us calculate p-values:
```R
pnorm(2, mean=0, sd=1, lower.tail=F)*2
pnorm(3, mean=0, sd=1, lower.tail=F)*2

# add a p-value column
f4$p <- pnorm(abs(f4$Z), mean=0, sd=1, lower.tail=F)*2
f4$p
hist(f4$p)
sum(f4$p < 0.05)
f4[f4$p < 0.05,]
```

Note, however, that here we have conducted **multiple tests** tests here, and at p<0.05, we expect 5 false positives among 100 tests, even if no true positive existed.

Because of this effect, the community currently prefers to use more stringent cutoffs, such as |Z|>3, instead of |Z|>2. Alternatively, you can perform multiple testing correction, and focus on values at p<0.05 after that.

```R
f4$p <- p.adjust(f4$p)
sum(f4$p < 0.05)
```


#### 2.3. Visualization in R

We may now plot the values. We will add 3*SE error bars, to indicate what would be considered approximately significant at nominal p<0.01 cutoff. We will also indicate significance after multiple testing correction with color.

In the first type of plot, we keep popX fixed in f4(popA,popB;popX,popO), and vary popA and popB:

```R
# plot the values
pdf("f4ancient_v1.pdf", width=10, height=8)
par(mfrow=c(2,2), mgp=c(1.5, 0.5, 0), mar=c(3,6,3,6), tck=0.01)
# for the chosen pops
for (myPop in c("Boncuklu", "Tepecik", "Levant_N", "Iran_N")) {
    x <- f4[f4$X == myPop,]
    # choose the subset where f4 values are negative
    x2 <- x[x$f4 < 0, ]
    # and reverse both the values and population order
    x2[,2:3] <- x2[,3:2]
    x2$f4 <- -1*x2$f4
    x2$Z <- -1*x2$Z
    x3 <- x[x[,"f4"] > 0, ]
    # bind the positive and negative f4 parts (but we need to use matrix objects, because data.frames complicate this operation)
    x <- data.frame(rbind(as.matrix(x2), as.matrix(x3)))
    # convert factor (data.frame columns) to numeric again
    x$f4 <- as.numeric(as.character(x$f4))
    x$Z <- as.numeric(as.character(x$Z))
    x$se <- as.numeric(as.character(x$se))
    x$p <- as.numeric(as.character(x$p))
    # order the values to improve visualisation
    x <- x[order(x$f4),]
	# create a color vector describing significance
    cols <- ifelse(x$p < 0.05, "tomato", "pink")
    # plot
    plot(x$f4, 1:nrow(x), yaxt = "n", ylab="", xlab="f4", xlim=c(-0.01,0.01), pch=19, col=cols)
    axis(side = 2, labels = as.character(x$A), at = 1:nrow(x), las=2)
    axis(side = 4, labels = as.character(x$B), at = 1:nrow(x), las=2)
    abline(v=0, lty=3)
    axis(side = 3, labels = myPop, at = 0, tick = F, font = 2)
    segments((x$f4 - 3*x$se), 1:nrow(x), (x$f4 + 3*x$se), 1:nrow(x) , col=cols, lwd=2)
}
dev.off()
```

In the second type of plot, we keep popA and popB fixed in f4(popA,popB;popX,popO), and vary popX:

```R
pdf("f4ancient_v2.pdf", width=8, height=7)
par(mfrow=c(3,1), mgp=c(1.5, 0.5, 0), mar=c(3,6,3,6), tck=0.01)
for (i in 1:3) {
    if (i == 1) { POP1 <- "Boncuklu"; POP2 <- "Levant_N"}
    if (i == 2) { POP1 <- "Boncuklu"; POP2 <- "Iran_N"}
    if (i == 3) { POP1 <- "Boncuklu"; POP2 <- "Tepecik"}

    x <- f4[(f4$A == POP1 & f4$B == POP2) | (f4$A == POP2 & f4$B == POP1),]
    x <- x[order(x$f4),]

    cols <- ifelse(x$p < 0.05, "tomato", "pink")
    plot(x$f4, 1:nrow(x), yaxt = "n", ylab="", xlab="f4", xlim=c(-0.01,0.01), pch=19, col=cols, frame.plot = F)
    axis(side = 2, labels = as.character(x[,"B"]), at = 1:nrow(x), las=2)
    axis(side = 4, labels = as.character(x[,"A"]), at = 1:nrow(x), las=2)
    abline(v=0, lty=3)
    segments((x$f4 - 3*x$se), 1:nrow(x), (x$f4 + 3*x$se), 1:nrow(x) , col=cols, lwd=2)
    text((x$f4 - 3*x$se), 1:nrow(x), labels = x$X, pos = 2)
}
dev.off()

```
