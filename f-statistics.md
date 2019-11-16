# f-statistics


## Getting Started

Reading materials: [Patterson et al 2012](https://www.genetics.org/content/192/3/1065) ; [Peter 2016](https://www.genetics.org/content/202/4/1485) ; [Harris and DeGiorgio](https://doi.org/10.13110/humanbiology.89.1.02)

## Prerequisites

```bash
#dataset
data=/mnt/NAS/workshop/dataset/ho

#software
admixtools=/usr/local/sw/AdmixTools/bin
```

## 1. Outgroup f3

The outgroup f3 method calculates shared genetic drift between two populations

### 1.1. Modern vs Ancient

#### 1.1.1. Run f3

```bash
cat > parf3.modern <<EOF
genotypename:    ${data}.geno
snpname:         ${data}.snp
indivname:       ${data}.ind
popfilename:     f3out.modern.poplist
EOF
```

``f3out.modern.poplist`` contains 3 populations on each line:
- an Anatolian Neolithic population, either Boncuklu or Tepecik
- a West Eurasian present-day population
- Yoruba, as outgroup



```bash
head f3out.modern.poplist
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
$admixtools/qp3Pop -p parf3.modern > f3out.modern.log
grep result f3out.modern.log > f3modern.out
```

#### 1.1.2. Visualization in R

```R
tmp<-read.table("f3modern.out")
head(tmp)
f3 <- tmp[,c(2,3,5,6)]
colnames(f3) <- c("B", "A", "f3", "se")
dim(f3)
head(f3)

par(mfrow=c(2,1), mgp=c(1.5, 0.5, 0), mar=c(10,6,3,6), tck=-0.01)
for (myPop in c("Boncuklu", "Tepecik")) {
  x <- f3[f3$B == myPop,]
  x <- x[order(x$f3),]
  # assign(paste("f3", myPop, "_"), value = x)
  plot(1:nrow(x), x$f3, xlab="", ylab="f3", pch=19, col="tomato", frame.plot = T, xaxt="n", main=myPop)
  axis(side = 1, labels = as.character(x$A), at = 1:nrow(x), las=2)
  segments(1:nrow(x), (x[,"f3"] - 3*x[,"se"]), 1:nrow(x), (x[,"f3"] + 3*x[,"se"]), col="tomato")
}

f3x <- f3[f3$B == "Boncuklu",][order(unique(f3$A)),]$f3
f3y <- f3[f3$B == "Tepecik",][order(unique(f3$A)),]$f3
plot(f3x, f3y, xlab = "Boncuklu", ylab="Tepecik")
abline(lm(f3y~f3x))
text(f3x, f3y, unique(f3$A), pos = 4)
```

### 1.2. Ancient vs Ancient

#### 1.2.1. Run f3
```bash
cat > parf3.ancient <<EOF
genotypename:    ${data}.geno
snpname:         ${data}.snp
indivname:       ${data}.ind
popfilename:     f3out.ancient.poplist
EOF
```

``f3out.ancient.poplist`` contains 3 populations on each line <Source1 (Ancient population)> <Source2 (Ancient population)> < Target (Yoruba)>

```bash
cat f3out.ancient.poplist
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

```bash
$admixtools/qp3Pop -p parf3.ancient > f3out.ancient.log
grep result f3out.ancient.log > f3ancient.out
```


#### 1.2.2. Visualization in R: Multidimensional scaling

```R
tmp<-read.table("f3ancient.out")
head(tmp)
f3a <- tmp[,c(2,3,5,6)]
colnames(f3a) <- c("B", "A", "f3", "se")
dim(f3a)
head(f3a)

par(mfrow=c(2,2), mgp=c(1.5, 0.5, 0), mar=c(10,6,3,6), tck=-0.01)
for (myPop in c("Boncuklu", "Tepecik", "Levant_N", "Iran_N")) {
  x <- f3a[f3a$B == myPop | f3a$A == myPop,]
  x <- x[order(x$f3),]
  otherPops <- apply(x[,1:2], 1, function(z) setdiff(z, myPop))
  # assign(paste("f3", myPop, "_"), value = x)
  plot(1:nrow(x), x$f3, xlab="", ylab="f3", pch=19, col="tomato", frame.plot = T, xaxt="n", main=myPop)
  axis(side = 1, labels = otherPops, at = 1:nrow(x), las=2)
  segments(1:nrow(x), (x[,"f3"] - 3*x[,"se"]), 1:nrow(x), (x[,"f3"] + 3*x[,"se"]), col="tomato")
}

pops <- unique(c(as.matrix(f3a[,1:2])))
pops
f3dist <- matrix(,length(pops),length(pops))
colnames(f3dist) <- pops
rownames(f3dist) <- pops
for (i in pops) for (j in pops) {
  if (i == j) f3dist[i,j] <- 0
  else {
    f3val <- (f3a$f3[(f3a$B == i & f3a$A == j) | (f3a$A == i & f3a$B == j)])
    f3dist[i,j] <- (1 - f3val)
  }
}
f3dist

fit <- cmdscale(f3dist, eig=TRUE, k=2)
cols <- c("slateblue1", "chocolate1", "slateblue1", "springgreen4", "springgreen4", "chocolate1", "chocolate1")
plot(fit$points[,1], fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2",
     main="MDS on f3", type="n", xlim=c(min(fit$points[,1]), max(fit$points[,1])+0.1))
points(fit$points[,1:2], pch=19, col=cols)
text(fit$points[,1:2], labels = rownames(f3dist), cex=.7, adj = -0.05)
```


## 2. f4 statistics

#### 2.1. Run f4

```bash
cat > parf4stat <<EOF
genotypename:    ${data}.geno
snpname:         ${data}.snp
indivname:       ${data}.ind
popfilename:     f4.poplist
f4mode:          YES
printsd:         YES
EOF
```

``f4.poplist`` contains 4 populations on each line <pop1> <pop2> <pop3> <pop4>) where pop1 is outgroup

```bash
head f4.poplist
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
$admixtools/qpDstat -p parf4stat > f4.log
grep result f4.log > f4.out
```

#### 2.2. Visualization in R

```R
tmp<-read.table("f4.out")
head(tmp)
f4 <- tmp[,c(3,4,5,6,7,8,11)]
colnames(f4) <- c("X", "B", "A", "f4", "se", "Z", "snps")
dim(f4)
head(f4)
unique(f4[,"X"])

pdf("~/f4_plots.pdf", width = 10)
    par(mfrow=c(2,2), mgp=c(1.5, 0.5, 0), mar=c(3,6,3,6), tck=0.01)
    for (myPop in c("Boncuklu", "Tepecik", "Iran_N", "Levant_N")) {
    x <- f4[f4$X == myPop,]
    dim(x)
    x <- x[order(x$f4),]
    x2 <- x[x$f4 < 0, ]
    x2 <- x2[,c(1,3,2,4:7)]
    x2$f4 <- -1*x2$f4
    x2$Z <- -1*x2$Z
    x3 <- x[x[,"f4"] > 0, ]
    x <- data.frame(rbind(as.matrix(x2), as.matrix(x3)))
    x$f4 <- as.numeric(as.character(x$f4))
    x$Z <- as.numeric(as.character(x$Z))
    x$se <- as.numeric(as.character(x$se))
    x <- x[order(x$f4),]
    plot(x$f4, 1:nrow(x), yaxt = "n", ylab="", xlab="f4", xlim=c(-0.01,0.01), pch=19, col="tomato")
    axis(side = 2, labels = as.character(x$A), at = 1:nrow(x), las=2)
    axis(side = 4, labels = as.character(x$B), at = 1:nrow(x), las=2)
    abline(v=0, lty=3)
    axis(side = 3, labels = myPop, at = 0, tick = F, font = 2)
    segments((x$f4 - 3*x$se), 1:nrow(x), (x$f4 + 3*x$se), 1:nrow(x) , col="tomato")
}
dev.off()
```


```R
par(mfrow=c(3,1), mgp=c(1.5, 0.5, 0), mar=c(3,6,3,6), tck=0.01)
for (i in 1:3) {
    if (i == 1) { POP1 <- "Boncuklu"; POP2 <- "Levant_N"}
    if (i == 2) { POP1 <- "Boncuklu"; POP2 <- "Iran_N"}
    if (i == 3) { POP1 <- "Boncuklu"; POP2 <- "Tepecik"}

    x <- f4[(f4$A == POP1 & f4$B == POP2) | (f4$A == POP2 & f4$B == POP1),]
    x <- x[order(x$f4),]

    plot(x$f4, 1:nrow(x), yaxt = "n", ylab="", xlab="f4", xlim=c(-0.01,0.01), pch=19, col="tomato", frame.plot = F)
    axis(side = 2, labels = as.character(x[,"B"]), at = 1:nrow(x), las=2)
    axis(side = 4, labels = as.character(x[,"A"]), at = 1:nrow(x), las=2)
    abline(v=0, lty=3)
    segments((x$f4 - 3*x$se), 1:nrow(x), (x$f4 + 3*x$se), 1:nrow(x) , col="tomato")
    text((x$f4 + 3*x$se), 1:nrow(x), labels = x$X, pos = 4)
}
```
