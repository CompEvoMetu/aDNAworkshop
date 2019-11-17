# Admixture: model-based ancestry estimation

[ADMIXTURE](https://genome.cshlp.org/content/19/9/1655) is a software tool for maximum likelihood estimation of individual ancestries from multilocus SNP genotype datasets.

## Getting Started

* [Admixture-manual](http://software.genetics.ucla.edu/admixture/admixture-manual.pdf)

## Prerequisites

```bash
# dataset
data=/mnt/NAS/workshop/dataset/ho

#softwares
admixture=NA
plink=/usr/local/sw/plink-1.90/plink
pong=NA
R-studio
```

The dataset has 2068 present-day individuals  and 21 ancient individuals, and  616 600 autosomal markers.

Genotypes of Affymetrix Human Origins present-day individuals downloaded from [here](https://reich.hms.harvard.edu/datasets) ([Lazaridis et al 2014](https://www.nature.com/articles/nature13673), [Lazaridis et al 2016](https://www.nature.com/articles/nature19310))


## Pipeline

### 1. Prepare dataset for admixture

```bash

cd
mkdir -p admixture

# To keep only ancient individuals in a file and remove all loci where more than 99.9% of genotypes are missing.
tail -n +2069 ${data}.pedind > ancient.indv
$plink --file ${data} --keep ancient.indv --geno 0.99 --make-bed --out ancient.geno99

# To keep only modern individuals (reference panel) in a file and remove variants where more than 99.9% of genotypes are missing in ancient data
head -n 2068 ${data}.pedind > modern.indv
awk '{print $2}' ancient.geno99.bim > snps
$plink --file ${data} --keep modern.indv --extract snps --make-bed --out modern

# Prune data
# using the --indep-pairwise option: removal each SNP that has an R2 value of greater than 0.4 with any other SNP within a 200-SNP sliding window (advanced by 25 SNPs each time).
$plink --bfile modern --indep-pairwise 200 25 0.4
$plink --bfile modern --extract plink.prune.in --make-bed --out modern.admixsnps
$plink --bfile ancient.geno99 --extract plink.prune.in --make-bed --out ancient.admixsnps
```

### 2. Run admixture

Assume the individuals derive their ancestry from 10 ancestral population:

```bash
mkdir -p k10/run1
cd  k10/run1
admixture -s $RANDOM /mnt/NAS/workshop/admixture/modern.admixsnps.bed 10
```
**Two output:**

> 1. modern.admixsnps.10.Q (the ancestry fractions)

```bash
head /mnt/NAS/workshop/admixture/k10/run1/modern.admixsnps.10.Q
```

```
0.158028 0.672972 0.002553 0.025055 0.089839 0.000010 0.035225 0.015474 0.000833 0.000010
0.000010 0.908173 0.000010 0.000010 0.038680 0.000010 0.053078 0.000010 0.000010 0.000010
0.000010 0.999910 0.000010 0.000010 0.000010 0.000010 0.000010 0.000010 0.000010 0.000010
0.000010 0.727451 0.000010 0.072127 0.093945 0.000869 0.050367 0.047225 0.000010 0.007985
0.093496 0.759203 0.000010 0.009948 0.036727 0.000010 0.054153 0.021743 0.008708 0.016003
0.068244 0.890481 0.000010 0.000010 0.004674 0.000010 0.036540 0.000010 0.000010 0.000010
0.000010 0.728189 0.000010 0.214340 0.000010 0.000011 0.040328 0.013217 0.003875 0.000010
0.000010 0.673473 0.005429 0.116532 0.000010 0.000010 0.094209 0.098506 0.000010 0.011812
0.000010 0.000010 0.000017 0.040367 0.000010 0.000010 0.000010 0.000010 0.701821 0.257735
0.000010 0.000010 0.000010 0.028053 0.000010 0.000010 0.000010 0.000010 0.729116 0.242761
```

> 2. modern.admixsnps.10.P (the allele frequencies of the inferred ancestral populations)


```bash
head /mnt/NAS/workshop/admixture/k10/run1/modern.admixsnps.10.P
```

```
0.330951 0.641307 0.795701 0.871557 0.518768 0.927323 0.769164 0.803492 0.630176 0.925985
0.876535 0.602634 0.820633 0.835990 0.727728 0.940231 0.795822 0.778324 0.453960 0.887941
0.704559 0.465145 0.707353 0.964257 0.804436 0.075585 0.904571 0.911557 0.590574 0.793412
0.266606 0.499776 0.643588 0.775814 0.602188 0.175487 0.695923 0.720671 0.621338 0.739245
0.461696 0.326021 0.999990 0.962489 0.405563 0.999990 0.855227 0.913788 0.962554 0.950173
0.462990 0.781608 0.867691 0.549678 0.400188 0.796729 0.591254 0.515443 0.570955 0.835824
0.999990 0.999990 0.999990 0.975038 0.999990 0.999990 0.992974 0.999989 0.999990 0.999990
0.999990 0.999990 0.831823 0.903551 0.989368 0.807794 0.843045 0.920314 0.610000 0.886243
0.636211 0.430156 0.798091 0.613824 0.656324 0.759256 0.560640 0.636532 0.603132 0.837607
0.502172 0.671637 0.837731 0.977087 0.925560 0.645656 0.906791 0.855188 0.932370 0.971443
```

**Projection:** estimating ancient individual ancestry is to “project” the new samples on to the population structure (allele frequencies) learned from the reference panels.

```bash
cp modern.admixsnps.10.P ancient.admixsnps.10.P.in
admixture -P /mnt/NAS/workshop/admixture/ancient.admixsnps.bed 10
```

```bash
cat /mnt/NAS/workshop/admixture/k10/run1/modern.admixsnps.10.Q /mnt/NAS/workshop/admixture/k10/run1/ancient.admixsnps.10.Q > k10.run1.all.Q
```
> How do I choose the correct value for K?

Use ADMIXTURE’s cross-validation procedure


### 3. Visualization

Reading recommendation: [Lawson et al 2018](https://www.nature.com/articles/s41467-018-05257-7)

#### 3.1 pong

* [Pong-manual](http://brown.edu/Research/Ramachandran_Lab/files/pong/pong-manual.pdf)
* [Pong-github](https://github.com/ramachandran-lab/pong)
* [Behr et al. 2016, Bioinformatics](https://academic.oup.com/bioinformatics/article/32/18/2817/1744074)

> Input to pong

#### 3.1.1. filemap file (required)

The filemap must be a three-column, tab-delimited file. Each Q matrix is represented by a three-column line in
the filemap with the following format:

```bash
cat > pong_filemap << EOF
k10r1 10  /mnt/NAS/workshop/admixture/k10/run1/k10.run1.all.Q
k10r2 10  /mnt/NAS/workshop/admixture/k10/run2/k10.run2.all.Q
k10r3 10  /mnt/NAS/workshop/admixture/k10/run3/k10.run3.all.Q
EOF
```

```bash
#to be sure filemap file is tab-delimited
awk '{print $1, "\t", $2, "\t", $3}' pong_filemap > pong_filemap2
```

#### 3.1.2. Population labels (optional, but recommended)

**Version 1**:

```bash
awk '{print $3}' /mnt/NAS/workshop/dataset/ho.ind > ind2pop.txt
```

```bash
head ind2pop.txt
```

```
Khomani
Khomani
Khomani
Khomani
Khomani
Khomani
Khomani
Khomani
Yukagir
Yukagir
```

```bash
tail ind2pop.txt
```

```
WHG
Boncuklu
Boncuklu
Boncuklu
Boncuklu
Tepecik
Tepecik
Tepecik
Tepecik
Tepecik
```

```bash
pong -m pong_filemap2 -i ind2pop.txt
```

**Version 2**

```bash
head /mnt/NAS/workshop/admixture/ind2pop_metaPops.txt
```

```
SouthAfrica
SouthAfrica
SouthAfrica
SouthAfrica
SouthAfrica
SouthAfrica
SouthAfrica
SouthAfrica
CentralAsiaSiberia_SIB
CentralAsiaSiberia_SIB
```

```bash
tail /mnt/NAS/workshop/admixture/ind2pop_metaPops.txt
```

```
European
Anatolian
Anatolian
Anatolian
Anatolian
Anatolian
Anatolian
Anatolian
Anatolian
```

```bash
pong -m pong_filemap2 -i /mnt/NAS/workshop/admixture/ind2pop_metaPops.txt
```
#### 3.1.3. Population order and detailed names (optional)

```bash
cat > pop_order.txt << EOF
EOF
```

#### 3.1 R

```R
adm_modern<-read.table("/mnt/NAS/workshop/admixture/k10/run1/modern.admixsnps.10.Q", header = F)
head(adm_modern)
popInfo<-read.table("/mnt/NAS/workshop/admixture/ind2pop.txt", stringsAsFactors = F, header = F) # including both present-day and ancient samples
metapopInfo<-read.table("/mnt/NAS/workshop/admixture/ind2pop_metaPops.txt", stringsAsFactors = F, header = F) # including both present-day and ancient samples


modernPopInfo<-cbind(as.character(popInfo[1:2068,]), metapopInfo[1:2068,])
colnames(modernPopInfo)<-c("pop", "metapop")
head(modernPopInfo)


adm_modern <- cbind(modernPopInfo, adm_modern )
colnames(adm_modern) <- c("pop", "metapop", "k1","k2","k3","k4","k5","k6","k7","k8","k9","k10")
head(adm_modern)

## ---------------------------------------------------------------

# install.packages("RColorBrewer")
library("RColorBrewer")
colors = brewer.pal(n = 10, name = "Paired")
colors

# How to not draw
barplot(t(as.matrix(adm_modern[,3:12])), col = colors, border = NA, cex.names=0.7,
        names.arg=adm_modern$pop, las=2)


## ---------------------------------------------------------------

# Dirty solution
pop_order <- t(sapply(unique(adm_modern$metapop), function(x){
  apply(adm_modern[adm_modern$metapop == x, 3:12], 2, median)
}))

pop_order=data.frame(pop_order)
colnames(pop_order)<-c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10")
rownames(pop_order)<-unique(adm_modern$metapop)
head(pop_order)

pop_ordered <- pop_order[with(pop_order, order(-pop_order$k1,
                                               -pop_order$k2,
                                               -pop_order$k3,
                                               -pop_order$k4,
                                               -pop_order$k5,
                                               -pop_order$k6,
                                               -pop_order$k7,
                                               -pop_order$k8,
                                               -pop_order$k9,
                                               -pop_order$k10)), ]

adm_modern$metapop <- factor(adm_modern$metapop, levels=rownames(pop_ordered))
ordered <- adm_modern[with(adm_modern, order(adm_modern$metapop,
                                             -adm_modern$k1,
                                             -adm_modern$k2,
                                             -adm_modern$k3,
                                             -adm_modern$k4,
                                             -adm_modern$k5,
                                             -adm_modern$k6,
                                             -adm_modern$k7,
                                             -adm_modern$k8,
                                             -adm_modern$k9,
                                             -adm_modern$k10)), ]


# Let's see what happened
View(ordered[1:30,])


## ---------------------------------------------------------------

pdf("admixture_modern.k10.pdf", width = 50, height = 5)
par(mar=c(10,5,5,5))
a <- barplot(t(as.matrix(ordered[,3:12])), col = colors, border = NA, cex.names=0.2,
             xaxt='n', las=2)
x <-  match(unique(ordered$metapop), ordered$metapop)
axis(side=1, at = a[x], labels =  ordered$metapop[x], las=2)

dev.off()

## ---------------------------------------------------------------

# adding ancient samples

adm_ancient<-read.table("/mnt/NAS/workshop/admixture/k10/run1/ancient.admixsnps.10.Q", header = F)

ancientPopInfo<-cbind(as.character(popInfo[2069:nrow(popInfo),]), metapopInfo[2069:nrow(popInfo),])
colnames(ancientPopInfo)<-c("pop", "metapop")
head(ancientPopInfo)


adm_ancient <- cbind(ancientPopInfo, adm_ancient )
colnames(adm_ancient) <- c("pop", "metapop", "k1","k2","k3","k4","k5","k6","k7","k8","k9","k10")
head(adm_ancient)

pdf("admixture_ancient.k10.pdf", width = 12, height = 6)
par(mar=c(10,5,5,5))
barplot(t(as.matrix(adm_ancient[,3:12])), col = colors, border = NA, cex.names=0.7,
        names.arg=adm_ancient$pop, las=2)
dev.off()

## ---------------------------------------------------------------

```

**NOTE:** For visualization you can also check [CLUMPP](https://web.stanford.edu/group/rosenberglab/clumpp.html) and [distruct](http://web.stanford.edu/group/rosenberglab/distruct.html)
