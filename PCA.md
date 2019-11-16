# PCA

## Getting Started

* [Patterson et al 2006](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020190)

## Prerequisites

The dataset has 2068 present-day individuals and 21 ancient individuals, and  616 600 markers.

```bash
# dataset
data=/mnt/NAS/workshop/dataset/ho

#softwares
smartpca=/usr/local/sw/EIG-7.2.0/bin/smartpca
```

## Pipeline

### 1. Prepare parameter file

You will need calculate.pops file

```bash
cat > calculate.pops << EOF
Abkhasian
Adygei
Albanian
Armenian
Balkar
Basque
BedouinA
BedouinB
Belarusian
Bulgarian
Canary_Islander
Chechen
Chuvash
Croatian
Cypriot
Czech
Druze
English
Estonian
Finnish
French
Georgian
Greek
Hungarian
Icelandic
Iranian
Italian_North
Italian_South
Jordanian
Kumyk
Lebanese
Lezgin
Lithuanian
Maltese
Mordovian
North_Ossetian
Norwegian
Orcadian
Palestinian
Polish
Russian
Sardinian
Saudi
Scottish
Sicilian
Spanish
Syrian
Turkish
Ukrainian
EOF
```

```bash
cat > par.pca << EOF
DIR:	/mnt/NAS/workshop/dataset
genotypename:    DIR/ho.geno
snpname:         DIR/ho.snp
indivname:       DIR/ho.ind
evecoutname:     hoWE.evec
evaloutname:     hoWE.eval
poplistname:     calculate.pops
lsqproject:      YES
EOF
```

### 2. Run PCA

```bash
cd
mkdir -p pca
cd pca
$smartpca -p par.pca &> pca.log &
```

### 3. Visualization in R

```{r}
setwd("pca/")

#install.packages("ggplot2")
#install.packages("plotly")
#install.packages("tidyverse")

library(ggplot2)
#library(plotly)
library(tidyverse)

pcadat <- read.table("hoWE.evec",comment.char = "#", sep = "")

modern <- scan("calculate.pops",character())

ancient <- c("Levant_N","Natufian","CHG","Iran_N",
             "EN","WHG","Boncuklu","Tepecik")

moderndat <- subset(pcadat, V12 %in% modern)
ancdat <- subset(pcadat, V12 %in% ancient)

p <- ggplot()+
  geom_point(data = moderndat,aes(V2,-V3, label = V12), color = "grey", size = 1) +
  geom_point(data = ancdat, aes(V2,-V3,color = V12), size = 3) +
  xlab("PC1")+
  ylab("PC2")+
  theme_bw()
p


#q <- ggplotly(p)
#q
```
