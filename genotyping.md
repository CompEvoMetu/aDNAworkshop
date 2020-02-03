**GENOTYPING**

**trimming bam biles**

```bash
## looking inside a bam bile
samtools view /mnt/NAS/workshop/genotyping/merged_bams/Tep001_031102_15.merged.hs37d5.fa.cons.90perc.bam | head
## trimming 10 bases at each end
cd
/usr/local/sw/bamUtil/bin/bam trimBAM /mnt/NAS/workshop/genotyping/merged_bams/Tep001_031102_15.merged.hs37d5.fa.cons.90perc.bam Tep001_031102_15.merged.hs37d5.fa.cons.90perc.trim.bam 10
ls
##indexing the bam file
samtools index Tep001_031102_15.merged.hs37d5.fa.cons.90perc.trim.bam
ls
## bam file after trimming
samtools view Tep001_031102_15.merged.hs37d5.fa.cons.90perc.trim.bam | head
```
**file formats**

*plink format*
```bash
less -S /mnt/NAS/workshop/dataset/ho.ped

less /mnt/NAS/workshop/dataset/ho.map
```
*eigenstrat format*
```bash
less /mnt/NAS/workshop/genotyping/modern/HumanOrigins.geno

less /mnt/NAS/workshop/genotyping/modern/HumanOrigins.snp

less /mnt/NAS/workshop/genotyping/modern/HumanOrigins.ind  



```

**converting one format to other**

```bash
less /mnt/NAS/workshop/genotyping/scripts/ped2eigen
less /mnt/NAS/workshop/genotyping/scripts/eigen2ped

```

*pedtoeigen*

AdmixTools/bin/convertf -p ped2eigen

*eigentoped*

AdmixTools/bin/convertf -p eigen2ped


**genotype calling**

```bash
/usr/local/sw/samtools-1.9/samtools mpileup -R -B -q30 -Q30 -l /mnt/NAS/workshop/genotyping/positions.bed -f /mnt/NAS/workshop/prep/hs37d5.fa Tep001_031102_15.merged.hs37d5.fa.cons.90perc.trim.bam > Tep001pileup.txt

/usr/local/sw/sequenceTools/pileupCaller --sampleNames Tep001 -f /mnt/NAS/workshop/genotyping/eigenfile.snp -o EigenStrat -e Tep001.out < Tep001pileup.txt
```

**merging ancient samples with modern**

```bash
cd

mkdir snpcall

cd snpcall

cat > /mnt/NAS/workshop/genotyping/scripts/parmerge <<EOF
geno1:            HumanOrigins.geno
snp1:             HumanOrigins.snp
ind1:             HumanOrigins.ind
geno2:            Tep001.out.geno.txt
snp2:             Tep001.out.snp.txt
ind2:             Tep001.out.ind.txt
genooutfilename:  ho.all.geno
snpoutfilename:   ho.all.snp
indoutfilename:   ho.all.ind
strandcheck:      NO
EOF

cp /mnt/NAS/workshop/genotyping/scripts/parmerge ./

/usr/local/sw/EIG-7.2.0/bin/mergeit -p parmerge
```
