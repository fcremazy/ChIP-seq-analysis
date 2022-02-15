# ChIP-seq-analysis

## 1. Trimming/QC des reads.

Je te conseille d'utiliser *fastp* pour ça.

La ligne de commande:
`fastp -A -h prefix_qc.html -i prefix_reads1.fastq.gz -I prefix_reads2.fastq.gz`

N.B.: par *prefix*, j'entends la partie spécifique à ton experience du nom de ton fichier. 

**Lien:** https://github.com/OpenGene/fastp

## 2. Alignement sur génome de référence

Comme souvent il y a différentes alternatives mais *bowtie* reste le plus utilisé.
La commande que j'utilise pour nos analyses est la suivante:

```
bowtie2 \
  --sensitive \
  -x nom_index_genome \
  -1 prefix_reads1.fastq.gz \
  -2 prefix_reads2.fastq.gz \
  -S prefix_.sam \
  -p 8
  2> prefix.log
```

N.B.:
- le caractère `\` permet de continuer la commande sur la ligne suivante. Tu peux aussi mettre tout à la suite sur la même ligne.
- la partie `2> prefix.log` permet de rediriger toutes les informations affichées à l'écran durant l'executions des différents programmes que tu vas utiliser vers un fichier "*log*" que tu pourras consulter tranquillement plus tard s'il y a un problème. Par la suite, il faudra utiliser `2>> prefix.log` par la suite afin que les nouvelles infos soient ajoutées à la suite du fichier déjà existant.

**Lien:** http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

## 3. Conversion de Sam en Bam et filtrages de l'alignement

### Converting Sam to Bam

`samtools view -bSo $FULL_PATH_RES/$PREFIX.bam $FULL_PATH_RES/$PREFIX.sam`

### Filtrage unpaired and unique reads

```
samtools view -b -f 2 -F 268 \
  prefix.bam \
  -o prefix.pairs.unique.bam \
  2>> prefix.log
```

### Sorting Bam file

```
samtools sort -T prefix.pairs.unique.bam \
  > prefix.sorted.pairs.unique.bam \
  2>> prefix.log
```

N.B.: Le *sorting* (tri) du fichier Bam peut être nécessaire avant d'utiliser certains programmes.

### Keeping the "good" chromosomes

```
samtools index prefix.sorted.pairs.unique.bam
samtools idxstats prefix.sorted.pairs.unique.bam |gawk '{if($1!="*"){print}}'|cut -f 1| \
   grep "chr"|grep -v "chrUn\|chrM\|random"| \
   xargs samtools view -b prefix.sorted.pairs.unique.bam \
   > $FULL_PATH_RES/prefix.noMTorRandom.pairs.unique.bam
   2>> prefix.log
```

## 4. Filtrage(s) des "duplicates"

Comme je te le disais hier, je fais ça en deux fois pour plus d'efficacité.

### First duplicate removal step using samtools

```
samtools collate -o prefix.namecollate.noMTorRandom.pairs.unique.bam prefix.noMTorRandom.pairs.unique.bam
samtools fixmate -m prefix.namecollate.noMTorRandom.pairs.unique.bam prefix.noMTorRandom.pairs.unique.mrk1.bam
samtools sort -o prefix.noMTorRandom.pairs.unique.sorted.mrk1.bam prefix.noMTorRandom.pairs.unique.mrk1.bam 2>/dev/null
samtools markdup -r -s prefix.noMTorRandom.pairs.unique.sorted.mrk1.bam prefix.noDup1.noMTorRandom.pairs.unique.bam 2>> prefix.log
```

### Marking duplicates with Picard

**ATTENTION:** il va falloir que nous installions *Picard* dans ton répertoire *bin* local pour que ça fonctionne!

```
java -jar $PICARD_FOLDER/picard.jar MarkDuplicates \
	I=prefix.noDup1.noMTorRandom.pairs.unique.bam \
	O=prefix.noDup1.noMTorRandom.pairs.unique.mrk2.bam \
	M=prefix.picard_marked_dup_metrics.txt \
	AS=TRUE 2>> prefix.log
```

### Removing Picard (or second) duplicates with Samtools

```
samtools view -b -F 0x400 \
	prefix.noDup1.noMTorRandom.pairs.unique.mrk2.bam \
	-o prefix.noDup2.noMTorRandom.pairs.unique.bam \
	2>> prefix.log
```

**Liens:**
- http://www.htslib.org/
- https://broadinstitute.github.io/picard/

## 5. Dernier filtrage: ne garder que les reads alignés de qualité "MAPQ > 30"

J'ai préféré faire ce filtrage car il sera plus facile de l'éviter si la qualité de tes reads n'est pas top mais que tu veux quand même tenter l'analyse.
Tu peux d'ailleurs chosir de baisser un peu la qualité en changeant la valeur après le paramètre `-q`.

```
samtools view -b -q 30 \
	prefix.noDup2.noMTorRandom.pairs.unique.bam \
	-o prefix.MAPQ30.noDup2.noMTorRandom.pairs.unique.bam \
	2>> prefix.log
```

## 6. Renommage et Indexage du Bam final
 
```
mv prefix.MAPQ30.noDup2.noMTorRandom.pairs.unique.bam prefix.analyzed.bam
samtools index prefix.analyzed.bam 2>> prefix.log
```

## 7. Génération d'un fichier graph de type *bigWig* (plus léger que le Bam)

Différents outils peuvent te permettre ça, mais perso je te conseille d'utiliser les *deepTools* que j'aprrécie particulièrement.
Les paramètres comme `--binSize` et `--smoothLength' te permetttront de jouer sur la résolution du graphe. Je te laisse consulter la doc des *deepTools* si tu souhaites modifier ces paramètres.

```
bamCoverage \
	--bam prefix.analyzed.bam \
	--outFileFormat bigwig \
	--outFileName prefix.analyzed.bw \
	--binSize 20 \
	--smoothLength 60 \
	--normalizeUsing RPKM \
	--effectiveGenomeSize 142573017\
	--numberOfProcessors 8 \
	2>> prefix.log
```

**Lien:** https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

## 8. "Peakcalling" avec MACS2

Encore une fois je te renvoie aux liens que j'ajoute à la fin de ce paragraphe pour plus de détails sur le model mis au point et les paramètres sue lesquels tu peux jouer. Il manque un paramètre `-c` après lequel il faudrait ajouter le nom du fichier *Input* que tu n'as pas pour l'instant.

Perso, je l'utilise comme ceci:

```
macs2 callpeak \
	-t prefix.analyzed.bam \
	-f BAMPE \
	-n prefix \
	-g 142573017 \
	-p 1e-7 \
	--keep-dup all \
	2>> prefix.log
```

Les fichiers résulats vont tous débuter encore une fois par "*préfix*".

**Liens:**
- https://github.com/macs3-project/MACS
- https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
