# ChIP-seq-analysis

## 1. Trimming/QC des reads.

Je te conseille d'utiliser *fastp* pour ça.

La ligne de commande:
`fastp -A -h prefix_qc.html -i prefix_reads1.fastq.gz -I prefix_reads2.fastq.gz`

N.B.: par *prefix*, j'entends la partie spécifique à ton experience du nom de ton fichier. 

Link: https://github.com/OpenGene/fastp

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

N.B.: tu constateras que cette opération sera effectuée ponctuellement durant l'analyse. En effet, le *sorting* (tri) du fichier Bam peut-être nécessaire avant d'utiliser certains programmes.

