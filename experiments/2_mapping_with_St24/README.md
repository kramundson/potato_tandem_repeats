Analysis of satellite polymorphism between potato genotypes Superior and Russet Burbank

## Aim

In light of results from Experiment 1, explore why exact matches weren't recovered along
the lengths of St24 for both Superior and Russet Burbank. This could be due to SNPs/indels
in the reads. To test this, align genomic reads from Superior and Russet Burbank to the
potato reference genome with St24 appended, then inspect read depth along St24 and SNPs
and/or short indels in the alignments.

Gong et al (2012) show that St24 is a tandemly arranged repeat with a long monomer size,
so I'll append a duplicated St24 monomer before mapping to catch junction-spanning reads.

## Procedure

1. Duplicate St24

```
cd /share/comailab/kramundson/bruno_move_to_share/tohru_reproduce/potato_tandem_repeats/experiments/2_mapping_with_St24/

ST24="/share/comailab/kramundson/bruno_move_to_share/tohru_reproduce/potato_tandem_repeats/experiments/1_superior_rb_comparison/data/St24.fasta"
echo $ST24

head -n 1 $ST24 > St24_duplicated.fasta && tail -n +2 $ST24 $ST24 | grep -v "St24.fasta" | tr "\n" " " | sed -e "s/ //g" >> St24_duplicated.fasta && echo "" >> St24_duplicated.fasta
```

2. Download reference genome sequences. I added trailing newlines to the chloroplast and
mitochondrion sequences for consistency with the DM1-3 assembly:

```
# DM1-3 base assembly
wget http://solanaceae.plantbiology.msu.edu/data/potato_dm_v404_all_pm_un.fasta.zip && unzip potato_dm_v404_all_pm_un.fasta.zip

# chloroplast
wget http://solanaceae.plantbiology.msu.edu/data/S_tuberosum_Group_Phureja_chloroplast_DM1-3-516-R44.fasta.zip \
    && unzip S_tuberosum_Group_Phureja_chloroplast_DM1-3-516-R44.fasta.zip \
    && echo "" >> S_tuberosum_Group_Phureja_chloroplast_DM1-3-516-R44.fasta

# mitochondrion
wget http://solanaceae.plantbiology.msu.edu/data/S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta.zip \
    && unzip S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta.zip \
    && echo "" >> S_tuberosum_Group_Phureja_mitochondrion_DM1-3-516-R44.fasta
```

3. Combine reference sequences and duplicated St24. Here, I also remove metacharacters and
shorten the FASTA headers from St24, chloroplast, and mitochondrion sequences.

```
cat *.fasta | sed 's/\?/qmark_segment/g' | sed -e 's/ .\+//g' -e 's/-/_/g' \
    > potato_ref_St24dup.fasta
```

4. Index the reference genome

```
bwa index potato_ref_St24dup.fasta
```

5. Align reads to reference. For Superior reads, use the individual SRA runs rather than
the combined output.

```
# one superior library as an example
# Using this as a template; do the rest of Superior and Russet Burbank on your own
# adapt path names as needed

# set up symbolic link to reads from experiment 1; can loop this to make all symlinks at once
ln -s /share/comailab/kramundson/bruno_move_to_share/tohru_reproduce/potato_tandem_repeats/data/reads/SRR2069941.fastq .

# aligns one library using 12 cores
bwa mem -t 12 -R '@RG\tID:SRR2069941\tSM:Superior' potato_ref_St24dup.fasta \
    -p SRR2069940.fastq 2> SRR2069941_aln.err | \
    samtools view -b -o SRR2069941.bam -
```

6. Sort bams

```
samtools sort -o SRR2069941.sort.bam SRR2069941.bam
```

7. Merge bams from Superior.

```
samtools merge superior_merged.bam SRR2069941.sort.bam SRR2069940.sort.bam SRR2069942.sort.bam SRR2070067.sort.bam SRR5349638.sort.bam
```

8. Inspect coverage of St24 in a genome browser. There are bound to be St24 sequences strewn
about chromosome 1 and perhaps the rest of the genome, so the mapping quality will be low.