## Temp repo for interview

Shaopeng Liu (sml6467@psu.edu)

Description: see email.

</br>

### ENV setup

---

```
git clone https://github.com/ShaopengLiu1/temp_interview_Fulgent.git
cd temp_interview_Fulgent
mamba env create -f src/conda_env.yml
conda activate bam_cov
```

</br>

### Data download

---

The full bam file is too large (~10G) for demo purpose. So I used chr11 here. By using `samtools view`, we can find it's paired-end data. 

```
mkdir -p result
cd result

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/exome_alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.exome.20121211.bam
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/exome_alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.exome.20121211.bam.bai
```

</br>

### Result 1: coverage by bed regions 

---

Return a region-based coverage, which is usually used for visualization in genome browser.

| File                                                | Contents                                                     |
| --------------------------------------------------- | ------------------------------------------------------------ |
| NA12878.chrom11.ILLUMINA.bwa.CEU.exome.20121211.bam | Downloaded alignment file in bam format                      |
| chr11_gene_regions.bed                              | Chr11 gene regions extracted from Gencode                    |
| **output_coverage_by_bed_regions.tsv**              | Output file with 6 columns:<br />chrom, start, end, sum_per_base_cov, total_reads, avg_cov |



```
# download chr11 gene regions
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
less gencode.v46.annotation.gtf.gz | grep "^chr11" | awk '($3=="gene") {print $1"\t"$4"\t"$5}' > chr11_gene_regions.bed

# the bam file uses "11" instead of "chr11" for annotation
sed -i 's/^chr//g' chr11_gene_regions.bed

# generate coverage
### bedcov docs here: https://www.htslib.org/doc/samtools-bedcov.html
samtools bedcov -c -Q 10 chr11_gene_regions.bed NA12878.chrom11.ILLUMINA.bwa.CEU.exome.20121211.bam > temp
### get region avg coverage
awk  -F "\t" '{$6 = sprintf("%4.2f", $4/($3-$2))}1' OFS="\t" temp > temp2
### add header
cat <(echo -e "chr\tstart\tend\tsum_per_base_cov\ttotal_reads\tavg_depth") temp2 > output_coverage_by_bed_regions.tsv
### clean file
rm temp temp2 gencode.v46.annotation.gtf.gz
```

</br>

### Result 2: per base coverage

---

Return a whole genome base-level coverage. This algorithm is almost identical to the `mpileup` methods in `samtools bedcov`.

| File                              | Contents                                               |
| --------------------------------- | ------------------------------------------------------ |
| filtered.bam                      | Remove alignements whose mapQ <10                      |
| filtered.bed                      | Extracted alignments of reads whose mapQ >= 10         |
| top10k_reads.bed                  | Top 10k reads (for faster demo purpose)                |
| **chr11_base_level_coverage.npy** | Numpy array of 136M indicating the base-level coverage |



1. A chrom size file would be required for whole genome analysis. For simplicity, I manually use the chr11 size of approximately 136M

```
# bam to bed
samtools view -bSq 10 NA12878.chrom11.ILLUMINA.bwa.CEU.exome.20121211.bam > filtered.bam
bamToBed -i filtered.bam | cut -f 1-3 > filtered.bed

# run on smaller data
head -10000 filtered.bed > top10k_reads.bed
python ../src/get_base_cov_on_chr11.py top10k_reads.bed
```

