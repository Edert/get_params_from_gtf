# get_params_from_gtf

**Python script to extract and calculate genomic parameters from GTF and FASTA files**

---

## Description

This script processes a GTF annotation file and a reference genome FASTA file to compute per-gene summary statistics, including:

- Mean GC content of transcripts  
- Average number of exons per gene  
- Average transcript length  
- Average gene span length  
- Mean lengths of 5' and 3' UTRs  

It outputs a tab-separated file summarizing these parameters for each gene.

---

## Requirements

- Python 3.x  
- [gffutils](https://pypi.org/project/gffutils/)  
- [Biopython](https://biopython.org/)  
- [pandas](https://pandas.pydata.org/)

Install dependencies with:

```bash
pip install gffutils biopython pandas
```

## Usage

```bash
python get_params_from_gtf.py -g input_annotations.gtf -f reference_genome.fasta -o output.tsv
```

## Arguments

| Option | Description                         | Required |
|--------|-----------------------------------|----------|
| `-g`   | Path to input GTF annotation file | Yes      |
| `-f`   | Path to reference genome FASTA file | Yes      |
| `-o`   | Path to output TSV file            | Yes      |


## Output

The output file contains one row per gene with the following columns:

- `gene_id`  
- `gene_name`  
- `mean_gc_content` (%)  
- `mean_num_exons`  
- `mean_transcript_length` (bp)  
- `mean_gene_length` (bp)  
- `mean_5utr_length` (bp)  
- `mean_3utr_length` (bp)  

## Example

```bash
python get_params_from_gtf.py -g example.gtf -f genome.fasta -o gene_parameters.tsv
```
