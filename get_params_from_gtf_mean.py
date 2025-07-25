import argparse
import gffutils
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd
from collections import defaultdict

# Command-line arguments
parser = argparse.ArgumentParser(description="Calculate gene parameters from GTF and FASTA files.")
parser.add_argument("-g", "--gtf", required=True, help="Input GTF file.")
parser.add_argument("-f", "--fasta", required=True, help="Reference genome FASTA file.")
parser.add_argument("-o", "--output", required=True, help="Output file to save parameters.")

args = parser.parse_args()

# GTF database
db = gffutils.create_db(args.gtf, dbfn=':memory:', force=True, keep_order=True,disable_infer_genes=True, disable_infer_transcripts=True)

# Genome reference
genome = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

# Dictionary to collect per-gene data
gene_data = defaultdict(list)

# Loop over transcripts
for transcript in db.features_of_type('transcript'):
    gene_id = transcript.attributes.get('gene_id', [''])[0]
    gene_name = transcript.attributes.get('gene_name', [''])[0]
    transcript_id = transcript.id

    exons = list(db.children(transcript, featuretype='exon', order_by='start'))
    if not exons:
        continue

    # Sequence construction
    gene_seq = ""
    try:
        for exon in exons:
            gene_seq += genome[exon.seqid].seq[exon.start - 1:exon.end].upper()
    except KeyError:
        continue

    # Compute metrics
    gc = gc_fraction(gene_seq) * 100
    num_exons = len(exons)
    transcript_length = sum([exon.end - exon.start + 1 for exon in exons])
    gene_span_length = exons[-1].end - exons[0].start + 1

    utr5s = list(db.children(transcript, featuretype='5UTR'))
    utr3s = list(db.children(transcript, featuretype='3UTR'))

    len_5utr = sum([utr.end - utr.start + 1 for utr in utr5s])
    len_3utr = sum([utr.end - utr.start + 1 for utr in utr3s])

    # Save per-transcript stats under its gene
    gene_data[gene_id].append({
        'gene_name': gene_name,
        'gc_content': gc,
        'num_exons': num_exons,
        'transcript_length': transcript_length,
        'gene_span_length': gene_span_length,
        'mean_5utr_length': len_5utr / len(utr5s) if utr5s else 0,
        'mean_3utr_length': len_3utr / len(utr3s) if utr3s else 0
    })

# Now compute mean values per gene
aggregated = []

for gene_id, records in gene_data.items():
    n = len(records)
    gene_name = records[0]['gene_name']  # assume same for all

    aggregated.append({
        'gene_id': gene_id,
        'gene_name': gene_name,
        'mean_gc_content': sum(r['gc_content'] for r in records) / n,
        'mean_num_exons': sum(r['num_exons'] for r in records) / n,
        'mean_transcript_length': sum(r['transcript_length'] for r in records) / n,
        'mean_gene_length': sum(r['gene_span_length'] for r in records) / n,
        'mean_5utr_length': sum(r['mean_5utr_length'] for r in records) / n,
        'mean_3utr_length': sum(r['mean_3utr_length'] for r in records) / n,
    })

# Save to output
df = pd.DataFrame(aggregated)
df.to_csv(args.output, sep='\t', index=False)
