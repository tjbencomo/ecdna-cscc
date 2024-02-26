"""
Extract genes in ecDNA amplicons from Kim 2020 NatGen ecDNA paper
41588_2020_678_MOESM2_ESM.xlsx is supplemental files from Kim 2020 Nature Genetics ecDNA paper
"""

import os
import pandas as pd
from collections import defaultdict
from collections import Counter

from interlap import InterLap


def parse_interval(x):
    xs = x.split(':')
    chrom = xs[0]
    start, end = xs[1].split('-')
    return (chrom, int(start), int(end))

def parse_interval_list(x):
    chroms = list(map(str, list(range(1, 22)))) + ['X', 'Y']
    xs = x.split(',')
    xs = [parse_interval(interval) for interval in xs]
    out = [interval for interval in xs if interval[0] in chroms]
    return out

def build_gene_db(genes):
    """
    genes is a list of tuples with format (gene, chrom, start, end)
    """
    inter = defaultdict(InterLap)
    for gene in genes:
        inter[gene[1]].add((gene[2], gene[3], gene[0]))
    return inter


natgen_fp = os.path.join('data', '41588_2020_678_MOESM2_ESM.xlsx')
natgen_df = pd.read_excel(natgen_fp)
# print(natgen_df.head())

ecdna_df = natgen_df[natgen_df['amplicon_classification'] == 'Circular']
# ecdna_df = natgen_df[natgen_df['amplicon_classification'] == 'Circular'].head(2)

genes_fp = os.path.join('data', 'grch37_gene_info.csv.gz')
genes_df = pd.read_csv(genes_fp)
genes_df = genes_df[genes_df['symbol'].notnull()]
print(genes_df.head())
# print(genes_df.dtypes)

genes = list(zip(genes_df['symbol'], genes_df['chr'], genes_df['start'], genes_df['end']))
gene_db = build_gene_db(genes) 


intervals = ecdna_df['amplicon_intervals'].apply(parse_interval_list)
interval_data = list(zip(ecdna_df['sample_barcode'].tolist(), intervals))
sampleGenes = {}
for x in interval_data:
    sampleID = x[0]
    genes = set()
    for interval in x[1]:
        overlap = list(gene_db[interval[0]].find((interval[1], interval[2])))
        for gene in overlap:
            genes.add(gene[2])
    if sampleID in sampleGenes:
        sampleGenes[sampleID] = sampleGenes[sampleID] | genes
    else:
        sampleGenes[sampleID] = genes


rows = [(sample_id, gene) for sample_id, genes in sampleGenes.items() for gene in genes]
out_df = pd.DataFrame(rows, columns=['sample_barcode', 'gene'])

output_fp = os.path.join('data', 'natgen_ecdna_genes.csv.gz')
out_df.to_csv(output_fp, index = False)

# masterList = []
# for sample in sampleGenes.keys():
    # print(f"{sample} has {len(sampleGenes[sample])} genes in ecDNA amplicons")
    # print(list(sampleGenes[sample])[0:5])
    # masterList = masterList + list(sampleGenes[sample])

# c = Counter(masterList)
# print(c.most_common()[0:2])
# print(len(c.most_common()))
# print(len(sampleGenes.keys()))

# print(masterList.count('MDM2'))
# print(masterList.count('MYC'))
# print(masterList.count('EGFR'))
# print(masterList.count('ERBB2'))
