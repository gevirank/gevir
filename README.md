# GeVIR: Gene Variation Intolerance Rank
Gene Variation Intolerance Rank (GeVIR) is a gene level metric which can aid dominant and recessive Mendelian disease genes discovery.

### Project Description

With large scale population sequencing projects gathering pace there is a need for strategies that progress disease gene prioritisation. Metrics that provide information about a gene and its ability to tolerate protein altering variation can aid clinical interpretation of human genomes and advance disease gene discovery. Previous methods analysed total variant load in a gene, but not their distribution pattern within a gene. Utilising data from 138,632 exome/genome sequences, we developed Gene Variation Intolerance Rank (GeVIR), to produce a continuous gene level metric for 19,361 genes that is able to prioritise both dominant and recessive Mendelian disease genes, outperforming missense constraint metrics and comparable, but complementary, to loss-of-function constraint metrics. GeVIR is also able to prioritise short genes, for which loss-of-function constraint cannot be confidently estimated. The majority of the most intolerant genes identified have no defined phenotype and are candidates for severe dominant disorders. 

### Required Datasets
GeVIR analysis requires local version of gnomAD v2.0.1 database, specifically following collections: _exome_coverage_, _exome_variants_, _genome_variants_, _exons_, _genes_, _transcripts_. Instructions how to install gnomAD database can be found here:
https://github.com/macarthur-lab/gnomad_browser

Additionally, it requires following datasets:
- The GERP++ scores, can be obtained from
http://mendel.stanford.edu/SidowLab/downloads/gerp/hg19.GERP_scores.tar.gz
- The Ensembl coding sequences, can be obtained from
http://grch37.ensembl.org/biomart/martview/e19a3419814f357a70908603f15f7bae?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.sequences.ensembl_gene_id|hsapiens_gene_ensembl.default.sequences.ensembl_transcript_id|hsapiens_gene_ensembl.default.sequences.coding&FILTERS=&VISIBLEPANEL=resultspanel
- The Ensembl peptide sequences, can be obtained from
http://grch37.ensembl.org/biomart/martview/e19a3419814f357a70908603f15f7bae?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.sequences.ensembl_gene_id|hsapiens_gene_ensembl.default.sequences.ensembl_transcript_id|hsapiens_gene_ensembl.default.sequences.peptide&FILTERS=&VISIBLEPANEL=resultspanel
- The gnomAD gene constraint metrics, can be obtained from 
https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz
- The OMIM genemap2.txt file can be found, **_AFTER REGISTRATION_**, at
https://omim.org/downloads/
- The ClinVar datasets, can be obtained from
ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz and 
ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt
- The Constraned Coding Regions (CCRs) datasets, can be obtained from
https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz and
https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz
- The cell essential genes, can be obtained from MacArthur repository
https://github.com/macarthur-lab/gene_lists/blob/master/lists/NEGv1_subset_universe.tsv
- The cell non-essential genes, can be obtained from MacArthur repository
https://github.com/macarthur-lab/gene_lists/blob/master/lists/homozygous_lof_tolerant_twohit.tsv
- The Loss-of-Function tolerant genes (i.e. Null), can be obtained from MacArthur repository
https://github.com/macarthur-lab/gene_lists/blob/master/lists/CEGv2_subset_universe.tsv
- The human-mouse ortholog mapping file can be obtained from
http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
-The mouse heterozygous lethal genes can be obtained from http://mousemine.org/ by quiring the database with following constraints: subject.zygosity = 'ht' (heterozygous) AND ontologyTerm.name CONTAINS 'lethal'. Alternatively, import_data.py contains script to query mousemine database programatically.

**Note 1**: All datasets have to be placed into _./source_data_ directory and have names as defined at the top of _import_data.py_ script (can be changed there). ClinVar dataset requires reannotation with Ensembl VEP script, instructions how to install and use it can be found here: https://github.com/Ensembl/ensembl-vep. 

**Note 2**: _./tables_ folder contains gnomAD gene constraint scores, mouse heterozygous lethal, cell essential, cell non-essential and loss-of-function tolerant datasets which were used in original GeVIR study.

### Usage
Download all required datasets, place them in _./source_data_ folder.
To install required python modules run:
```
pip install -r requirements.txt
```
To import datasets into the local database run: 
```
python import_data.py
```
To calculate GeVIR scores run:
```
python gevir.py
```
To create gene sets used in the figures production run:
```
python gene_sets.py
```
To produce all the figures run:
```
python figures.py
```
To export gene scores from the database run:
```
python export_data.py
```

**Note**: Individual operations can be enabled/disabled via comments in the _main_ method in all scripts, check them before performing the analysis.

### Code Description

- **import_data.py** - imports data from various datasets into the local database

- **gevir.py** - computes Variant Intolerant Regions (VIRs) from gnomAD variant data and creates GeVIR gene scores

- **gene_sets.py** - combines gene scores (e.g. GeVIR, LOEUF, Missense z-score) into a single dataset, loads disease and essential gene lists from the database

- **export_data.py** - exports VIR and GeVIR gene scores data as csvs

- **figures.py** - performs gene scores evaluation, draws figures and reports statistics

- **gnomad_utils.py** - modified version of gnomAD browser code, obtained from:
https://github.com/macarthur-lab/gnomad_browser/blob/master/utils.py

- **common.py** - contains methods commonly used by other scripts

- **csv_reader.py** - custom csv reader, used to import data into local database
