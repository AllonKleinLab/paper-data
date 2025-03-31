For a reference genome, the [HT2019 assembly with KY21 gene models from the Ghost Database](http://ghost.zool.kyoto-u.ac.jp/download_ht.html) was used in combination with the mitochondrial genes from the [Ensembl KH genome assembly](https://useast.ensembl.org/Ciona_intestinalis/Info/Index). The available GFF3 file of the KY21 gene models was reformatted as a GTF file using [GffRead](http://dx.doi.org/10.12688/f1000research.23297.1) in order to be compatible with 10X Cell Ranger.

This folder contains files associated with this genome version:
- `HT_mt.fasta` - genome assembly fasta
- `HT_KY21_mt.genes.gtf` - gene annotations
- `KY21_with_Ens_mito_transcripts.fasta` - transcript sequences
- `KY21_with_Ens_mito_protein.fasta` - translated sequences

The first two files, `HT_mt.fasta` and `HT_KY21_mt.genes.gtf` were together used to geenerate a CellRanger genome, in the `cellranger_reference/` folder.
