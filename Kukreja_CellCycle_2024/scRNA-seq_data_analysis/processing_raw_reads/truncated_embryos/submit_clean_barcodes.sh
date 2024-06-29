
# Mismatch threshold
max_hamming_dist=1

# Output directory
outdir=output_tags
mkdir -p $outdir

# Be sure to have expected_barcodes.txt in folder

# Run get_cell_hashing_barcodes_flex.py
f=TAG_1

echo $f
cat $f/filtered_parts/*sorted.fastq.gz | gunzip -c | \
    python get_cell_hashing_barcodes_flex.py expected_barcodes.txt $outdir/$f.counts.csv $outdir/$f.reads.csv $outdir/$f.pickle ${max_hamming_dist}
echo done


f=TAG_2

echo $f
cat $f/filtered_parts/*sorted.fastq.gz | gunzip -c | \
    python get_cell_hashing_barcodes_flex.py expected_barcodes.txt $outdir/$f.counts.csv $outdir/$f.reads.csv $outdir/$f.pickle ${max_hamming_dist}
echo done


f=TAG_3

echo $f
cat $f/filtered_parts/*sorted.fastq.gz | gunzip -c | \
    python get_cell_hashing_barcodes_flex.py expected_barcodes.txt $outdir/$f.counts.csv $outdir/$f.reads.csv $outdir/$f.pickle ${max_hamming_dist}
echo done


f=TAG_4

echo $f
cat $f/filtered_parts/*sorted.fastq.gz | gunzip -c | \
    python get_cell_hashing_barcodes_flex.py expected_barcodes.txt $outdir/$f.counts.csv $outdir/$f.reads.csv $outdir/$f.pickle ${max_hamming_dist}
echo done
