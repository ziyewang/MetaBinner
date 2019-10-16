#eg: ./run.sh /mnt/data4/wzy/preprocess/strmg_megahit_assembly/test_input/strmgCAMI2_short_read_pooled_megahit_assembly.fasta 1000 4
# ./run.sh fasta_file  contig_length_threshold kmer_number


input_dir=$(dirname $1)

# prpocess.sh $1s
python Filter_tooshort.py $1 $2
python join.py ${input_dir} $2
python gen_kmer.py $1 $2 $3
python split_coverage.py ${input_dir} $2
