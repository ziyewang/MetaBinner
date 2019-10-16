# MetaBinner: imporving large scale binning performance with ensemble K-means by considering multiple types of features.

## <a name="started"></a>Getting Started

### <a name="docker"></a>Conda

We recommend using conda to run Metabinner. Download [here](https://www.continuum.io/downloads)

After installing Anaconda (or miniconda), fisrt obtain Metabinner:

```sh
git clone https://github.com/ziyewang/MetaBinner
```
Then simply create a metabinner environment 

```sh
conda env create -f metabinner_env.yaml
conda activate metabinner_env
```

Set checkM data path as described in [here](https://github.com/Ecogenomics/CheckM/wiki/Installation)
```sh
checkm data setRoot <checkm_data_dir>
```

You may need to run these commands to make the files executable
```sh
chmod +x ~path_to_MetaBinner/auxiliary/test_getmarker.pl
chmod +x ~path_to_MetaBinner/auxiliary/FragGeneScan1.19/run_FragGeneScan.pl
chmod +x ~path_to_MetaBinner/auxiliary/hmmer-3.1b1/bin/hmmsearch
chmod +x ~path_to_MetaBinner/auxiliary/FragGeneScan1.19/FragGeneScan
```
### Coverage Profile
For the coverage profile we use minimap, since it can address both short read and long read samples.

You can put the files into your input directory, then slightly modify gen_cov.sh and run it.
minimap2, samtools and bedtools are need to be installed to run gen_cov.sh

Your input directory should look like this:

```
.
+-- assembly.fasta
+-- sr
|   +-- short_read_sample_1
|   +-- short_read_sample_2
+-- pb
|   +-- pacbio_sample_1
|   +-- pacbio_sample_2
|   +-- pacbio_sample_3
```

For conda environment, you should check whether perl is installed.

### Input files
You need to run like this to obtain the composition profiles and the coverage profiles (only short read samples; only long read samples)
```sh
scripts/run.sh input_path/contigs.fasta contigs_length_threshold kmer_length
```

"contigs_length_threshold" means only the contigs longer than the value are kept for binning. For example,
```sh
./run.sh input_path/contigs.fasta 1000 4
```
Note:The coverage_new.tsv file should be put into the "input_path" and if the sample numbers of the long read samples and the short read samples are different, you need to change the "split_coverage.py"

## <a name="usage"></a>Usage


> - Usage:         [--contig_file CONTIG_FILE]
                   [--coverage_profiles COVERAGE_PROFILES]
                   [--composition_profiles COMPOSITION_PROFILES]
                   [--output OUTPUT] [--log LOG] [--clusters CLUSTERS]
                   [--use_hmm] [--hmm_icm_path HMM_ICM_PATH]
                   [--hmm_file HMM_FILE]
                   [--pacbio_read_profiles PACBIO_READ_PROFILES]
                   [--binscore BINSCORE]

> - arguments

  	--contig_file CONTIG_FILE: 
              The contigs file.
	
  	--coverage_profiles COVERAGE_PROFILES: 
              The coverage_profiles, containing a table where each
              row correspond to a contig, and each column correspond
              to a sample. All values are separated with tabs. (the coverage profiles of short read samples)
  	--composition_profiles COMPOSITION_PROFILES: 
              The composition profiles, containing a table where
              each row correspond to a contig, and each column
              correspond to the kmer composition of particular kmer.
              All values are separated with comma.
	
  	--output OUTPUT:
              The output file, storing the binning result.
              
> - optional

    --log LOG_LOCATION:
              The log file, storing the binning process and parameters.
              
    --clusters CLUSTERS: 
              Specify the number of clusters. If not specified, the
              cluster number is estimated by single-copy genes.
              
    --pacbio_read_profiles PACBIO_READ_PROFILES:
              The coverage_profiles, containing a table where each
              row correspond to a contig, and each column correspond
              to a sample. All values are separated with tabs. (the coverage profiles of long read samples)
    
    --binscore BINSCORE:
              The value of Das_tool bin score threshold. Default is 0.3 in the tool.
              
    --use_hmm:
              Use hidden markov model score.
    
    --hmm_icm_path HMM_ICM_PATH:
             Metabinner/hmm_data/hmm/

```sh
OPENBLAS_NUM_THREADS=1 python Metabinner.py --contig_file input_path/contigs.fasta --coverage_profiles input_path/coverage_sr_new.tsv --composition_profiles input_path/kmer_4_f0.csv --output output_path/result.tsv --log output_path/result.log --pacbio_read_profiles input_path/coverage_pb_new.tsv --use_hmm --hmm_icm_path path_to_MetaBinner/hmm_data/hmm/
```

## <a name="preprocessing"></a>Contacts and bug reports
Please send bug reports or questions (such as the appropriate modes for your datasets) to
Ziye Wang: zwang17@fudan.edu.cn and Dr. Shanfeng Zhu: zhusf@fudan.edu.cn

## <a name="preprocessing"></a>References
         

[1] Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome Research, 25: 1043–1055.

[2] Graham ED, Heidelberg JF, Tully BJ. (2017) "BinSanity: unsupervised clustering of environmental microbial assemblies using coverage and affinity propagation." PeerJ 5:e3035

[3] Wang Z., Wang Z., et al.(2019) "SolidBin: Improving Metagenome Binning with Semi-supervised Normalized Cut." Bioinformatics. 2019 Apr 12. pii: btz253. 

[4] Christian M. K. Sieber, Alexander J. Probst., et al. (2018). "Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy". Nature Microbiology. https://doi.org/10.1038/s41564-018-0171-1.

[5] Kelley,D.,and Salzberg,S.(2010). "Clustering metagenomic sequences with interpolated markov models". BMC Bioinformatics 11:544. doi:10.1186/1471-2105-11-544

