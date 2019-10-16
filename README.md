# MetaBinner：imporving large scale binning performance with K-means results of multiple features.

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
