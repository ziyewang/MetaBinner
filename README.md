# Metabinner
GitHub repository for the manuscript "Metabinner: a stand-alone ensemble binning method to recover individual genomes from complex microbial communities"

## <a name="started"></a>Getting Started

### <a name="docker"></a>Conda

We recommend using conda to run MetaBinner.

### <a name="docker"></a>Obtain codes and create an environment
After installing Anaconda (or miniconda), fisrt obtain MetaBinner:

```sh
git clone https://github.com/ziyewang/MetaBinner.git
```
Then simply create a metabinner environment 

```sh
cd metabinner_path
conda env create -f metabinner_env.yaml
conda activate metabinner_env
```

### <a name="docker"></a>Install checkM (python3 version) like this

(please make sure you have installed openssl)

```sh
cd CheckM-1.0.18
python setup.py install
```
Install checkM database:

CheckM relies on a number of precalculated data files which can be downloaded from https://data.ace.uq.edu.au/public/CheckM_databases/. (More details are available at https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm):

```sh
mkdir <checkm_data_dir>
cd <checkm_data_dir>
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xzf checkm_data_2015_01_16.tar.gz 
checkm data setRoot .
```

CheckM requires the following programs to be added to your system path:

HMMER (>=3.1b1)

prodigal (2.60 or >=2.6.1)
executable must be named prodigal and not prodigal.linux

pplacer (>=1.1)
;guppy, which is part of the pplacer package, must also be on your system path;
pplacer binaries can be found on the pplacer GitHub page

or you can install the dependencies as follows:
```sh
conda install -c bioconda prodigal
conda install -c bioconda hmmer 
conda install -c bioconda pplacer
```

## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate coverage profile and composition profile as input to our program.

There are several binning methods that can generate these two types of information (such as CONCOCT and MetaWRAP) and we provide one method to generate the input files as follows.
### Coverage Profile
The coverage profiles of the contigs for the results in the manuscript were obtained via MetaWRAP 1.2.1 script: ``binning.sh".

If users have obtained the coverage (depth) file generated for MaxBin (mb2_master_depth.txt) using MetaWRAP, they can run the following command to generate the input coverage file for MeatBinner:
```sh
cat mb2_master_depth.txt | cut -f -1,4- > ${out}/coverage_profile.tsv
```
or remove the contigs no longer than 1000bp like this:
```sh
cat mb2_master_depth.txt | awk '{if ($2>1000) print $0 }' | cut -f -1,4- > coverage_profile_f1k.tsv

```

If users would like to generate coverage from sequencing reads directly, then can run the script slightly modified from the "binning.sh" of MetaWRAP. The script support different types of sequencing reads, and the defalut type is "paired" ([readsX_1.fastq readsX_2.fastq ...]).
```
bash gen_coverage_file.sh -a contig_file \
-o output_dir_of_coveragefile \
path_to_sequencing_reads/*fastq

Options:

        -a STR          metagenomic assembly file
        -o STR          output directory
        -t INT          number of threads (default=1)
        -m INT          amount of RAM available (default=4)
        -l INT          minimum contig length to bin (default=1000bp).
        --single-end    non-paired reads mode (provide *.fastq files)
        --interleaved   the input read files contain interleaved paired-end reads

```

### Composition Profile

Composition profile is the vector representation of contigs and we use kmer (k=4 in the example) to generate this information. Users can keep the contigs longer than contig_length_threshold, such as 1000, for binning as follows:

```
python scripts/gen_kmer.py test_data/final.contigs_f1k.fa 1000 4 
```
Here we choose k=4. By default we usually keep contigs longer than 1000, you can specify a different number. The kmer_file will be generated in the /path/to/contig_file. 

And the user can run the following command to keep the contigs longer than 1000bp for binning.

```
python scripts/filter_tooshort.py test_data/final.contigs_f1k.fa 1000
```


## <a name="started"></a>An example to run MetaBinner:
Test data is available at https://drive.google.com/file/d/1a-IOOpklXQr_C4sgNxjsxGEkx-n-0aa4/view?usp=sharing
```sh

#path to MetaBinner
metabinner_path=/home/wzy/MetaBinner
##test data
#path to the input files for MetaBinner and the output dir:
contig_file=/home/wzy/MetaBinner/test_data/final_contigs_f1k.fa
output_dir=/home/wzy/MetaBinner/test_data/output
coverage_profiles=/home/wzy/MetaBinner/test_data/coverage_profile_f1k.tsv
kmer_profile=/home/wzy/MetaBinner/test_data/kmer_4_f1000.csv


bash run_metabinner.sh -a ${contig_file} -o ${output_dir} -d ${coverage_profiles} -k ${kmer_profile} -p ${metabinner_path}


#The file "metabinner_result.tsv" in the "${output_dir}/metabinner_res" is the final output.
```

## <a name="preprocessing"></a>Contacts and bug reports
Please send bug reports or questions to
Ziye Wang: zwang17@fudan.edu.cn and Dr. Shanfeng Zhu: zhusf@fudan.edu.cn

## <a name="preprocessing"></a>References

[1] Lu, Yang Young, et al. "COCACOLA: binning metagenomic contigs using sequence COmposition, read CoverAge, CO-alignment and paired-end read LinkAge." Bioinformatics 33.6 (2017): 791-798.

[2] https://github.com/dparks1134/UniteM.

[3] Parks, Donovan H., et al. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome research 25.7 (2015): 1043-1055.

[4] Christian M. K. Sieber, Alexander J. Probst., et al. (2018). "Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy". Nature Microbiology. https://doi.org/10.1038/s41564-018-0171-1.

[5] Uritskiy, Gherman V., Jocelyne DiRuggiero, and James Taylor. "MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis." Microbiome 6.1 (2018): 1-13.
