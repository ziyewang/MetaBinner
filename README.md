# MetaBinner
GitHub repository for the manuscript "MetaBinner: a high-performance and stand-alone ensemble binning method to recover individual genomes from complex microbial communities". We are glad that overall Metabinner achieve top performance on the CAMI II Challenge. Please refer to Meyer, F., et al[1] for the results of CAMI II Challenge.

MetaBinner consists of two modules: 1) “Component module” includes steps 1-4, developed for generating high-quality, diverse component binning results; and 2) “Ensemble module” includes step 5, developed for recovering individual genomes from the component binning results. MetaBinner is an ensemble binning method, but it does not need the outputs of other individual binners. Instead, MetaBinner generates multiple high-quality component binning results based on the proposed “partial seed” method, for further integration. Please see our manuscript for the details.

<p align="center">
<img src="https://github.com/ziyewang/MetaBinner/blob/master/figures/framework.png" width="550"/>
</p>

## <a name="started"></a>Getting Started

### <a name="docker"></a>Install MetaBinner via bioconda
```sh
conda create -n metabinner_env python=3.7.6
```
```sh
conda activate metabinner_env
conda install -c bioconda metabinner
```

### <a name="docker"></a>or Install MetaBinner via source code


Obtain codes and create an environment:
After installing Anaconda (or miniconda), fisrt obtain MetaBinner:

```sh
git clone https://github.com/ziyewang/MetaBinner.git
```
Then simply create a environment to run MetaBinner.

```sh
cd MetaBinner
conda env create -f metabinner_env.yaml
conda activate metabinner_env
```

## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate coverage profile and composition profile as input to our program.

There are several binning methods that can generate these two types of information (such as CONCOCT and MetaWRAP) and we provide one method to generate the input files as follows.
### Coverage Profile
The coverage profiles of the contigs for the results in the manuscript were obtained via MetaWRAP 1.2.1 script: ``binning.sh".

If users have obtained the coverage (depth) file generated for MaxBin (mb2_master_depth.txt) using MetaWRAP, they can run the following command to generate the input coverage file for MetaBinner:
```sh
cat mb2_master_depth.txt | cut -f -1,4- > coverage_profile.tsv
```
or remove the contigs no longer than 1000bp like this:
```sh
cat mb2_master_depth.txt | awk '{if ($2>1000) print $0 }' | cut -f -1,4- > coverage_profile_f1k.tsv

```

To generate coverage from sequencing reads directly, run the following script slightly modified from the "binning.sh" of MetaWRAP. The script support different types of sequencing reads, and the defalut type is "paired" ([readsX_1.fastq readsX_2.fastq ...]). If MetaBinner is installed via bioconda, users can obtain path_to_MetaBinner via running this command: $(dirname $(which run_metabinner.sh))

```sh
cd path_to_MetaBinner
cd scripts

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

Composition profile is the vector representation of contigs and we use kmer (k=4 in the example) to generate this information. To generate the composition profile and keep the contigs longer than contig_length_threshold, such as 1000, for binning, run the script as follows:

```
cd path_to_MetaBinner
cd scripts

python gen_kmer.py test_data/final.contigs_f1k.fa 1000 4 
```
Here we choose k=4. By default we usually keep contigs longer than 1000, you can specify a different number. The kmer_file will be generated in the /path/to/contig_file. 

And the users can run the following command to keep the contigs longer than 1000bp for binning.

```
cd path_to_MetaBinner
cd scripts

python Filter_tooshort.py test_data/final.contigs_f1k.fa 1000
```


## <a name="started"></a>An example to run MetaBinner:
Test data is available at https://drive.google.com/file/d/1a-IOOpklXQr_C4sgNxjsxGEkx-n-0aa4/view?usp=sharing
```sh

#path to MetaBinner
metabinner_path=/home/wzy/MetaBinner
Note: If users intall MetaBinner via bioconda, they can set metabinner_path as follows: metabinner_path=$(dirname $(which run_metabinner.sh))

##test data
#path to the input files for MetaBinner and the output dir:
contig_file=/home/wzy/MetaBinner/test_data/final_contigs_f1k.fa
output_dir=/home/wzy/MetaBinner/test_data/output
coverage_profiles=/home/wzy/MetaBinner/test_data/coverage_profile_f1k.tsv
kmer_profile=/home/wzy/MetaBinner/test_data/kmer_4_f1000.csv


bash run_metabinner.sh -a ${contig_file} -o ${output_dir} -d ${coverage_profiles} -k ${kmer_profile} -p ${metabinner_path}

Options:

        -a STR          metagenomic assembly file
        -o STR          output directory
        -d STR          coverage_profile.tsv; The coverage profiles, containing a table where each row correspond
                            to a contig, and each column correspond to a sample. All values are separated with tabs.
        -k STR          kmer_profile.csv; The composition profiles, containing a table where each row correspond to a contig,
                            and each column correspond to the kmer composition of particular kmer. All values are separated with comma.
        -p STR          path to MetaBinner; e.g. /home/wzy/MetaBinner
        -t INT          number of threads (default=1)
        -s STR          Dataset scale; eg. small,large,huge (default="huge"); Users can choose "huge" to run MetaBinner on huge datasets
                        with lower memory requirements.

#The file "metabinner_result.tsv" in the "${output_dir}/metabinner_res" is the final output.
```

## <a name="contact"></a>Contacts and bug reports
Please feel free to send bug reports or questions to
Ziye Wang: zwang17@fudan.edu.cn and Prof. Shanfeng Zhu: zhusf@fudan.edu.cn


## <a name="References"></a>References
[1] Meyer, F., Fritz, A., Deng, ZL. et al. Critical Assessment of Metagenome Interpretation: the second round of challenges. Nat Methods (2022). https://doi.org/10.1038/s41592-022-01431-4

[2] Lu, Yang Young, et al. "COCACOLA: binning metagenomic contigs using sequence COmposition, read CoverAge, CO-alignment and paired-end read LinkAge." Bioinformatics 33.6 (2017): 791-798.

[3] https://github.com/dparks1134/UniteM.

[4] Parks, Donovan H., et al. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome research 25.7 (2015): 1043-1055.

[5] Christian M. K. Sieber, Alexander J. Probst., et al. (2018). "Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy". Nature Microbiology. https://doi.org/10.1038/s41564-018-0171-1.

[6] Uritskiy, Gherman V., Jocelyne DiRuggiero, and James Taylor. "MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis." Microbiome 6.1 (2018): 1-13.

## <a name="Citation"></a>Citation

Wang, Z., Huang, P., You, R., Sun, F., & Zhu, S. (2021). MetaBinner: a high-performance and stand-alone ensemble binning method to recover individual genomes from complex microbial communities. bioRxiv. https://doi.org/10.1101/2021.07.25.453671
