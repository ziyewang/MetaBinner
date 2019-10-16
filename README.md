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
conda env create -f environment.yml
conda activate metabinner_env
```

You may need to run these commands to make the files executable
```sh
chmod +x ~path_to_SolidBin/auxiliary/test_getmarker.pl
chmod +x ~path_to_SolidBin/auxiliary/FragGeneScan1.19/run_FragGeneScan.pl
chmod +x ~path_to_SolidBin/auxiliary/hmmer-3.1b1/bin/hmmsearch
chmod +x ~path_to_SolidBin/auxiliary/auxiliary/FragGeneScan1.19/FragGeneScan
```
