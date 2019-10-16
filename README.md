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
