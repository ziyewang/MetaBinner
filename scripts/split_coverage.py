
from __future__ import division
import os
import click
import pandas as pd


@click.command()
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('k', type=click.INT)
def main(input_dir, k):
    data_file = os.path.join(input_dir, 'coverage_new_f' + str(k) + '.tsv')
    pb_file = os.path.join(input_dir, 'coverage_new_f' + str(k) + '_pb.tsv')
    sr_file = os.path.join(input_dir, 'coverage_new_f' + str(k) + '_sr.tsv')

    data = pd.read_table(data_file, header=0, index_col=0, sep='\t')

    data_pb=data.copy()

    k = len(data.columns) // 2

    data_pb.drop(data_pb.columns[k:k*2], axis=1, inplace=True)
    data_pb.to_csv(pb_file, sep='\t', header=True)

    #data = pd.read_table(datafile, header=0, index_col=0, sep='\t')
    data.drop(data.columns[:k], axis=1, inplace=True)
    data.to_csv(sr_file, sep='\t', header=True)


if __name__ == "__main__":
    # datafile = r'/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_megahit_assembly/input/coverage_f1k.tsv'
    # pb_outfile='/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_megahit_assembly/input/coverage_f1k_pb.tsv'
    # sr_outfile='/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_megahit_assembly/input/coverage_f1k_sr.tsv'
    main()
