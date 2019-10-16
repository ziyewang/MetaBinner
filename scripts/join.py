import os
import sys
import click
import pandas as pd


@click.command()
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('k', type=click.INT)
def main(input_dir, k):
    # datafile = r'/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_megahit_assembly/input/coverage_new.tsv'
    # labelfile = r'/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_megahit_assembly/input/contig_length_filter500bp.txt'
    # outfile = '/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_megahit_assembly/input/coverage_new_f500bp.tsv'
    
    data_file = os.path.join(input_dir, 'coverage_new.tsv')
    label_file = os.path.join(input_dir, 'contig_length_filter' + str(k) + '.txt')
    output_file = os.path.join(input_dir, 'coverage_new_f' + str(k) + '.tsv')

    data = pd.read_table(data_file, header=0, index_col=0, sep='\t')
    contig_id = pd.read_table(label_file, header=None, index_col=0, sep='\t', names=['contig_id', 'a'])

    joined = data.join(contig_id, how="inner")
    #joined.drop(joined.columns[0], axis=1, inplace=True)#
    k = len(data.columns)
    joined.drop(joined.columns[k:k+1], axis=1, inplace=True)#
    joined.to_csv(output_file, sep='\t', header=True)


if __name__ == "__main__":
    main()
