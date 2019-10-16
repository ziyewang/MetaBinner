# input_dir: the directory containts all data we need
# assembly: fasta file consisting of contigs
# mapdir: a temp directory to save sam/bam files
# short_read_dir: contains short read samples if any
# pb_read_dir: contains pacbio read samples if any

# for non-docker env, please install samtools,minimap,bedtools and add them into path

# input_dir="/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_megahit_assembly/input"
# assembly="${input_dir}/strmgCAMI2_short_read_pooled_megahit_assembly.fasta"
# mapdir="${input_dir}/map"
# short_read_dir="/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_gold_assembly/input/sr"
# pb_read_dir="/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/strmg_gold_assembly/input/pb"

assembly=$1
input_dir=dirname ${assembly}
mapdir=${input_dir}/map
short_read_dir=${input_dir}/sr
pb_read_dir=${input_dir}/pb

if [ ! -d $mapdir ]; then
mkdir $mapdir
fi

samtools faidx $assembly;
awk -v OFS='\t' {'print $1,$2'} ${assembly}.fai > ${input_dir}/length.txt;

cnt=0;

for file in ${short_read_dir}/*;
do echo $file;
let cnt=cnt+1;
echo $cnt;
minimap2 -t 45 -ax sr $assembly $file > "${mapdir}/sr_${cnt}.sam";

done

for file in ${pb_read_dir}/*;
do echo $file;
let cnt=cnt+1;
echo $cnt;
minimap2 -t 45 -ax map-pb $assembly $file > "${mapdir}/pb_${cnt}.sam";

done

for file in ${mapdir}/*.sam
do
    stub=${file%.sam}
    echo $stub  
    samtools view -h -b -S $file > ${stub}.bam; samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam;samtools sort -m 1000000000 ${stub}.mapped.bam -o ${stub}.mapped.sorted.bam; bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g ${input_dir}/length.txt > ${stub}_cov.txt 
done
wait

for i in ${mapdir}/*_cov.txt
do
   echo $i
   stub=${i%_cov.txt}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > ${stub}_cov.csv
done


/mnt/data4/wzy/scripts/Collate.pl $mapdir > ${input_dir}/coverage.tsv

#下面是我自己的代码
perl -pe "s/,/\t/g;" ${input_dir}/coverage.tsv > ${input_dir}/coverage_new.tsv
