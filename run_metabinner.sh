#!/usr/bin/env bash

##############################################################################################################################################################
# This script is meant to be run MetaBinner after obtaining the coverage file and the kmer profile.
# Author of pipeline: Ziye Wang.
# For questions, bugs, and suggestions, contact me at zwang17@fudan.edu.cn
##############################################################################################################################################################

help_message () {
	echo "Options:"
	echo ""
	echo "	-a STR          metagenomic assembly file"
	echo "	-o STR          output directory"
	echo "  -d STR          coverage_profile.tsv; The coverage profiles, containing a table where each row correspond
	                        to a contig, and each column correspond to a sample. All values are separated with tabs."
	echo "  -k STR          kmer_profile.csv; The composition profiles, containing a table where each row correspond to a contig,
	                        and each column correspond to the kmer composition of particular kmer. All values are separated with comma."
	echo "  -p STR          path to metabinner; e.g. /home/wzy/MetaBinner"
	echo "	-t INT          number of threads (default=10)"
	echo "";}

num_threads=10

while getopts a:o:d:k:p:t: OPT; do
 case ${OPT} in
  a) contig_file=${OPTARG}
    ;;
  o) output_dir=${OPTARG}
    ;;
  d) coverage_profiles=${OPTARG}
    ;;
  k) kmer_profile=${OPTARG}
    ;;
  p) metabinner_path=${OPTARG}
    ;;
  t) num_threads=${OPTARG}
    ;;
  \?)
#    printf "[Usage] `date '+%F %T'` -i <INPUT_FILE> -o <OUTPUT_DIR> -o <P
#RODUCT_CODE> -s <SOFTWARE_VERSION> -t <TYPE>\n" >&2
    exit 1
 esac
done

# check parameter
if [ -z "${contig_file}" -o -z "${output_dir}" -o -z "${coverage_profiles}" -o -z "${kmer_profile}" -o -z "${metabinner_path}" ]; then
  help_message
  exit 1
fi

contig_file_name=`basename ${contig_file}`

########################################################################################################
########################            generate component binning results          ########################
########################################################################################################


mkdir ${output_dir}
mkdir ${output_dir}/metabinner_res

cd ${metabinner_path}/scripts/
time ./component_binning.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_profile} \
--output ${output_dir}/metabinner_res/result.tsv \
--log ${output_dir}/metabinner_res/result.log \
--threads ${num_threads}

if [[ $? -ne 0 ]] ; then echo "Something went wrong with generating component binning results. Exiting.";exit 1; fi

#unitem/checkm_ms/
########################################################################################################
#########   generate marker gene information using one component binning result and checkm       #######
#########   note: Some of the codes in this part are modified from UniteM.                       #######
#########   https://github.com/dparks1134/UniteM                                                 #######
########################################################################################################


mkdir ${output_dir}/metabinner_res/unitem_profile
cp -r ${output_dir}/metabinner_res/intermediate_result/kmeans_length_weight_X_t_logtrans_result.tsv ${output_dir}/metabinner_res/unitem_profile
res_path=${output_dir}/metabinner_res/unitem_profile/kmeans_length_weight_X_t_logtrans_result.tsv

${metabinner_path}/scripts/gen_bins_from_tsv.py \
-f ${contig_file} \
-r ${res_path} \
-o ${res_path}_bins

echo -e "X_t_logtrans_ori\t${res_path}_bins" > ${output_dir}/metabinner_res/unitem_profile/bins_dir.tsv

${metabinner_path}/scripts/unitem_profile.py --threads ${num_threads} --bin_dirs ${output_dir}/metabinner_res/unitem_profile/bins_dir.tsv \
--output_dir ${output_dir}/metabinner_res/unitem_profile

if [[ $? -ne 0 ]] ; then echo "Something went wrong with running unitem_profile.py. Exiting.";exit 1; fi

########################################################################################################
# Step 4: Split bins with high contamination and completeness according to the
# single-copy marker genes.
########################################################################################################

bac_mg_table=${output_dir}/metabinner_res/unitem_profile/binning_methods/X_t_logtrans_ori/checkm_bac/marker_gene_table.tsv
ar_mg_table=${output_dir}/metabinner_res/unitem_profile/binning_methods/X_t_logtrans_ori/checkm_ar/marker_gene_table.tsv

if [ ! -s ${bac_mg_table} ] ; then echo "Something went wrong with running unitem_profile.py. Please check CheckM installation. Exiting.";exit 1; fi
if [ ! -s ${ar_mg_table} ] ; then echo "Something went wrong with running unitem_profile.py. Please check CheckM installation. Exiting.";exit 1; fi

logtrans_namelist_file=${metabinner_path}/scripts/logtrans_namelist.tsv
components_path=${output_dir}/metabinner_res/intermediate_result

cd ${components_path}
ls *tsv > ${components_path}.res_namelist.tsv

cat ${components_path}.res_namelist.tsv | while read LINE
do
echo $LINE;
${metabinner_path}/scripts/gen_bins_from_tsv.py \
-f ${contig_file} \
-r ${components_path}/${LINE} \
-o ${components_path}/${LINE}_bins
done

cd ${metabinner_path} #for hmmpath

cat ${logtrans_namelist_file} | while read LINE
do
echo $LINE;
ori_path=${components_path}/${LINE}_bins

${metabinner_path}/scripts/split_hhbins.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_profile} \
--bac_mg_table ${bac_mg_table} \
--ar_mg_table ${ar_mg_table} \
--ori_result_path ${ori_path} \
--log ${output_dir}/metabinner_res/post_result.log \
--threads ${num_threads}
done

if [[ $? -ne 0 ]] ; then echo "Something went wrong with running split_hhbins.py. Exiting.";exit 1; fi


###postprocess again

cat ${logtrans_namelist_file} | while read LINE
do
echo $LINE;
ori_path=${components_path}/${LINE}_bins_post_process_mincomp_70_mincont_50_bins

cd ${metabinner_path}

${metabinner_path}/scripts/split_hhbins.py \
--contig_file ${contig_file} \
--coverage_profiles ${coverage_profiles} \
--composition_profiles ${kmer_profile} \
--bac_mg_table ${bac_mg_table} \
--ar_mg_table ${ar_mg_table} \
--ori_result_path ${ori_path} \
--log ${output_dir}/metabinner_res/post_result.log \
--threads ${num_threads}
done

if [[ $? -ne 0 ]] ; then echo "Something went wrong with running split_hhbins.py. Exiting.";exit 1; fi



########################################################################################################
# obtain bins_dir.tsv
########################################################################################################

rm -rf ${components_path}/2postprocess_X_t_logtrans_bins_dir.tsv
rm -rf ${components_path}/2postprocess_X_cov_logtrans_bins_dir.tsv
rm -rf ${components_path}/2postprocess_X_com_logtrans_bins_dir.tsv

echo -e "kmeans_length_weight_X_t_logtrans_result.tsv\t${components_path}/kmeans_length_weight_X_t_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_t_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_1quarter_X_t_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_1quarter_X_t_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_t_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_2quarter_X_t_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_2quarter_X_t_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_t_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_t_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_t_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_t_logtrans_bins_dir.tsv

echo -e "kmeans_length_weight_X_cov_logtrans_result.tsv\t${components_path}/kmeans_length_weight_X_cov_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_cov_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_1quarter_X_cov_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_1quarter_X_cov_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_cov_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_2quarter_X_cov_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_2quarter_X_cov_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_cov_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_cov_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_cov_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_cov_logtrans_bins_dir.tsv

echo -e "kmeans_length_weight_X_com_logtrans_result.tsv\t${components_path}/kmeans_length_weight_X_com_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_com_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_1quarter_X_com_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_1quarter_X_com_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_com_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_2quarter_X_com_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_2quarter_X_com_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_com_logtrans_bins_dir.tsv
echo -e "partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_com_logtrans_result.tsv\t${components_path}/partial_seed_kmeans_bacar_marker_seed_length_weight_3quarter_X_com_logtrans_result.tsv_bins_post_process_mincomp_70_mincont_50_bins_post_process_mincomp_70_mincont_50_bins" >>${components_path}/2postprocess_X_com_logtrans_bins_dir.tsv


mkdir ${output_dir}/metabinner_res/ensemble_res

mkdir ${output_dir}/metabinner_res/ensemble_res/X_t_logtrans_2postprocess
mkdir ${output_dir}/metabinner_res/ensemble_res/X_cov_logtrans_2postprocess
mkdir ${output_dir}/metabinner_res/ensemble_res/X_com_logtrans_2postprocess

########################################################################################################
# generate first-stage ensemble result
########################################################################################################

path=${output_dir}/metabinner_res/ensemble_res
Method_name=greedy_cont_weight_3_mincomp_50.0_maxcont_15.0_bins

${metabinner_path}/scripts/ensemble.py ${bac_mg_table} \
${ar_mg_table} \
${components_path}/2postprocess_X_t_logtrans_bins_dir.tsv \
${output_dir}/metabinner_res/ensemble_res/X_t_logtrans_2postprocess &
${metabinner_path}/scripts/ensemble.py ${bac_mg_table} \
${ar_mg_table} \
${components_path}/2postprocess_X_cov_logtrans_bins_dir.tsv \
${output_dir}/metabinner_res/ensemble_res/X_cov_logtrans_2postprocess &
${metabinner_path}/scripts/ensemble.py ${bac_mg_table} \
${ar_mg_table} \
${components_path}/2postprocess_X_com_logtrans_bins_dir.tsv \
${output_dir}/metabinner_res/ensemble_res/X_com_logtrans_2postprocess &

wait

if [[ $? -ne 0 ]] ; then echo "Something went wrong with running first-stage ensemble (ensemble.py). Exiting.";exit 1; fi

########################################################################################################
# generate second-stage ensemble result
########################################################################################################

binsA=${path}/X_t_logtrans_2postprocess/${Method_name}
binsB=${path}/X_cov_logtrans_2postprocess/${Method_name}
binsC=${path}/X_com_logtrans_2postprocess/${Method_name}

mkdir ${path}/${Method_name}
cd ${path}/${Method_name}

mkdir ensemble_3logtrans
cd ensemble_3logtrans

cp -r ${binsA} X_t_logtrans_${Method_name}
cp -r ${binsB} X_cov_logtrans_${Method_name}
cp -r ${binsC} X_com_logtrans_${Method_name}

rm -rf bins_dir.tsv
rm -rf addrefined3comps_bins_dir.tsv
rm -rf addrefined2and3comps_bins_dir.tsv

echo -e "X_t_logtrans\t${binsA}" >> bins_dir.tsv
echo -e "X_cov_logtrans\t${binsB}" >> bins_dir.tsv
echo -e "X_com_logtrans\t${binsC}" >> bins_dir.tsv

binsA=X_t_logtrans_${Method_name}
binsB=X_cov_logtrans_${Method_name}
binsC=X_com_logtrans_${Method_name}


${metabinner_path}/scripts/binning_refiner.py -1 ${binsA} -2 ${binsB} -3 ${binsC} -o Refined_ABC > Refined_ABC.out &
${metabinner_path}/scripts/binning_refiner.py -1 ${binsA} -2 ${binsB} -o Refined_AB > Refined_AB.out &
${metabinner_path}/scripts/binning_refiner.py -1 ${binsA} -2 ${binsC} -o Refined_AC > Refined_AC.out &
${metabinner_path}/scripts/binning_refiner.py -1 ${binsB} -2 ${binsC} -o Refined_BC > Refined_BC.out &

wait

if [[ $? -ne 0 ]] ; then echo "Something went wrong with running binning_refiner.py. Exiting.";exit 1; fi


mv Refined_ABC/Refined Refined_ABC/Refined_ABC
mv Refined_AB/Refined Refined_AB/Refined_AB
mv Refined_AC/Refined Refined_AC/Refined_AC
mv Refined_BC/Refined Refined_BC/Refined_BC


cp -r bins_dir.tsv addrefined3comps_bins_dir.tsv
echo -e "Refined_ABC\t${path}/${Method_name}/ensemble_3logtrans/Refined_ABC/Refined_ABC" >> addrefined3comps_bins_dir.tsv

cp -r addrefined3comps_bins_dir.tsv addrefined2and3comps_bins_dir.tsv
echo -e "Refined_AB\t${path}/${Method_name}/ensemble_3logtrans/Refined_AB/Refined_AB" >> addrefined2and3comps_bins_dir.tsv
echo -e "Refined_AC\t${path}/${Method_name}/ensemble_3logtrans/Refined_AC/Refined_AC" >> addrefined2and3comps_bins_dir.tsv
echo -e "Refined_BC\t${path}/${Method_name}/ensemble_3logtrans/Refined_BC/Refined_BC" >> addrefined2and3comps_bins_dir.tsv

mkdir addrefined2and3comps

${metabinner_path}/scripts/ensemble.py ${bac_mg_table} \
${ar_mg_table} \
${output_dir}/metabinner_res/ensemble_res/${Method_name}/ensemble_3logtrans/addrefined2and3comps_bins_dir.tsv \
${output_dir}/metabinner_res/ensemble_res/${Method_name}/ensemble_3logtrans/addrefined2and3comps

if [[ $? -ne 0 ]] ; then echo "Something went wrong with running ensemble.py. Exiting.";exit 1; fi

cp -r ${output_dir}/metabinner_res/ensemble_res/${Method_name}/ensemble_3logtrans/addrefined2and3comps/${Method_name}_res.tsv ${output_dir}/metabinner_res/metabinner_result.tsv

echo "Binning Finished!"


