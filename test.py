contig_file='/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim40_20_original/input/final_contigs_f1k.fa'
path = '/mnt/data3/wzy/ncbi_genomes_metawatt_process/filter_small_genome/'  # 需要改路径
files = os.listdir(path)
par = 40
output='/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim40_20_original/output/metabinner_pipeline/aaa_kmeans.csv'
# train_imm(path, files, par)  # 直接用训练好的就行
score_reads(path, files, contig_file, par,output)
likelihood, read_probs = get_read_probs(path, files,output)
df = pd.DataFrame(read_probs)
df = df.T
df[df < 1e-8] = 0
hmm_file = os.path.join(os.path.abspath(os.path.dirname('/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim40_20_original/output/metabinner_pipeline/aaa.kmeans.csv')), 'hmm_profile.tsv')
df.to_csv(hmm_file, sep='\t', header=True)

cov_file = '/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim40_20_original/input/Coverage_f1k.tsv'

nClass = sum((np.sum(X_hmm, axis=0) >= 100).astype(int))  # 100 can be reset.
km = KMeans(n_clusters=nClass, init='k-means++', n_jobs=-1, n_init=30, random_state=7)
km.fit(X_hmm)
idx = km.labels_
kmeans_length_weight_hmm_output = os.path.dirname(
    args.output) + '/kmeans_length_weight_hmm_result.csv'
save_result(idx, kmeans_length_weight_hmm_output, namelist)
kmeans_length_weight_hmm_output_dir = os.path.dirname(
    args.output) + '/kmeans_length_weight_hmm_result'
os.mkdir(kmeans_length_weight_hmm_output_dir)
gen_bins(contig_file, kmeans_length_weight_hmm_output, kmeans_length_weight_hmm_output_dir, "hmm_result")

save_result(idx, '/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim40_20_original/output/metabinner_pipeline/hmm_result.csv', namelist)

logger.info("Recluster the contigs from high_com_p_high_cont bins")
kmeans_partial_seed_initial_output_dir='/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim200_20_original/output/metabinner_pipeline_latest/kmeans_partial_seed_initial_ori_result'
high_com_p_high_cont_path = kmeans_partial_seed_initial_output_dir + "/High_completion_high_contamination"
recluster_hh_bins(high_com_p_high_cont_path, mapObj, X_t, length_weight, namelist)

# recluster other contigs
logger.info("Recluster other contigs.")
not_clustered_path = kmeans_partial_seed_initial_output_dir + "/others"
recluster_other_contigs(not_clustered_path, X_t, namelist, mapObj, length_weight)

convert(kmeans_partial_seed_initial_output_dir + '/good_bins', os.path.dirname(kmeans_partial_seed_initial_output_dir) + '/kmeans_seed_partial_with_postprocess.csv')

lengths = get_length(contig_file)
length_weight = []
for seq_id in namelist:
    length_weight.append(lengths[seq_id])

high_com_p_high_cont_path = kmeans_partial_seed_initial_output_dir + "/High_completion_high_contamination"
hh_files = os.listdir(high_com_p_high_cont_path)
for file_ in hh_files:
    if file_.endswith('.bin'):
        hh_contigs_id = []
        hh_contig_file = high_com_p_high_cont_path + '/' + file_
        for seq_record in SeqIO.parse(hh_contig_file, "fasta"):
            hh_contigs_id.append(seq_record.id)
        hh_contigs_id_number = [mapObj[x] for x in hh_contigs_id]
        X_t_hh_unclustered = X_t[hh_contigs_id_number]
        bin_number = gen_bestk(hh_contig_file, 0)

        hh_weight = []
        for i in range(len(hh_contigs_id_number)):
            hh_weight.append(length_weight[hh_contigs_id_number[i]])
        # seedurl不一定存在
        seedURL = hh_contig_file + ".seed"
        global seed_idx
        if os.path.exists(seedURL):
            seed_list = []
            with open(seedURL) as f:
                for line in f:
                    seed_list.append(line.rstrip('\n'))
            name_map = dict(zip(hh_contigs_id, range(len(hh_contigs_id))))
            seed_idx = [name_map[seed_name] for seed_name in seed_list]
            km = KMeans(n_clusters=bin_number, n_init=30, random_state=7)
        else:
            km = KMeans(n_clusters=bin_number, n_jobs=-1, n_init=30, random_state=7)

        km.fit(X_t_hh_unclustered, sample_weight=hh_weight)
        idx = km.labels_
        save_result_refine(idx, hh_contig_file + ".reclustered.csv",
                           namelist, hh_contigs_id_number)
        gen_bins(hh_contig_file, hh_contig_file + ".reclustered.csv",
                 os.path.dirname(high_com_p_high_cont_path) + '/good_bins', file_ + "_reclustered")



#############################
contig_file='/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim200_20_original/input/final_contigs_1k.fa'
path = '/mnt/data3/wzy/ncbi_genomes_metawatt_process/filter_small_genome/'  # 需要改路径
files = os.listdir(path)
par = 40
output='/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim200_20_original/output/metabinner_pipeline_latest/aaa_kmeans.csv'
# train_imm(path, files, par)  # 直接用训练好的就行
score_reads(path, files, contig_file, par,output)
likelihood, read_probs = get_read_probs(path, files,output)
df = pd.DataFrame(read_probs)
df = df.T
df[df < 1e-8] = 0
hmm_file = os.path.join(os.path.abspath(os.path.dirname('/mnt/data2/wzy/binning/simulation_sample/no_cutup_dataset/Sim200_20_original/output/metabinner_pipeline/aaa.kmeans.csv')), 'hmm_profile.tsv')
df.to_csv(hmm_file, sep='\t', header=True)


##########
#得到cami_length的字典
#

count = 0
for k in cami_length:
        if cami_length[k] > 2000:
            count += 1
fastx_file='/mnt/data4/wzy/CAMI2019/2nd_CAMI_Challenge_Marine_Dataset_Assembly/assembly/marmgCAMI2_short_read_pooled_megahit_assembly.fasta.gz'
cami_length_assembly = get_length(fastx_file)

fastx_file='/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/Strain_Madness_Dataset_Assembly/assembly/strmgCAMI2_short_read_pooled_gold_standard_assembly.fasta.gz'
strain_gold_length=get_length(fastx_file)

count = 0
for k in strain_gold_length:
        if strain_gold_length[k] > 500:
            count += 1
print(count)

fastx_file='/mnt/data4/wzy/CAMI2019/Strain_Madness_Dataset/Strain_Madness_Dataset_Assembly/assembly/strmgCAMI2_short_read_pooled_megahit_assembly.fasta.gz'
strain_assembly_length=get_length(fastx_file)

count = 0
for k in strain_assembly_length:
        if strain_assembly_length[k] > 2000:
            count += 1
print(count)

logger.info("start calculate contig length")
lengths = get_length(contig_file)
length_weight = []
for seq_id in namelist:
    length_weight.append(lengths[seq_id])
logger.info("Recluster the contigs from high_com_p_high_cont bins")
high_com_p_high_cont_path = "/mnt/data2/wzy/CAMI_2018/test/short_read/output_for_review/metabinner_pipeline/kmeans_partial_seed_length_weight_result/hh_new"
recluster_hh_bins(high_com_p_high_cont_path, mapObj, X_t, length_weight, namelist)

def recluster_hh_bins(high_com_p_high_cont_path,mapObj,X_t,length_weight,namelist):
    hh_files = os.listdir(high_com_p_high_cont_path)
    for file_ in hh_files:
        if file_.endswith('.bin'):
            hh_contig_file = high_com_p_high_cont_path + '/' + file_
            gen_bins(hh_contig_file, hh_contig_file+".reclustered.csv",
                     os.path.dirname(high_com_p_high_cont_path) + '/good_bins',file_+ "_reclustered")

############
#用kmeans_partial_seed_ori来做后处理
checkm_file='/mnt/data2/wzy/CAMI_2018/test/short_read/output_for_review/metabinner_pipeline/kmeans_partial_seed_initial_ori_result/checkm_out/checkm_result.txt'
checkm_out_dir='/mnt/data2/wzy/CAMI_2018/test/short_read/output_for_review/metabinner_pipeline/kmeans_partial_seed_initial_ori_result/checkm_out'
suffix_str = '.bin'
checkm_analysis(checkm_file, suffix_str, checkm_out_dir)

kmeans_partial_seed_length_weight_output_dir='/mnt/data2/wzy/CAMI_2018/test/short_read/output_for_review/metabinner_pipeline/kmeans_partial_seed_initial_ori_result'
logger.info("Recluster the contigs from high_com_p_high_cont bins")
high_com_p_high_cont_path = kmeans_partial_seed_length_weight_output_dir + "/High_completion_high_contamination"
recluster_hh_bins(high_com_p_high_cont_path, mapObj, X_t, length_weight, namelist)

# recluster other contigs
logger.info("Recluster other contigs.")
not_clustered_path = kmeans_partial_seed_length_weight_output_dir + "/others"
recluster_other_contigs(not_clustered_path, X_t, namelist, mapObj, length_weight)

convert(kmeans_partial_seed_length_weight_output_dir + '/good_bins', kmeans_partial_seed_length_weight_output_dir + '/kmeans_seed_partial_ori_with_postprocess_change_hh.csv')

def recluster_hh_bins(high_com_p_high_cont_path,mapObj,X_t,length_weight,namelist):
    hh_files = os.listdir(high_com_p_high_cont_path)
    for file_ in hh_files:
        if file_.endswith('.bin'):
            hh_contigs_id = []
            hh_contig_file = high_com_p_high_cont_path + '/' + file_
            for seq_record in SeqIO.parse(hh_contig_file, "fasta"):
                hh_contigs_id.append(seq_record.id)
            hh_contigs_id_number = [mapObj[x] for x in hh_contigs_id]
            X_t_hh_unclustered = X_t[hh_contigs_id_number]
            bin_number = gen_bestk(hh_contig_file,0)

            hh_weight = []
            for i in range(len(hh_contigs_id_number)):
                hh_weight.append(length_weight[hh_contigs_id_number[i]])
             #seedurl不一定存在
            seedURL = hh_contig_file + ".seed"
            global seed_idx
            if os.path.exists(seedURL):
                seed_list = []
                with open(seedURL) as f:
                    for line in f:
                        seed_list.append(line.rstrip('\n'))
                name_map = dict(zip(hh_contigs_id, range(len(hh_contigs_id))))
                seed_idx = [name_map[seed_name] for seed_name in seed_list]
                km = KMeans(n_clusters=bin_number, n_init=30, random_state=7)
            else:
                km = KMeans(n_clusters=2, n_jobs=-1,n_init=30, random_state=7)#no seed,k=2

            km.fit(X_t_hh_unclustered, sample_weight=hh_weight)
            idx = km.labels_
            save_result_refine(idx, hh_contig_file + ".reclustered.csv",
                               namelist, hh_contigs_id_number)
            gen_bins(hh_contig_file, hh_contig_file+".reclustered.csv",
                     os.path.dirname(high_com_p_high_cont_path) + '/good_bins',file_+ "_reclustered")


