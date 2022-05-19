import os
import sys
import argparse
import re

def parse_args():
	AP = argparse.ArgumentParser("Collect 2019nCoV info")
	AP.add_argument('-report_dir',help='report dir [/results/analysis/output/Home/Auto_xxx]',dest='report_dir')
	AP.add_argument('-outdir',help='output dir [/results/analysis/output/Home/Auto_xxx/plugin_out/ReportSummary_2019nCoV.xxx/]',dest='outdir')
	return AP.parse_args()


def get_bams(report_dir):
	all_bams = glob.glob("%s/*_rawlib.bam" % (report_dir)) # IonXpress_053_rawlib.bam
	bams = []
	for bam in all_bams:
		name = os.path.basename(bam)
		bams.append(bam)

	return(sort(bams))

def get_basic_info(infile):
	'''
	get below info from 'ion_params_00.json' file
		* seq time
		* expName
		* chip type
	'''
	with open(infile,'r') as json_file:
		json_str = json_load(json_file)
		seq_date = json_str.get('date','NA')
		expName = json_str.get('expName','NA')
		chipType = json_str['exp_json'].get('chipType','NA')

	return([seq_date,expName,chipType])

def barcode_to_sampleName(infile): # ion_params_00.json
	bc2name = {}
	with open(infile,'r') as json_file:
		json_str = json_load(json_file)
		for sample in json_str['experimentAnalysisSettings']['barcodedSamples'].keys()
			name = sample
			bc   = json_str['experimentAnalysisSettings']['barcodedSamples'][name]['barcodes']
			bc2name[bc] = name

	return(bc2name)

def get_bc_reads_num(infile): # IonCode_0301_rawlib.ionstats_alignment.json
	with open(infile,'r') as json_file:
		json_str = json_load(json_file)
		reads_num = json_str['full']['num_reads']
		mean_len  = json_str['full']['mean_read_length']
		#max_len   = json_str['full']['max_read_length']
		#total_base_num = json_str['full']['num_bases']
		#q20_base_num = ''

	return([reads_num,mean_len])

def get_first_run_plugin_result(report_dir,plugin_name):
	'''
	SARS_CoV_2_coverageAnalysis
	SARS_CoV_2_variantCaller
	generateConsensus
	
	'''
	start_json_list = glob.glob("%s/%s_out.*/startplugin.json" % (report_dir,plugin_name))
	if len(start_json_list) == 0:
		return([])
	else:
		xxx_int = []
		for file in start_json_list:
			plugin_full_path = os.path.dirname(file)
			plugin_name = os.path.basename(plugin_full_path)
			xxx = int(plugin_name.split('.')[1])
			xxx_int.append(xxx)

		first_xxx = xxx_int.sort()[0]
		first_name = "%s_.%s" % (plugin_name,first_xxx)
		return(first_name) # SARS_CoV_2_coverageAnalysis_out.xxx

def get_cov_stat(infile,barcode): # plugin_out/SARS_CoV_2_coverageAnalysis_out.xxx/*.bc_summary.xls
	'''
	Barcode ID      Sample Name     Mapped Reads    Filtered Reads  Target Reads    Mean Depth      Uniformity
	IonXpress_001   001 FluB        6578122 0.32%   98.04%  28274   97.10%
	IonXpress_003   002 FluB        635     0.00%   100.00% 4.5     49.97%
	IonXpress_004   003 FluB        53592   0.15%   99.56%  332.6   73.31%

	'''
	with open(infile,'r') as cov_summary:
		for line in cov_summary.readline():
			if line.startwith('Barcode ID'):
				continue
			else: 
				vals = line.split('\t')
				if vals[0] == barcode: # this sample line
					mapped_reads_n = vals[2]
					on_target = vals[4]
					mean_depth = vals[-2]
					uni = vals[-1]

	return([mapped_reads_n,on_target,mean_depth,uni])

def reads_per_pool(infile):
	'''
	/plugin_out/SARS_CoV_2_coverageAnalysis_out.1723/IonXpress_001/*.amplicon.cov.xls

	contig_id       contig_srt      contig_end      region_id       attributes      gc_count        overlaps        fwd_e2e rev_e2e total_reads     fwd_reads       rev_reads       cov20x  cov100x cov500x
	2019-nCoV       3604    3824    r1_1.4.1477602  GENE_ID=r1;Pool=2       79      1494    0       0       0       0       0       209
     209     172
	2019-nCoV       4020    4239    r1_1.5.1289446  GENE_ID=r1;Pool=1       83      1859    0       0       0       0       0       187
     187     90
	2019-nCoV       6805    7017    r1_1.8.592180   GENE_ID=r1;Pool=2       64      918     0       0       0       0       0       189
     175     27
	'''
	p1_num = []
	p2_num = []
	
	with open(infile,'r') as amplicon_info:
		amplicon_info.readline() # skip header line
		for line in amplicon_info.readline():
			vals = line.split('\t')
			pool_info = vals[4] # GENE_ID=LDLRAP1;POOL=1;CNV_HS=0;CNV_ID=LDLRAP1
			total_reads = vals[-6]
			for v in pool_info.split(';'):
				if re.match('POOL',v):
					if v == 'POOL=1':
						p1_num.append(total_reads)
					if v == 'POOL=2':
						p2_num.append(total_reads)

	avg_p1 = float(mean(p1_num),2)
	avg_p2 = float(mean(p2_num),2)

	return([avg_p1,avg_p2])

def get_tvc_info(infile,barcode):
	'''
	variantCaller
	SARS_CoV_2_coverageAnalysis
	
	/plugin_out/SARS_CoV_2_variantCaller_out.1733/results.json
	'''
	with open(infile,'r') as json_file:
		json_str = json_load(json_file)
		var_num = json_str['barcodes'][barcode]['variants']['variants']
		return(var_num)

def get_cons_info(infile,barcode):
	'''
	/plugin_out/generateConsensus_out.1782/results.json
	
	"Launch_Mode": "Manual",
  	"Major_Allele_Only": "Yes",
  	"Maximum_Percent_N": "1",
  	"Minimum_Read_Depth": "20",
  	"Minimum_Variant_Frequency": "0.5",
  	"Minimum_Variant_Frequency_HPInDels": "0.6",
  	"barcodes": {
  		"IonXpress_001": {
  			"Chromosome": "2019-nCoV",
  			"Contig Length": "29753",
  			"Het indels": "5",
  			"Het snps": "15",
  			"Homo indels": "1",
  			"Homo snps": "30",
  			"Other": "4",
  			"Percent N": "0.2185",
  			"Variants": "55"
		}
	}

	return value:
		* 组装N比例
		* 一致性序列变异位点个数
		* 一致性序列杂合SNP个数 
	'''
	with open(infile,'r') as json_file:
		json_str = json_load(json_file)
		pct_N = json_str['barcodes'][barcode].get('Percent N','NA')
		var_num = json_str['barcodes'][barcode].get('Variants','NA')
		het_snp = json_str['barcodes'][barcode].get('Het snps','NA')
		#het_indel = json_str['barcodes'][barcode]['Het indels']
		return([pct_N,var_num,het_snp])

def get_pangolin_info(infile,barcode):
	'''
	plugin_out/SARS_CoV_2_lineageID_out.1553/results.json
	'''
	with open(infile,'r') as json_file:
		json_str = json_load(json_file)
		Lineage = json_str['barcodes'][barcode].get('Lineage','NA')
		return(Lineage)

def main():
	args = parse_args()
	report_name = os.path.basename(args.report_dir)
	outfile = "%s/%s.summary.xls" % (args.outdir,report_name)
	of = open(outfile,'w')
	header = "\t".join("建库日期","测序日期","expName","报告名称","芯片类型","芯片总数据量","Barcode","样本名","Pangolin分型","Nextclade分型","样本数据量","均一性","组装N比例","TVC变异位点个数","一致性序列变异位点个数","一致性序列杂合SNP个数","Q20碱基百分比","Reads平均长度","比对Reads数","Ontarget率","平均测序深度","Pool1-Mean Reads per Amplicon","Pool2-Mean Reads per Amplicon","是否提交","提交日期","备注")
	of.write(header+'\n')
	
	ion_params_00_json = os.path.join(args.report_dir,'ion_params_00.json')
	# check ion_params_00.json exists
	if os.path.exists(ion_params_00_json):
		print(ion_params_00_json)
		continue
	else:
		sys.exit('[Error: can not find ion_params_00.json file, will exit]')

	# get libary date
	chef_date = 'NA'

	# get sequencing date / expName / chip type
	basic_info = get_basic_info(ion_params_00_json)
	seq_date = basic_info[0]
	expName  = basic_info[1]
	chipType = basic_info[2]

	# get report name
	report_name = os.path.basename(args.report_dir)

	# get total reads number
	total_reads = 'NA'

	############################## for each barcode ##############################
	#bams = get_bams(args.report_dir) # [IonXpress_003_rawlib.bam,...,]
	all_barcodes = barcode_to_sampleName(ion_params_00_json)
	for bc in all_barcodes:
		# IonCode_0303
		sample_name = all_barcodes[bc] # maybe None
		
		# pangolin result
		pangolin_dir = get_first_run_plugin_result(args.report_dir,'SARS_CoV_2_lineageID')
		if len(pangolin_dir) == 0:
			# no SARS_CoV_2_lineageID plugin run
			pangolin_result = 'NA'
		else:
			p_infile = "%s/%s/results.json" % (args.report_dir,pangolin_dir)
			pangolin_result = get_pangolin_info(p_infile,bc)
		
		# nextclade result
		nextclade_result = 'NA'

		# 样本数据量
		aln_stat = "%s/%s_rawlib.ionstats_alignment.json" % (args.report_dir,bc)
		if os.path.exists(aln_stat):
			readsNum_meanLen = get_bc_reads_num(aln_stat)
		else:
			readsNum_meanLen = ['NA','NA']
		
		reads_num = readsNum_meanLen[0]
		
		# 平均长度
		read_mean_len = readsNum_meanLen[1]

		# mapped_reads_n,on_target,mean_depth,uni [均一性]
		cov_dir = get_first_run_plugin_result(args.report_dir,'SARS_CoV_2_coverageAnalysis')
		if len(cov_dir) == 0:
			# no SARS_CoV_2_coverageAnalysis plugin run
			cov_stat = ['NA','NA','NA','NA']
		else:
			cov_files = glob.glob("%s/%s/*.bc_summary.xls" % (args.report_dir,cov_dir))
			if len(cov_files) == 1:
				# exist one file
				cov_file = cov_files[0]
				cov_stat = get_cov_stat(cov_file,bc)
			else:
				# do not exists OR has more than one *.bc_summary.xls file
				cov_stat = ['NA','NA','NA','NA']

		# 均一性
		uniformity = cov_stat[3]

		# 比对Reads数
		mapped_reads_num = cov_stat[0]

		# Ontarget率
		on_target_pct = cov_stat[1]

		# 平均测序深度
		mean_depth = cov_stat[2]

		# 一致性序列信息
		cons_dir = get_first_run_plugin_result(args.report_dir,'generateConsensus')
		if len(cons_dir) == 0:
			# no generateConsensus plugin run
			cons_info = ['NA','NA','NA']
		else:
			cons_file = "%s/%s/results.json" % (args.report_dir,cons_dir)
			cons_info = get_cons_info(cons_file,bc)
		
		# 组装N比例
		cons_N_pct = cons_info[0]

		# 一致性序列变异位点个数
		cons_var_num = cons_info[1]

		# 一致性序列杂合SNP个数
		cons_het_snp_num = cons_info[2]

		# TVC变异位点个数
		tvc_dir = get_first_run_plugin_result(args.report_dir,'SARS_CoV_2_variantCaller')
		if len(tvc_dir) == 0:
			# no SARS_CoV_2_variantCaller plugin run
			tvc_var_num = 'NA'
		else:
			tvc_file = "%s/%s/results.json" % (args.report_dir,tvc_dir)
			tvc_var_num = get_tvc_info(tvc_file,bc)

		# Q20碱基百分比
		q20_base_pct = 'NA'

		# Reads平均长度
		# pass

		# Pool1-Mean Reads per Amplicon
		reads_per_amp_p1 = ''

		# Pool2-Mean Reads per Amplicon
		reads_per_amp_p2 = ''
		
		# 是否提交
		if_submit = 'NA'

		# 提交日期
		submit_date = 'NA'

		# 备注
		note = 'NA'
		
		val = "\t".join(chef_date,seq_date,expName,report_name,chipType,total_reads,bc,sample_name,pangolin_result,nextclade_result,reads_num,uniformity,cons_N_pct,tvc_var_num,cons_var_num,cons_het_snp_num,q20_base_pct,read_mean_len,mapped_reads_num,on_target_pct,mean_depth,reads_per_amp_p1,reads_per_amp_p2,if_submit,submit_date,note)
		of.write(val+'\n')

	of.close()


if __name__ == '__main__':
	main()