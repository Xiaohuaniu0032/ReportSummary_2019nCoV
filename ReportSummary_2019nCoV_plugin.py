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
	all_bams = glob.glob("%s/*_rawlib.bam")
	bams = []
	for bam in all_bams:
		name = os.path.basename(bam)
		bams.append(bam)

	bams_sort = sort(bams)
	return(bams_sort)

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
		max_len   = json_str['full']['max_read_length']
		#total_base_num = json_str['full']['num_bases']
		#q20_base_num = ''

	return([reads_num,mean_len,max_len])

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
		return(first_name)

def get_uniformity(infile): # plugin_out/SARS_CoV_2_coverageAnalysis_out.xxx/*.bc_summary.xls
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
				mapped_reads_n = vals[2]
				on_target = vals[4]
				mean_depth = vals[-2]
				uni = vals[-1]

	return([mapped_reads_n,on_target,mean_depth,uni])

def reads_per_pool(infile):
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

	return(avg_p1,avg_p2)

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
	header = "\t".join("建库日期","测序日期","expName","报告名称","芯片类型","芯片总数据量","Barcode","样本名","Pangolin分型","Nextclade分型","样本数据量","均一性","组装N比例","TVC变异位点个数","一致性序列变异位点个数","一致性序列杂合SNP个数","Q20碱基百分比","Reads平均长度","比对Reads数","Ontarget率","平均测序深度","Pool1-Mean Reads per Amplicon","Pool2-Mean Reads per Amplicon","是否提交国家疾控","提交日期","备注")
	of.write(header+'\n')
	of.close()

	# get libary date

	# get sequencing date

	# get expName

	# get report name
	report_name = os.path.basename(args.report_dir)

	# get chip type

	# get total reads number

	############################## for each barcode ##############################
	bams = get_bams(args.report_dir) # [IonXpress_003_rawlib.bam,]
	ion_params_00_json = os.path.join(args.report_dir,'ion_params_00.json')
	
	# check ion_params_00.json exists
	if os.path.exists(ion_params_00_json):
		print(ion_params_00_json)
		continue
	else:
		sys.exit('[Error: can not find ion_params_00.json file, will exit]')
	
	all_barcodes = barcode_to_sampleName(ion_params_00_json)

	for bc in all_barcodes:
		# IonCode_0303
		sample_name = all_barcodes[bc] # maybe None
		
		# check pangolin plugin
		pangolin_result = ''
		
		nextclade_result = 'NA'

		aln_stat = "%s/%s_rawlib.ionstats_alignment.json" % (args.report_dir,bc)
		
		# 
		uniformity = get_uniformity()
		cons_N_pct = 
		pool1_avg_reads_per_amplicon = ''
		pool2_avg_reads_per_amplicon = ''
