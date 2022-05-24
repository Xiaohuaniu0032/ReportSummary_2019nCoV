# -*- coding: UTF-8 -*-
import os
import sys
import argparse
import re
import json
import glob
import io
import codecs

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

def get_loading_info(infile):
	'''
	usable reads
	loading rate
	enrichment

	/Auto_xxx/serialized_Auto_user_GSS5PR-0070-144-ludaopei-530chip2_502.json
	'''

	with open(infile,'r') as json_file:
		json_str = json.load(json_file)
		addressable = json_str['fields'].get('adjusted_addressable','NA')
		final_usable_reads = json_str['fields'].get('libFinal','NA')
		
		load_rate = json_str['fields'].get('loading','NA') # bead / adjusted_addressable
		
		# enrichment = live / bead
		live = int(json_str['fields'].get('live')) # 34182286
		bead = int(json_str['fields'].get('bead')) # 34183815
		
		if bead != 0:
			enrichment = float(live/bead)
		else:
			enrichment = 'NA'

		# live = lib_live + tf_live
		lib_live = json_str['fields'].get('lib') # 33914533
		tf_live  = json_str['fields'].get('tf')  # 267753

	return([final_usable_reads,load_rate,enrichment])

def get_filter_info(infile):
	'''
	/Auto_xxx/basecaller_results/BaseCaller.json
	
	return values:
		poly clonal
		primer dimer
		low quality
	'''
	with open(infile,'r') as json_file:
		json_str = json.load(json_file)
		poly         = json_str['Filtering']['LibraryReport'].get('filtered_polyclonal','NA') # 9284994
		primer_dimer = json_str['Filtering']['LibraryReport'].get('filtered_primer_dimer','NA') # 15865
		low_qual     = json_str['Filtering']['LibraryReport'].get('filtered_low_quality','NA') # 3797720

	return([poly,primer_dimer,low_qual])

def get_analysis_date(infile):
	# expMeta.dat
	'''
	Run Name = R_2022_03_10_17_29_05_user_GSS5PR-0070-144-ludaopei-530chip2
	Run Date = 2022-03-10 09:33:46+00:00
	Run Flows = 800
	...
	Instrument = GSS5PR-0070
	Flow Order = TACGTACGTCTGAGCATCGATCGATGTACAGC
	Analysis Date = 2022-03-10
	Analysis Flows = 0
	runID = D2100
	'''
	exp_info = {}
	with open(infile,'r') as exp_mata:
		for line in exp_mata:
			vals = line.strip().split('=')
			k = vals[0].strip()
			v = vals[1].strip()
			exp_info[k] = v

	seq_date = exp_info.get('Analysis Date','NA')
	return(seq_date)



def get_basic_info(infile):
	'''
	get below info from 'ion_params_00.json' file
		* seq time
		* expName
		* chip type
	'''
	with open(infile,'r') as json_file:
		json_str = json.load(json_file)
		seq_date = json_str['log'].get('start_time','NA')
		expName = json_str.get('expName','NA')
		chipType = json_str['exp_json'].get('chipType','NA')

	return([seq_date,expName,chipType])

def barcode_to_sampleName(infile): # ion_params_00.json
	bc2name = {}
	with open(infile,'r') as json_file:
		json_str = json.load(json_file)
		#sss = json_str['experimentAnalysisSettings']['barcodedSamples']
		#print(sss)
		samples = json_str['experimentAnalysisSettings']['barcodedSamples'].keys()
		for sample in samples:
			name = sample
			bc   = json_str['experimentAnalysisSettings']['barcodedSamples'][name]['barcodes'][0]
			bc2name[bc] = name

	return(bc2name)

def get_bc_reads_num(infile): # IonCode_0301_rawlib.ionstats_alignment.json
	with open(infile,'r') as json_file:
		json_str = json.load(json_file)
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
	#print(plugin_name)
	start_json_list = glob.glob("%s/plugin_out/%s_out.*/startplugin.json" % (report_dir,plugin_name))
	if len(start_json_list) == 0:
		return([])
	else:
		xxx_int = []
		for file in start_json_list:
			plugin_full_path = os.path.dirname(file) # /results/analysis/output/Home/Auto_xxx/plugin_out/vcMerge_out.1759
			plugin_basename = os.path.basename(plugin_full_path) # vcMerge_out.1759
			xxx = plugin_basename.split('.')[1]
			# SARS_CoV_2_variantCaller_out.qinghaiCDC.20220403
			# this dir is for test purpose, and its name is not correct (just for manually test)
			# need to skip this dir
			if xxx[0].isdigit():
				xxx_int.append(xxx)
			else:
				continue

		xxx_int.sort()
		print(xxx_int)
		first_xxx = xxx_int[0]
		first_name = "%s_out.%s" % (plugin_name,first_xxx)
		return(first_name) # SARS_CoV_2_coverageAnalysis_out.xxx

def get_cov_stat(infile,barcode): # plugin_out/SARS_CoV_2_coverageAnalysis_out.xxx/*.bc_summary.xls
	'''
	Barcode ID      Sample Name     Mapped Reads    Filtered Reads  Target Reads    Mean Depth      Uniformity
	IonXpress_001   001 FluB        6578122 0.32%   98.04%  28274   97.10%
	IonXpress_003   002 FluB        635     0.00%   100.00% 4.5     49.97%
	IonXpress_004   003 FluB        53592   0.15%   99.56%  332.6   73.31%

	'''
	
	#print("args are: %s %s" % (infile,barcode))
	if_this_sample_exists = 0
	return_val = []
	with open(infile,'r') as cov_summary:
		for line in cov_summary:
			#print(line)
			if line.startswith('Barcode ID'):
				continue
			else: 
				vals = line.rstrip().split('\t')
				#print(vals[0])
				if vals[0] == barcode: # this sample line
					print(line)
					if_this_sample_exists = 1
					mapped_reads_n = vals[2]
					on_target = vals[4]
					mean_depth = vals[-2]
					uni = vals[-1]
					
					return_val.append(mapped_reads_n)
					return_val.append(mean_depth)
					return_val.append(uni)
	if if_this_sample_exists:
		# exists this sample info
		return(return_val)
	else:
		# not exists
		return_val = ['NA','NA','NA']
		return(return_val)

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
		lines = amplicon_info.readlines()
		for line in lines:
			#print(line)
			vals = line.split('\t')
			#print(vals)
			
			pool_info = vals[4].upper()
			# GENE_ID=LDLRAP1;POOL=1;CNV_HS=0;CNV_ID=LDLRAP1
			# GENE_ID=r1;Pool=2
			total_reads = int(vals[-6])
			for v in pool_info.split(';'):
				if re.match('POOL',v):
					if v == 'POOL=1':
						p1_num.append(total_reads)
					if v == 'POOL=2':
						p2_num.append(total_reads)
	print(p1_num)
	print(p2_num)			
	avg_p1 = sum(p1_num)/len(p1_num)
	avg_p2 = sum(p2_num)/len(p2_num)
	#avg_p1 = float(mean(p1_num),2)
	#avg_p2 = float(mean(p2_num),2)

	return([avg_p1,avg_p2])

def get_tvc_info(infile,barcode):
	'''
	variantCaller
	SARS_CoV_2_coverageAnalysis
	
	/plugin_out/SARS_CoV_2_variantCaller_out.1733/results.json
	'''
	
	with open(infile,'r') as json_file:
		json_str = json.load(json_file)
		all_bc = json_str['barcodes'].keys()
		print("all bc are: %s" % (all_bc))
		if barcode in all_bc:
			var_num = json_str['barcodes'][barcode]['variants']['variants']
			return(var_num)
		else:
			# do not contain this barcode in results.json
			var_num = 'NA'
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
		json_str = json.load(json_file)
		all_bc = json_str['barcodes'].keys()
		# check if barcode in all_bc
		if barcode in all_bc:
			pct_N = json_str['barcodes'][barcode].get('Percent N','NA')
			var_num = json_str['barcodes'][barcode].get('Variants','NA')
			het_snp = json_str['barcodes'][barcode].get('Het snps','NA')
		#het_indel = json_str['barcodes'][barcode]['Het indels']
			return([pct_N,var_num,het_snp])
		else:
			return(['NA','NA','NA'])

def get_pangolin_info(infile,barcode):
	'''
	plugin_out/SARS_CoV_2_lineageID_out.1553/results.json
	'''
	with open(infile,'r') as json_file:
		json_str = json.load(json_file)
		all_bc = json_str['barcodes'].keys()
		if barcode in all_bc:
			Lineage = json_str['barcodes'][barcode].get('Lineage','NA')
			return(Lineage)
		else:
			return('NA')

def main():
	args = parse_args()
	report_name = os.path.basename(args.report_dir)
	outfile = "%s/%s.summary.xls" % (args.outdir,report_name)
	of = io.open(outfile,'w',encoding='utf-8')
	#header = "\t".join(["建库日期","测序日期","expName","报告名称","芯片类型","Barcode","样本名","Pangolin分型","Nextclade分型","样本数据量","均一性","组装N比例","TVC变异位点个数","一致性序列变异位点个数","一致性序列杂合SNP个数","Reads平均长度","平均测序深度","Pool1-Mean Reads per Amplicon","Pool2-Mean Reads per Amplicon","是否提交","提交日期","备注"])
	header = "\t".join(['seqDate','expName','reportName','chipType','Barcode','sampleName','Pangolin','Nextclade','totalReads (0.5~1M)','Uniformity','consensusN (<1)','tvcVarNum','consVarNum','consHetSnpNum','readMeanLength (>200bp)','meanDepth','Pool1-Mean Reads per Amplicon','Pool2-Mean Reads per Amplicon','P1/P2_Ratio','Loading','Enrichment','Polyclonal','Low Quality','Adapter Dimer','ifSubmit','submitDate','Note'])
	#of.write(codecs.BOM_UTF8)
	of.write(header.decode('utf-8')+'\n')
	
	plugin_info = {}
	# make a log file, this file will record which plugin's result was used by this plugin
	'''
	SARS_CoV_2_coverageAnalysis	SARS_CoV_2_coverageAnalysis_out.1212
	SARS_CoV_2_variantCaller	SARS_CoV_2_variantCaller_out.1567
	generateConsensus		generateConsensus_out.1532		
	SARS_CoV_2_lineageID		SARS_CoV_2_lineageID_out.1553
	'''
	plugin_log = "%s/plugin.log" % (args.outdir)
	of_plugin = open(plugin_log,'w')
	
	ion_params_00_json = os.path.join(args.report_dir,'ion_params_00.json')
	# check ion_params_00.json exists
	if os.path.exists(ion_params_00_json):
		print(ion_params_00_json)
		pass
	else:
		sys.exit('[Error: can not find ion_params_00.json file, will exit]')

	# get libary date
	chef_date = 'NA'
	#print("chef date is: %s" % (chef_date)

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
	
	# first print all barcode and its name
	print("########## sample info is ##########")
	print("Barcode\tSampleName")
	for bc in sorted(all_barcodes.keys()):
		sample_name = all_barcodes[bc]
		print("%s\t%s" % (bc,sample_name))
	print('\n')

	for bc in sorted(all_barcodes.keys()):
		print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>stat %s..." % bc)
		# IonCode_0303
		sample_name = all_barcodes[bc] # maybe None
		
		# pangolin result
		print("check SARS_CoV_2_lineageID plugin results...")
		pangolin_dir = get_first_run_plugin_result(args.report_dir,'SARS_CoV_2_lineageID')
		if len(pangolin_dir) == 0:
			# no SARS_CoV_2_lineageID plugin run
			print("no SARS_CoV_2_lineageID plugin run")
			pangolin_result = 'NA'
			plugin_info['SARS_CoV_2_lineageID'] = 'NA'
			#of_plugin.write("SARS_CoV_2_lineageID = %s" % ('NA'))
		else:
			print("using %s result" % (pangolin_dir))
			p_infile = "%s/plugin_out/%s/results.json" % (args.report_dir,pangolin_dir)
			# if plugin executed but failed, then you will not find results.json
			print(p_infile)
			if os.path.exists(p_infile):
				pangolin_result = get_pangolin_info(p_infile,bc)
				plugin_info['SARS_CoV_2_lineageID'] = pangolin_dir
				#of_plugin.write("SARS_CoV_2_lineageID = %s" % (pangolin_dir))
			else:
				pangolin_result = 'NA'
				print("[Warning]:can not find SARS_CoV_2_lineageID results.json file, will skipped")
				plugin_info['SARS_CoV_2_lineageID'] = 'NA'
		print("pangolin result is: %s" % (pangolin_result))
		print('\n')

		# nextclade result
		print("check Nextclade_2019nCoV plugin results")
		print("skip Nextclade results")
		nextclade_result = 'NA'
		print('\n')

		# 样本数据量
		print("check sample reads num and mean_len info...")
		aln_stat = "%s/%s_rawlib.ionstats_alignment.json" % (args.report_dir,bc)
		if os.path.exists(aln_stat):
			print("find %s file" % (aln_stat))
			readsNum_meanLen = get_bc_reads_num(aln_stat)
		else:
			print("[Warning: can not find %s file" % (aln_stat))
			readsNum_meanLen = ['NA','NA']
		
		reads_num = readsNum_meanLen[0]
		
		# 平均长度
		read_mean_len = readsNum_meanLen[1]
		print("reads num is: %s" % (reads_num))
		print("men_len is: %s" % (read_mean_len))
		print("\n")

		print("check mean_depth and uniformity info...")
		# mean_depth,uni [均一性]
		cov_dir = get_first_run_plugin_result(args.report_dir,'SARS_CoV_2_coverageAnalysis')
		if len(cov_dir) == 0:
			# no SARS_CoV_2_coverageAnalysis plugin run
			cov_stat = ['NA','NA','NA','NA']
			plugin_info['SARS_CoV_2_coverageAnalysis'] = 'NA'
			#of_plugin.write("SARS_CoV_2_coverageAnalysis = %s" % ('NA'))
		else:
			print("using %s results" % (cov_dir))
			cov_files = glob.glob("%s/plugin_out/%s/*.bc_summary.xls" % (args.report_dir,cov_dir))
			if len(cov_files) == 1:
				# exist one file
				cov_file = cov_files[0]
				print(cov_file)
				cov_stat = get_cov_stat(cov_file,bc)
				plugin_info['SARS_CoV_2_coverageAnalysis'] = cov_dir
				#of_plugin.write("SARS_CoV_2_coverageAnalysis = %s" % (cov_dir))
			else:
				# do not exists OR has more than one *.bc_summary.xls file
				cov_stat = ['NA','NA','NA']
				plugin_info['SARS_CoV_2_coverageAnalysis'] = 'NA'
				#of_plugin.write("SARS_CoV_2_coverageAnalysis = %s" % ('NA'))

		#print(cov_stat)
		# 均一性
		uniformity = cov_stat[2]

		# 比对Reads数
		#mapped_reads_num = cov_stat[0]

		# Ontarget率
		#on_target_pct = cov_stat[1]

		# 平均测序深度
		mean_depth = cov_stat[1]
		
		print("mean depth is: %s" % (mean_depth))
		print("uniformity is: %s" % (uniformity))
		print("\n")

		print("check generateConsensus plugin results...")
		# 一致性序列信息
		cons_dir = get_first_run_plugin_result(args.report_dir,'generateConsensus')
		if len(cons_dir) == 0:
			# no generateConsensus plugin run
			print("no generateConsensus plugin run")
			cons_info = ['NA','NA','NA']
			plugin_info['generateConsensus'] = 'NA'
			#of_plugin.write("generateConsensus = %s" % ('NA'))
		else:
			print("using %s results" % (cons_dir))
			cons_file = "%s/plugin_out/%s/results.json" % (args.report_dir,cons_dir)
			print(cons_file)
			if os.path.exists(cons_file):
				cons_info = get_cons_info(cons_file,bc)
				plugin_info['generateConsensus'] = cons_dir
				#of_plugin.write("generateConsensus = %s" % (cons_dir))
			else:
				cons_info = ['NA','NA','NA']
				print("[Warnings: can not find generateConsensus results.json file, will skip")
				plugin_info['generateConsensus'] = 'NA'
				#of_plugin.write("generateConsensus = %s" % ('NA'))
		
		# 组装N比例
		cons_N_pct = cons_info[0]

		# 一致性序列变异位点个数
		cons_var_num = cons_info[1]

		# 一致性序列杂合SNP个数
		cons_het_snp_num = cons_info[2]
		
		print("cons_N_pct is: %s" % (cons_N_pct))
		print("cons_var_num is: %s" % (cons_var_num))
		print("cons_het_snp_num is: %s" % (cons_het_snp_num))

		print("\n")
		
		print("check SARS_CoV_2_variantCaller plugin results...")
		# TVC变异位点个数
		tvc_dir = get_first_run_plugin_result(args.report_dir,'SARS_CoV_2_variantCaller')
		if len(tvc_dir) == 0:
			# no SARS_CoV_2_variantCaller plugin run
			print("no SARS_CoV_2_variantCaller plugin run")
			tvc_var_num = 'NA'
			plugin_info['SARS_CoV_2_variantCaller'] = 'NA'
			#of_plugin.write("SARS_CoV_2_variantCaller = %s" % ('NA'))
		else:
			print("using %s results" % (tvc_dir))
			tvc_file = "%s/plugin_out/%s/results.json" % (args.report_dir,tvc_dir)
			print(tvc_file)
			if os.path.exists(tvc_file):
				tvc_var_num = get_tvc_info(tvc_file,bc)
				plugin_info['SARS_CoV_2_variantCaller'] = tvc_dir
				#of_plugin.write("SARS_CoV_2_variantCaller = %s" % (tvc_dir))
			else:
				tvc_var_num = 'NA'
				print("[Warning]: can not find SARS_CoV_2_variantCaller results.json file, will skip")
				plugin_info['SARS_CoV_2_variantCaller'] = 'NA'
				#of_plugin.write("SARS_CoV_2_variantCaller = %s" % ('NA'))
		
		print("tvc num is: %s" % (tvc_var_num))
		print("\n")
		# Q20碱基百分比
		#q20_base_pct = 'NA'

		# Reads平均长度
		# pass

		print("check avg reads num per pool info...")
		# reads_per_pool
		cov_dir = get_first_run_plugin_result(args.report_dir,'SARS_CoV_2_coverageAnalysis')
		if len(cov_dir) == 0:
			# no SARS_CoV_2_coverageAnalysis plugin run
			avg_reads_pool = ['NA','NA']
		else:
			print("using %s results" % (cov_dir))
			amplicon_cov_file = glob.glob("%s/plugin_out/%s/%s/*.amplicon.cov.xls" % (args.report_dir,cov_dir,bc))
			if len(amplicon_cov_file) == 1:
				# exist one file
				cov_file = amplicon_cov_file[0]
				print(cov_file)
				avg_reads_pool = reads_per_pool(cov_file)
			else:
				# do not exists OR has more than one *.bc_summary.xls file
				print("did not see *.amplicon.cov.xls file")
				avg_reads_pool = ['NA','NA']

		# Pool1-Mean Reads per Amplicon
		reads_per_amp_p1 = avg_reads_pool[0]

		# Pool2-Mean Reads per Amplicon
		reads_per_amp_p2 = avg_reads_pool[1]
		
		# p1/p2
		if reads_per_amp_p2 != 0:
			p1_vs_p2 = round(reads_per_amp_p1/float(reads_per_amp_p2),3)
		else:
			p1_vs_p2 = 'NA'

		print("reads_per_amp_p1 is: %s" % (reads_per_amp_p1))
		print("reads_per_amp_p2 is: %s" % (reads_per_amp_p2))
		print("reads_per_amp_p1 / reads_per_amp_p2 is: %s" % (p1_vs_p2))

		
		print("\n\n\n")
		# 是否提交
		if_submit = 'NA'

		# 提交日期
		submit_date = 'NA'

		# 备注
		note = 'NA'
			
		h = (str(chef_date),str(seq_date),expName,report_name,str(chipType),bc,sample_name,pangolin_result,nextclade_result,str(reads_num),str(uniformity),str(cons_N_pct),str(tvc_var_num),str(cons_var_num),str(cons_het_snp_num),str(read_mean_len),str(mean_depth),str(reads_per_amp_p1),str(reads_per_amp_p2),if_submit,submit_date,note)
		val = "\t".join(h)
		of.write(val.decode('utf-8')+'\n')
	of.close()
	
	plugin_list = ['SARS_CoV_2_coverageAnalysis','SARS_CoV_2_variantCaller','generateConsensus','SARS_CoV_2_lineageID']
	SARS_CoV_2_coverageAnalysis = plugin_info.get('SARS_CoV_2_coverageAnalysis','NA')
	SARS_CoV_2_variantCaller    = plugin_info.get('SARS_CoV_2_variantCaller','NA')
	generateConsensus           = plugin_info.get('generateConsensus','NA')
	SARS_CoV_2_lineageID        = plugin_info.get('SARS_CoV_2_lineageID','NA')
	of_plugin.write(SARS_CoV_2_coverageAnalysis+'\n')
	of_plugin.write(SARS_CoV_2_variantCaller+'\n')
	of_plugin.write(generateConsensus+'\n')
	of_plugin.write(SARS_CoV_2_lineageID+'\n')
	of_plugin.close()


if __name__ == '__main__':
	main()
