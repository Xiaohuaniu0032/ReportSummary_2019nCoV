import os
import sys
import argparse

def parse_args():
	AP = argparse.ArgumentParser("Collect 2019nCoV info")
	AP.add_argument('-report_dir',help='report dir [Auto_xxx]',dest='report_dir')
	AP.add_argument('-outdir',help='output dir [plugin_out/ReportSummary_2019nCoV.xxx/]',dest='outdir')
	return AP.parse_args()


def main():
	args = parse_args()
	report_name = os.path.basename(args.report_dir)
	outfile = "%s/%s.summary.xls" % (args.outdir,report_name)
	of = open(outfile,'w')
	
	of.close()



'''
output file header line

建库日期
测序日期
expName
报告名称
芯片类型
芯片总数据量
Barcode
样本名
Pangolin分型
Nextclade分型
样本数据量
均一性
组装N比例
TVC变异位点个数
一致性序列变异位点个数
一致性序列杂合SNP个数
Q20碱基百分比
Reads平均长度
比对Reads数
Ontarget率
平均测序深度
Pool1-Mean Reads per Amplicon
Pool2-Mean Reads per Amplicon
是否提交国家疾控
提交日期
备注
'''

my ($report_dir,$this_plugin_out_dir) = @ARGV;

my $report_name = basename($report_dir);
my $outfile = "$this_plugin_out_dir/$report_name\.ReportSummary.xls";
open O, ">$outfile" or die;
print O "建库日期\t测序日期\texpName\t报告名称\t芯片类型\t芯片总数据量\tBarcode\t样本名\tPangolin分型\tNextclade分型\t样本数据量\t均一性\t组装N比例\tTVC变异位点个数\t一致性序列变异位点个数\t一致性序列杂合SNP个数\tQ20碱基百分比\tReads平均长度\t比对Reads数\tOntarget率\t平均测序深度\tPool1-Mean Reads per Amplicon\tPool2-Mean Reads per Amplicon\t是否提交国家疾控\t提交日期\t备注\n";
# status.txt [1 for complete]
# serialized_S5yanzheng-20211223-chip2-MeanAccuracy_v2.json

# 建库日期？
# 测序日期？
# expName
# 报告名称
# 芯片类型
# 总数据量

# For each barcode