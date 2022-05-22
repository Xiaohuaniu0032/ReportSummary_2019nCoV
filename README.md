## Summary `Ion AmpliSeq SARS-CoV-2 Insight Research` results into one xls file.

## Output file header line
	建库日期 : chefDate
	测序日期 : seqDate
	expName : expName
	报告名称 : reportName
	芯片类型 : chipType 
	Barcode : Barcode
	样本名 : sampleName
	Pangolin分型 : Pangolin
	Nextclade分型 : Nextclade
	样本数据量 : totalReads
	均一性 : Uniformity
	组装N比例 : consensusN
	TVC变异位点个数 : tvcVarNum
	一致性序列变异位点个数 : consVarNum
	一致性序列杂合SNP个数 : consHetSnpNum
	Reads平均长度 : readMeanLength
	平均测序深度 : meanDepth
	Pool1-Mean Reads per Amplicon : Pool1-Mean Reads per Amplicon
	Pool2-Mean Reads per Amplicon : Pool2-Mean Reads per Amplicon
	是否提交 : ifSubmit
	提交日期 : submitDate
	备注 : Note


```
chefDate	seqDate	expName	reportName	chipType	Barcode	sampleName	Pangolin	Nextclade	totalReads	Uniformity	consensusN	tvcVarNum	consVarNum	consHetSnpNum	readMeanLength	meanDepth	Pool1-Mean Reads per Amplicon	Pool2-Mean Reads per Amplicon	ifSubmit	submitDate	Note
NA	NA	R_2020_01_20_09_37_25_user_GSS5PR-0070-71-20200117-Ampliseq_Flu	2019-nCoV-map2hg19-exon-virus_241	510	IonXpress_001	001 FluB	B.1.617.2	NA	NA	NA	0.4264	49	49	0	NA	NA	28293	25927	NA	NA	NA
NA	NA	R_2020_01_20_09_37_25_user_GSS5PR-0070-71-20200117-Ampliseq_Flu	2019-nCoV-map2hg19-exon-virus_241	510	IonXpress_003	002 FluB	B.1.617.2	NA	NA	NA	0.3021	45	49	0	NA	NA	3	2	NA	NA	NA

```

## Plugin example

`SARS_CoV_2_coverageAnalysis_out.1212` means that in this run, `ReportSummary_2019nCoV_out.1868`, the `SARS_CoV_2_coverageAnalysis_out.1212` results were used.

![ReportSummary_2019nCoV.png](https://github.com/Xiaohuaniu0032/ReportSummary_2019nCoV/blob/master/ReportSummary_2019nCoV.png)
