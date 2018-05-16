# Genomics_Research

## Step 1
```
hgrefFile = "refGene_hg19_sorted.txt" #prebuilt file  
combine_patient_data_hgRef.build(hgrefFile)
```
first 5 lines of hg19:
```
585	NR_024540	chr1	-	14361	29370	29370	29370	11	14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,	14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,	0	WASH7P	unk	unk	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  
585	NR_026818	chr1	-	34610	36081	36081	36081	3	34610,35276,35720,	35174,35481,36081,	0	FAM138A	unk	unk	-1,-1,-1,  
585	NR_026820	chr1	-	34610	36081	36081	36081	3	34610,35276,35720,	35174,35481,36081,	0	FAM138F	unk	unk	-1,-1,-1,  
585	NM_001005484	chr1	+	69090	70008	69090	70008	1	69090,	70008,	0	OR4F5	cmpl	cmpl	0,  
587	NR_028322	chr1	+	323891	328581	328581	328581	3	323891,324287,324438,	324060,324345,328581,	0	LOC100132287	unk	unk	-1,-1,-1,
```
Description of columns can be read up on here  
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=genes&hgta_track=refGene&hgta_table=refGene&hgta_doSchema=describe+table+schema  
We only need the geneId, chromosome number, strand, transcription start and stop points  
columns 1,2,3,4,5 respectively  

combine_patient_data_hgRef.build() uses the prebuilt patient files (inside the directory "patient_data"), labeled as   
```
ff870342-f0d6-4450-8f9c-344c046a0baf.pcawg_consensus_1.6.161116.somatic.sv.bedpe
ffa976f0-aa60-4867-842e-361afa7d68ac.pcawg_consensus_1.6.161116.somatic.sv.bedpe
ffad9288-c622-11e3-bf01-24c6515278c0.pcawg_consensus_1.6.161116.somatic.sv.bedpe
.
.
.
```  
(those file prefixes (up to .pcawg) will be used later)  
and we create a new directory, called "Referenced_patient_data" inside your current working directory. 
After this process runs (for the size of this data, it takes approximately 10-12 minutes), the "Referenced_patient_data" directory will contain individual files with the same file prefix as the corresponding patient data, but with only the valid gene fusions listed. For example:  
```
first 5 lines (of ~50 lines) of 0bd3a230-531e-44aa-9999-b8ed8da0176b.pcawg_consensus_1.6.161116.somatic.sv.bedpe
chrom1	start1	end1	chrom2	start2	end2	sv_id	pe_support	strand1	strand2	svclass	svmethod
11	82798867	82798868	11	82936093	82936094	SVMERGE10	37	+	-	DEL	SNOWMAN_BRASS_dRANGER_DELLY
12	115046363	115046364	4	56334945	56334946	SVMERGE34	35	-	+	TRA	SNOWMAN_BRASS_dRANGER_DELLY
12	115046403	115046404	4	72591306	72591307	SVMERGE42	32	+	-	TRA	SNOWMAN_BRASS_dRANGER_DELLY
12	118147063	118147064	13	54442327	54442328	SVMERGE38	12	-	-	TRA	SNOWMAN_BRASS_DELLY
```
produces 
```
geneA	geneA orientation	geneB	geneB orientation	chr1	position	chr1 orientation	chr2	position	chr2 orientation
./patient_data/0bd3a230-531e-44aa-9999-b8ed8da0176b.pcawg_consensus_1.6.161116.somatic.sv.bedpe									
KSR2	-	PCDH9	-	chr12	118147303	+	chr13	67542094	-
CLDN14	-	DYRK1A	+	chr21	37873660	-	chr21	38850309	-
DYRK1A	+	TMPRSS2	-	chr21	38846232	-	chr21	42868379	-
DYRK1A	+	TMPRSS2	-	chr21	38846936	+	chr21	42869068	+
DYRK1A	+	ERG	-	chr21	38850232	+	chr21	39893645	+
AFF1	+	INTS12	-	chr4	88026639	-	chr4	106603789	-
INTS12	-	KIAA1429	-	chr4	106603808	+	chr8	95541749	-
```
That's the entire file. First line column headers, second line is the original source file, followed by the data  

## Step 2
