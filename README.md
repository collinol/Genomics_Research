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
