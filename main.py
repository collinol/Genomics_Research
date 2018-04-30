
import NMinusOneGroup
import Evaluate
import Genefinder
import file_creation
import combine_patient_data_hgRef
import summary_matching


ExonAnnotationFile = open("gencode.v19.annotation.gtf_withproteinids","r")
ExonScores = open("counts_all_2015_combined_exon_dex.STAR.whitelist.tsv", "r")
patient_fusion_directory = "Referenced_patient_data"
RNAData = open("counts_all_2015_combined_exon_dex.STAR.whitelist.tsv", "r")
SummaryFile = open("pcawg_summary.tsv","r")
RNASummary = open("./exonRefs/pcawg.rnaseq.extended.metadata.aliquot_id.V4.tsv","r")
hgrefFile = "refGene_hg19_sorted.txt"



print "enter step number to start from"
print "if you have a directory called CancerType_Groups/ begin from step 3"
step_number = raw_input()

def options(step_number):
	step_number = int(step_number)
	if step_number == 1:
		#Step 1 - This only needs to be run once
		print "step 1"
		combine_patient_data_hgRef.build(hgrefFile)
		print "step 1 complete"
		#Step 2 - Run Once
		print "step 2"
		file_creation.handleFileBuilding(ExonAnnotationFile,ExonScores,patient_fusion_directory,
			RNAData, SummaryFile, RNASummary) 
		print "step 2 complete"
		#Step 3 - Needs to be run for each cancer type 
		#build UI for selecting cancers to build
		print "step 3"
		file_creation.BuildResults(ExonScores)
		print "step 4a"
		 
		#Step 4a - will run for all existing cancer types that were build in step 3
		#Builds list of best test results for specific cancer
		print "step 4a complete"
		Evaluate.pvalueStats()
		print "step 4b complete"
		#Step 4b - will aggregate all individual results into one "Best scores" file
		Evaluate.MasterList()
		print "complete"

	elif step_number == 2:
		print "step 2"
		file_creation.handleFileBuilding(ExonAnnotationFile,ExonScores,patient_fusion_directory,
			RNAData, SummaryFile, RNASummary) 
		print "step 2 complete"
		print "step 3"
		file_creation.BuildResults(ExonScores)
		print "step 3 complete"
		print "step 4a"
		Evaluate.pvalueStats()
		print "step 4a complete"
		print "step 4b"
		Evaluate.MasterList()
		print "complete"

	elif step_number == 3:
		print "step 3"
		file_creation.BuildResults(ExonScores)
		print "step 3 complete"
		print "step 4a"
		Evaluate.pvalueStats()
		print "step 4a complete"
		print "step 4b"
		Evaluate.MasterList()
		print "complete"

	elif step_number == 4:
		print "step 4a"
		Evaluate.pvalueStats()
		print "step 4a complete"
		print "step 4b"
		Evaluate.MasterList()
		print "complete"

options(step_number)