import evaluate
import filecreation
import combinepatientdatahgRef

if __name__ == '__main__':

    ExonAnnotationFile = open("./source_files/gencode.v19.annotation.gtf_withproteinids", "r")
    ExonScores = open("./source_files/counts_all_2015_combined_exon_dex.STAR.whitelist.tsv", "r")
    patient_fusion_directory = "referenced_patient_data"
    RNAData = open("./source_files/counts_all_2015_combined_exon_dex.STAR.whitelist.tsv", "r")
    SummaryFile = open("./source_files/pcawg_summary.tsv", "r")
    RNASummary = open("./exon_refs/pcawg.rnaseq.extended.metadata.aliquot_id.V4.tsv", "r")
    hgrefFile = "./source_files/refGene_hg19_sorted.txt"

    print "enter step number to start from"
    print "if you have a directory called cancer_type_groups/ begin from step 3"


    def options(step_number):
        step_number = int(step_number)
        if step_number == 1:
            # Step 1 - This only needs to be run once
            print "step 1"
            combinepatientdatahgRef.build(hgrefFile)
            print "step 1 complete"
            # Step 2 - Run Once
            print "step 2"
            filecreation.handleFileBuilding(ExonAnnotationFile, ExonScores, patient_fusion_directory,
                                            RNAData, SummaryFile, RNASummary)
            print "step 2 complete"
            # Step 3 - Needs to be run for each cancer type
            # build UI for selecting cancers to build
            print "step 3"
            filecreation.BuildResults(ExonScores)
            print "step 4a"

            # Step 4a - will run for all existing cancer types that were build in step 3
            # Builds list of best test results for specific cancer
            print "step 4a complete"
            evaluate.pvalueStats()
            print "step 4b complete"
            # Step 4b - will aggregate all individual results into one "Best scores" file
            evaluate.MasterList()
            print "complete"

        elif step_number == 2:
            print "step 2"
            filecreation.handleFileBuilding(ExonAnnotationFile, ExonScores, patient_fusion_directory,
                                            RNAData, SummaryFile, RNASummary)
            print "step 2 complete"
            print "step 3"
            filecreation.BuildResults(ExonScores)
            print "step 3 complete"
            print "step 4a"
            evaluate.pvalueStats()
            print "step 4a complete"
            print "step 4b"
            evaluate.MasterList()
            print "complete"

        elif step_number == 3:
            print "step 3"
            filecreation.BuildResults(ExonScores)
            print "step 3 complete"
            print "step 4a"
            evaluate.pvalueStats()
            print "step 4a complete"
            print "step 4b"
            evaluate.MasterList()
            print "complete"

        elif step_number == 4:
            print "step 4a"
            evaluate.pvalueStats()
            print "step 4a complete"
            print "step 4b"
            evaluate.MasterList()
            print "complete"


    options(raw_input())
