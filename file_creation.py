import glob
import numpy as np
from scipy import stats
import Genefinder
from NMinusOneGroup import *
import summary_matching
import re
import pandas as pd
import os


#Matches the file prefixes for tumor files and rna data files to their
#donor ids and creates reference files
def generateDonorData(SummaryFile,patient_fusion_directory,RNASummary,RNAData):
	D1 = summary_matching.patientDonorId(SummaryFile,patient_fusion_directory)
	D2 = summary_matching.RNADonorId(RNASummary,RNAData)

	os.system("mkdir DonorData")
	with open("DonorData/PatientDonorIds.tsv",'w') as P1:
		for k in sorted(D1.keys()):
			P1.write(k+'\t'+D1[k][0]+'\t'+D1[k][1]+'\n')

	with open("DonorData/TumorIds.tsv",'w') as P2:
		for k in sorted(D2.keys()):
			P2.write(k+'\t'+D2[k][0]+'\t'+D2[k][1]+'\n')





#Take the previously held patient_fusion files and append to their filename
#the column from the RNA data that the fileheader exists in.
#When getting each patient's scores from the ~3GB RNA Data file,
#we'll know which scores to assign to which patient

def GroupCancers():
	exonValuesFileHeader = ExonScores.readline()
	Cancers = []
	os.system("mkdir CancerType_Groups")
	PDID = open("./DonorData/PatientDonorIds.tsv",'r')
	TID = open("./DonorData/TumorIds.tsv",'r')
	patientDonorId = {}
	tumorDonorId = {}
	for line in PDID:
		patientDonorId[line.split('\t')[1]] = (line.split('\t')[0],line.split('\t')[2].split('\n')[0])
	for line in TID:
		tumorDonorId[line.split('\t')[1]] = (line.split('\t')[0],line.split('\t')[2].split('\n')[0])


	patientfilenames = []
	for name in glob.glob("./patient_fusions_with_exons/*"):
	 	patientfilenames.append(name.split("exons/")[1].split('.')[0])

	tumorType = []

	inspect = [v[0] for v in patientDonorId.values() if v[0] in [x[0] for x in tumorDonorId.values()] ]

	for name in glob.glob("./patient_fusions_with_exons/*"):
		
		fileheader = name.split('exons/')[1].split('.')[0]
		try:
			patientDonorId[fileheader]
		except Exception:
			continue

		if patientDonorId[fileheader][0] in inspect:
			cancerType = patientDonorId[fileheader][1].split('-')[0]
		
			filename = open(name,'r')
			column = 0
			
			for filenameColumnMatch,j in enumerate(exonValuesFileHeader.split('\t')):
				if j in tumorDonorId.keys():
					if patientDonorId[fileheader][0]==tumorDonorId[j][0]:
						column = filenameColumnMatch
						break

			
			if cancerType not in Cancers:
				os.system("mkdir CancerType_Groups/"+cancerType)
			os.system("cp "+name+" CancerType_Groups/"+cancerType+"/"+fileheader+"_"+str(column)+".tsv")
			Cancers.append(cancerType)



def buildPatientDictionary(cancerType,file):

	filename = open(file,'r')
	patientDict = {}
	GeneA = 0
	GeneB = 0
	countGenes = 0
	start = 0
	exonsToSort = []
	seen_exon_numbers = np.zeros(100)
	for i,line in enumerate(filename):
		if line.split('\t')[0] == 'geneA':
			prevLine = i
		if i == prevLine+1:
			GeneA = line.split('\t')[0]
			GeneB = line.split('\t')[1].split('\n')[0]
		else:
			line = line.split('\t')
			

			if len(line)>2:
				if line[1] == GeneA:
					exonsToSort = []
					exonsToSort.append((line[0],line[1],line[2],int(line[3].split('\n')[0])))
				
				if line[1] == GeneB:
					exonsToSort = []
					exonsToSort.append((line[0],line[1],line[2],int(line[3].split('\n')[0])))
			if "chr" in line[0]:
				exonsToSort.append(line[1].split('\n')[0])
				countGenes+=1
				start = i+1
		if i >= start:
			try:
				int(line[0])
				
				if seen_exon_numbers[int(line[0])]!=1:
					exonsToSort.append((int(line[0]),int(line[1]),int(line[2].split('\n')[0])))
					seen_exon_numbers[int(line[0])] = 1
				else:
					pass
			except Exception as e:
				pass
		if line[0]=='\n':

			patientDict[countGenes]= sorted(exonsToSort, key = lambda x: x[0])
			seen_exon_numbers = np.zeros(100)
	return patientDict


#Here's where the statistical analysis needs to be changed.
#This is also the slowest part
#Deals with exon scores file which is ~3GB
def EvaluatePatientGroup(cancerType):
	checklist = []
	columnsToExamine = []
	patientfiles = []
	for name in glob.glob("./CancerType_Groups/"+cancerType+"/*"):
		fileheader = name.split("Groups/"+cancerType+"/")[1].split("_")[0]
		patientfiles.append(fileheader)
		column = int(name.split("Groups/"+cancerType+"/")[1].split("_")[1].split(".")[0])
		columnsToExamine.append((column,fileheader))
	
	scores = []
	print "Building Results for ",cancerType
	for name in glob.glob("./CancerType_Groups/"+cancerType+"/*"):
		patientDict = buildPatientDictionary(cancerType,name)
		if len(patientDict.keys())%2 != 0: continue
		fileheader = name.split(cancerType+"/")[1]
		main = fileheader.split("_")[0]
		print 'analzying file:',fileheader
		if cancerType not in checklist:
			os.system("mkdir ./Results/"+cancerType+"/")
			checklist.append(cancerType)
		with open("./Results/"+cancerType+"/"+name.split(cancerType+"/")[1].split("_")[0]+".tsv",'w') as FileWithScore:
			for k,v in patientDict.iteritems():

				patientsGene = int(re.sub('[^0-9]','', v[-1][0]))

				FileWithScore.write(v[-1][0]+'\t'+v[-1][1]+'\t'+v[-1][2]+'\t'+str(v[-1][3])+'\t'+v[-2]+'\n')
				#see NMinusOneGroup.py for description of UniqueGenes
				CompareTo = UniqueGenes("./CancerType_Groups/"+cancerType+"/*",cancerType)
				columnsForPatientGene = [int(i[0].split('_')[1].split('.')[0]) for i in CompareTo[fileheader] if i[1] == patientsGene]
				
				patientColumn = int(fileheader.split("_")[1].split('.')[0])
				ExonScores.seek(0)
				ExonScores.readline()
				scores = []

				'''
				Looks at the RNA Scores for each exon across 
				all patients of similar cancer types and calculates each exon's respective z score
				and stores them in new files under ./Results directory

				each line under the specific gene in the new .tsv file is saved as
				start-position, end-position, z score

				Several problems with this
				1) splicing leads to multiple repeat exon numbers
				Within the genecode annotation file, the exon_ids are unique, but the 
				RNA scores are only listed by exon_number, which isn't unique.
				-- for the exon numbers used, just the first occurrence of that number was stored
				-- more reading about the file needs to be done to know which one / how to use it
				-- The exon finding/storing is done in GeneFinder.geneIdRef()
				'''


				for line in ExonScores:
					gene_id = line.split('\t')[0].split('.')[0]
					exonNumber = int(line.split('\t')[0].split(':')[1])
					gene_id = int(re.sub('[^0-9]','', gene_id))
					
					if gene_id > patientsGene:break
					if exonNumber > max([int(x[0]) for x in v[0:-2]]) and gene_id == patientsGene:break
					if gene_id == patientsGene:
						scores.append(int(line.split('\t')[patientColumn]))
					
						for x in columnsForPatientGene:
							scores.append(int(line.split('\t')[x]))
				
						exonsPairs = [(int(exon[1]),int(exon[2])) for exon in v[0:-2] if exon[0] == exonNumber]

						Zscores = stats.mstats.zscore(scores)[0:len(exonsPairs)]
						
						for x,z in zip(exonsPairs,Zscores):
							
							FileWithScore.write(str(exonNumber)+'\t'+str(x[0])+','+str(x[1])+','+str(z)+'\n')
					
						scores = []
	
def handleFileBuilding(ExonAnnotationFile_input,ExonScores_input,patient_fusion_directory_input,
	RNAData_input, SummaryFile_input, RNASummary_input):

	ExonAnnotationFile = ExonAnnotationFile_input
	ExonScores = ExonScores_input
	patient_fusion_directory = patient_fusion_directory_input
	RNAData = RNAData_input
	SummaryFile = SummaryFile_input
	RNASummary = RNASummary_input
	
	generateDonorData(SummaryFile,patient_fusion_directory,RNASummary,RNAData)
	Genefinder.RNADonorId(ExonAnnotationFile,patient_fusion_directory)

	GroupCancers()


def BuildResults():
	os.system("mkdir Results")

	Cancers = [name.split("Groups/")[1] for name in glob.glob("./CancerType_Groups/*") ]
	for cname in Cancers:
		EvaluatePatientGroup(cname)
	#EvaluatePatientGroup('LIRI')
	#EvaluatePatientGroup('BLCA')
	#EvaluatePatientGroup('LUAD')
	#EvaluatePatientGroup('BRCA')
	#EvaluatePatientGroup("LIRI")

