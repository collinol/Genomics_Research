import glob
import os

##Match original patient file headers to cancer file headers via metadata
def prepExon(filename):
	exonIds = []
	ExonFile = filename
	exonheader = ExonFile.readline()
	for i in exonheader.split('\t'):
		exonIds.append (i)
	return exonIds


def prepPatients(directory):
	patientfusionIds = []
	for name in glob.glob("./"+directory+"/*"):
		patientfusionIds.append(name.split("data/")[1].split(".")[0])
	return patientfusionIds



def patientDonorId(PatientSummary,fromDirectory):
	patientfusionIds = prepPatients(fromDirectory)
	DonorIds_pcawg = {}
	SummaryFile = PatientSummary
	SummaryFile.seek(0)
	count = 0
	print "Matching patient files to donor id"
	for fusionid in patientfusionIds:
		if count%100 == 0:
			print int(100*float(count)/len(patientfusionIds)),"%"

		for line in SummaryFile:
			line = line.split('\t')
			
			for i in [line[27],line[41],line[46],line[51],line[56],line[60],line[25]]:

				if len(i.split('-')) > 1 or len(i.split(',')) > 1:
					if fusionid == i or fusionid in i.split(","):
						if line[4] not in DonorIds_pcawg:
							DonorIds_pcawg[line[4]] = (fusionid,line[2])
		SummaryFile.seek(0)
		count += 1
	
	return DonorIds_pcawg


def RNADonorId(RNASummary,fromRNAData):
	exonref = RNASummary
	exonIds = prepExon(fromRNAData)
	DonorIds_exons = {}
	count = 0
	print "Matching RNA summary data to donor id"
	for column,fusionid in enumerate(exonIds):
		if count%100 == 0:
			print int(100*float(count)/len(exonIds)),"%"
		exonref.readline()
		for line in exonref:
			line = line.split('\t')

			idformat = line[17]
			if fusionid == idformat:
				
				if line[3] not in DonorIds_exons:
				
					DonorIds_exons[line[3]] = (fusionid, line[2])
			
		exonref.seek(0)
		count += 1
	return DonorIds_exons
