import glob
from scipy import stats
import os

'''
Uses scipy stats to run T-test on the left and right half arrays of exon zscores
Some other statistical test should be done here because
1) zscores are most likely not uniformly distributed
2) breakpoints may be highly skewed to one side or another leading to array of 
zscores of size ~n*10^0 being compared to an array of size ~n*10^1

'''
def Ttests(cancertype):
	Totals = []
	for name in glob.glob("./Results/"+cancertype+"/*"):
		print name
		CancerList = []
		data = []
		left = []
		right = []
		Breakpoint = ''
		Dict = {}
		count = 0

		for line in open(name,'r'):

			line = line.split('\t')
			
			if "ENSG" in line[0]:

				Breakpoint = int(line[3])
				
				if len(left) != 0 and len(right)!= 0:
					data.append((saved_name,saved_orig,saved_oric,left,right))
					Dict[saved_id+'_'+str(count)] = data
					count += 1
				saved_id = line[0]
				saved_name = line[1]
				saved_orig = line[2]
				saved_oric = line[4].split('\n')[0]
				
				data = []
				left = []
				right = []

				
			else:
				try:
					exStart = int(line[1].split(',')[0])
					exStop = int(line[1].split(',')[1])

					if exStart < Breakpoint and exStop < Breakpoint:
						left.append(float(line[1].split(',')[2]))
					elif exStart > Breakpoint and exStop > Breakpoint:
						right.append(float(line[1].split(',')[2]))
					
				except Exception as e:
					pass


		Dict[saved_id+"_"+str(len(Dict.keys()))] = (saved_name,saved_orig,saved_oric,left,right)

		mod = 0
		best_pvalue = 1.0
		fusion = []
		for k in sorted(Dict.keys(), key = lambda x: int(x.split('_')[1])):
			name =  Dict[k][0][0]
			left = Dict[k][0][-2]
			right = Dict[k][0][-1]
			fusion.append((name,cancertype))
			
			try:
				best_pvalue = min(best_pvalue, float(stats.ttest_ind(left,right)[1]))
			except Exception:
				continue
			mod +=1 
			if mod%2 == 0:
				fusion.append(best_pvalue)
				Totals.append(fusion)
				fusion = []
				best_pvalue = 1.0

			
	return Totals


#build individual 
def pvalueStats():
	Cancers = ["BLCA"
	,"LUAD"
	,"BRCA"
	,"CESC"
	,"CLLE"
	,"COAD"
	,"DLBC"
	,"GBM"
	,"HNSC"
	,"KICH"
	,"KIRC"
	,"KIRP"
	,"LGG"
	,"LUSC"
	,"PRAD"
	,"THCA"
	,"MALY"
	,"LIRI"

	]
	#this takes a long time, comment out cancers you don't want to include
	os.system("mkdir ./Pvalues")
	for cname in Cancers:
		Totals = Ttests(cname)
		with open("./Pvalues/"+cname+"_pvalues.tsv",'w') as pvalues:
			for t in sorted(Totals,key = lambda x: x[2]):
				pvalues.write(t[0][0]+'\t'+t[1][0]+'\t'+str(t[2])+'\t'+t[0][1]+'\n')


## collect the best overall scores
def MasterList():
	masterList = []
	count = 0
	for pfile in glob.glob("./Pvalues/*"):
		for line in open(pfile,'r'):
			masterList.append(line.split('\t'))
	with open("Complete_PVlaues_Results.txt",'w') as temp:
		for i in sorted(masterList,key = lambda x:x[2]):
			if float(i[2]) <= .05:
				count += 1

			temp.write(i[0]+' '+i[1]+' '+i[2]+' '+i[3])
