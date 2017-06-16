

###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### lucifanalysis.py
### Copyright 2016 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 2017-06-16\
### Version 2.0

import numpy as np
import os
import pandas as pd
import sys
from outliers import smirnov_grubbs as grubbs
import scipy as sp 
import argparse




########################################################################################################################################################
########################################################## USER INPUT ##################################################################################
cwd = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description="""Takes Dual Glow luciferase readings from a file, also containing the legend and the plate plan, 
							consolidates all replicates no matter their number into a report for each condition, and generate a PDF report with the figures""")
parser.add_argument("-lu","--luciffile",
					help="Luciferase readings in a .txt file. 'Defaults to docs/Input/template_AP170428.txt'", default= str(cwd)+"/docs/Input/template_AP170428.txt")
parser.add_argument("-g", "--grubbs", action="store_true",
					help="Use Grubbs Test to automatically remove potential outliers within replicates. Useful for high-throughput. However, it is not recommended for less than 6 replicates")
parser.add_argument("-p", "--pdf", action="store_true",
					help="Export figures to PDF")

args = parser.parse_args()
luciffile = args.luciffile
dogrubbs = args.grubbs
pdf = args.pdf

print "\nlucifanalysis.py script for automated analysis of luciferase experiments  \n"
print " ** use the -h flag for detailed help on how to use this script **\n"

print "Using %s for Input luciferase readings ..." % luciffile
if dogrubbs == True:
	print "Removing outliers using Grubbs test..."

print "\n\n"

########################################################################################################################################################
########################################################################################################################################################


def grubbsfunc(ratiodff,columns):
	""" Take a dataframe of RLU values for replicates, find any outliers based on the Grubb's test and return modified """
	ratiodict = {}
	for i in xrange(len(ratiodff)):
		dfempty = pd.DataFrame.transpose(pd.DataFrame(np.nan, index=[i+1], columns=columns ))
		data= ratiodff.loc[int(i+1)]
		
		if len(list(data)) > 2:  ##impossible to detect outliers in 2 or less. 
		
			dfempty.update(grubbs.test(data, alpha=0.05))
			z = dfempty.to_dict()
			ratiodict.update(z)
			
		else:
			dfempty.update(data)
			z = dfempty.to_dict()
			ratiodict.update(z)
	ratiodf = pd.DataFrame.from_dict(ratiodict).transpose()

	return ratiodf



#### Start here ## Open input file
rawtemp = open(luciffile, "rw+")
raw = []
for i in rawtemp.readlines():
	a = i
	a = a.replace("\n", "").replace("\r", "")
	raw.append(a)

########################Read METADATA
inputfile = raw[0].split("\t")[1]
expdate = raw[1].split("\t")[1]
cellline = raw[2].split("\t")[1]
outputfile = raw[3].split("\t")[1]


####################Section to retrieve lucif and renilla reads
fireflylist = []
renillalist = []

for i in raw[12:20]:
	splitlist = i.split("\t")
	fireflylist.append(splitlist[1:13])
	renillalist.append(splitlist[16:28])
ffnp = np.array(fireflylist) # make numpy arrays of firefly and renilla readings
rnnp = np.array(renillalist)


###############section to retrieve plan and legend
plan = []
for i in raw[22:]:
	plan.append(i)

planlist = []
for j in plan[2:10]: #Remove peripheral wells. You should never have actual readings in peripheral wells, because of the edge effect. 
	a = j.replace("X", "0")

 
	planlist.append(a.split("\t")[1:12])
plannp = np.array(planlist, dtype="int32") ## make an array for the plan
plandict = {}
#######legend

legtemp= plan[11:len(plan)]  
leg = {} 
for i in legtemp: 
	leg[int(i.split("\t")[0])] = i.split("\t")[1:6]

legendcols = ["CONDITION", "TF", "PROMOTER", "EXP_OR_CONTROL", "ControlID"]

legend = pd.DataFrame.from_dict(leg).T
legend.columns = legendcols  #legend dataframe 


########






### count how many each wells there are for each condition ID
unique,pos = np.unique (plannp,return_inverse=True)
counts = np.bincount(pos)
counts2 = np.delete(counts,0)
maxpos = counts2.argmax()
maxcount= counts2[maxpos] #highest number of wells for all conditions: will define number of columns for report. 
uniquelist = list(unique) #
uniquelist.pop(0) #remove 0 from list since, that includes all empty wells and will necesarily be the first hit. 


legdict = {}
legdictrn = {}

for i in uniquelist:
	legdict[i] = []
	legdictrn[i] = []
		
	indexcount = 0
	fflist = maxcount* list("-") # make lists of "-" of the maximum number or replicates any given condition has. 
	rnlist = maxcount * list("-") # If another condition has elss replicate, the next block of code will simply change the "-" for the actual values, and leave the rest unchanged. 
	
	for j in xrange(maxcount):
		try:
			fflist[indexcount] = ffnp[list(np.where(plannp == i))[0][indexcount],list(np.where(plannp == i))[1][indexcount]]
			rnlist[indexcount] = rnnp[list(np.where(plannp == i))[0][indexcount],list(np.where(plannp == i))[1][indexcount]]
			
			try:
				x = fflist[indexcount].astype(np.float)
				y = rnlist[indexcount].astype(np.float)
			except ValueError:
				pass
			indexcount +=1
		except IndexError:
			pass
	
	legdict[i] = fflist ##  make a dictionary of ID:[fireflyvalues]
	legdictrn[i] = rnlist ## or [Renilla values]


	
## Convert Firefly (FF) or Renilla (RN) dictionaries into panda dataframes for easier calculations
legdff = pd.DataFrame.transpose(pd.DataFrame.from_dict(legdict))
legdff = legdff.apply(lambda x: x.str.strip() if isinstance(x, str) else x).replace('', np.nan).replace("-", np.nan)
legdff = legdff.apply(pd.to_numeric,errors="ignore")
legdff.where(legdff<0,0) # if negatives values are found, convert to zero. 

legdfrnf = pd.DataFrame.transpose(pd.DataFrame.from_dict(legdictrn))
legdfrnf = legdfrnf.apply(lambda x: x.str.strip() if isinstance(x, str) else x).replace('', np.nan).replace("-", np.nan)
legdfrnf = legdfrnf.apply(pd.to_numeric,errors="ignore")
legdfrnf.where(legdfrnf<0,0)

headfflist = legdff.columns.values.tolist() # extract column titles
ratiodff = legdff / legdfrnf  # calculate Relative Luciferase Units (RLU) , by normalising FF on RN values for each well. 


if dogrubbs == True: ## if the Grubbs flag was activated, identify potential outliers and remove them from the RLU values
	ratiodff = grubbsfunc(ratiodff,headfflist)



ratiomask = pd.isnull(ratiodff) # inverted mask to remove corresponding removed values from ratiodff from FF and RN dataframes as well
legdff = legdff.where(~ratiomask, other=np.nan)
legdfrnf = legdfrnf.where(~ratiomask, other=np.nan)

#Compute Mean and Average of all 3 metrics
legdff["Firefly_Average"] = legdff.mean(axis=1)
legdff["Firefly_St.Dev"] = legdff.std(axis=1)
legdfrnf["Renilla_Average"] = legdfrnf.mean(axis=1)
legdfrnf["Renilla_St.Dev"] = legdfrnf.std(axis=1)
ratiodff["RLU_Average"] = ratiodff.mean(axis=1)
ratiodff["RLU_St.Dev"] = ratiodff.std(axis=1)



#Merge all into one single DataFrame
legdf = pd.concat([legend,legdff,legdfrnf,ratiodff], axis=1)


#Next step is to calculate the fold changes. In order to do so, we must attribute to each row its associated control and its RLU value. 
controldict = {}
for i in list(legdf["ControlID"].unique()):
    controldict[i] = legdf["RLU_Average"].loc[int(i)]
legdf["Control_RLU"] = legdf["ControlID"].map(controldict)

#Calculate Fold changes versus control
fcdf = ratiodff.iloc[:,:3].div(legdf["Control_RLU"], axis=0)
fcdf["FC_Average"] = fcdf.mean(axis=1)
fcdf["FC_St.Dev"] = fcdf.std(axis=1)

#Merge Fold change dataframe and delete temporary Control_RLU column
legdf = pd.concat([legdf, fcdf], axis =1)
del legdf["Control_RLU"]



# Organize column headers
ffheader = maxcount * ["Firefly_Repl"] + ["Firefly_Average"] + ["Firefly_St.Dev"]
rnheader = maxcount * ["Renilla_Repl"] + ["Renilla_Average"] + ["Renilla_St.Dev"]
ratioheader = maxcount * ["RLU"] + ["RLU_Average"] + ["RLU_St.Dev"]
fcheader = maxcount * ["Fold_Change"] + ["FC_Average"] + ["FC_St.Dev"]
header = ["Condition", "TFBS", "Promoter", "Exp_vs_Control", "Reference_Control"] + ffheader + rnheader + ratioheader + fcheader
legdf.columns = header
legdf["Date"] = expdate
legdf["Cell line"] = cellline



#Export dataframe to file for text report. 
legdf.to_csv(outputfile,sep="\t")


if pdf == True:
	pdfreport = outputfile.replace(".txt", ".pdf")
	print "MAKING PDF..."
	import matplotlib
	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	from matplotlib.backends.backend_pdf import PdfPages
	sup = pdfreport  +"  " + cellline


	pdfopen = True
	while pdfopen == True:  # export pdf does not work if the target PDF is currently open (from a previous analysis). This prevents the script from crashing, asks the user to clsoe it and then retries
		try:
			pp = PdfPages(pdfreport)
			pdfopen = False
		
		except IOError:
			print "%s file is currently open. Close it and press any key" % pdfreport
			inputvalue = raw_input("")
			pdfopen = True


	pp = PdfPages(pdfreport)
	fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,10))
	# fig.set_canvas(plt.gcf().canvas)
	plt.style.use('seaborn-paper')
	plt.rcParams['errorbar.capsize']=3

	# plt.rcParams.update({'font.size': 5})
	# plt.rc('font', size=5)          # controls default text sizes
	# plt.rc('axes', titlesize=7)     # fontsize of the axes title
	# plt.rc('axes', labelsize=7)    # fontsize of the x and y labels
	# plt.rc('xtick', labelsize=5)    # fontsize of the tick labels
	# plt.rc('ytick', labelsize=5)    # fontsize of the tick labels
	# plt.rc('legend', fontsize=5)    # legend fontsize
	# plt.rc('figure', titlesize=10)  # fontsize of the figure title


	plt.suptitle(sup)

	ax1 = axes[0,0]
	ax2 = axes[1,0]
	ax3 = axes[0,1]
	ax4 = axes[1,1]

	proms = list(set(legdf["Promoter"]))

	promlist = []
	for i in proms:
		promlist.append(legdf[legdf["Promoter"] == i])
		

	legdf.plot(kind="bar", x="Condition", y="Firefly_Average" , yerr= "Firefly_St.Dev",ax=ax1,legend=None, color= "#fff68f", ylim=0) 

	ax1.set_xlabel("Condition")
	ax1.set_ylabel("Firefly")
	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)
	ax1.yaxis.set_ticks_position('left')
	ax1.xaxis.set_ticks_position('none')


	legdf.plot(kind="bar", x="Condition", y="Renilla_Average" , yerr="Renilla_St.Dev",ax=ax2,legend=None, color= "#ffed7c", ylim=0)  
	ax2.set_xlabel("Condition")
	ax2.set_ylabel("Renilla")
	ax2.ticklabel_format(style="sci", axis="y")
	ax2.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('none')


	legdf.plot(kind="bar", x="Condition", y="RLU_Average" , yerr="RLU_St.Dev",ax=ax3,legend=None, color= "#00ff00", ylim=0)
	ax3.set_xlabel("Condition")
	ax3.set_ylabel("RLU")
	ax3.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.yaxis.set_ticks_position('left')
	ax3.xaxis.set_ticks_position('none')


	legdf.plot(kind="bar", x="Condition", y="FC_Average" , yerr="FC_St.Dev",ax=ax4,legend=None,color= "#00ffff", ylim=0 )
	ax4.set_xlabel("Condition")
	ax4.set_ylabel("Fold Change")
	ax4.spines['right'].set_visible(False)
	ax4.spines['top'].set_visible(False)
	ax4.yaxis.set_ticks_position('left')
	ax4.xaxis.set_ticks_position('none')


	


	plt.tight_layout()
	fig.subplots_adjust(top=0.95)
	
	
	pp.savefig()
	pp.close()
			
	

	print "DONE"


