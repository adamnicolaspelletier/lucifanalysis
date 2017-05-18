

###                                            -*- Mode: Python -*-
###                                            -*- coding UTF-8 -*-
### lucifanalysis.py
### Copyright 2016 Institut de Recherche en Immunologie et Cancerologie (IRIC)
### Author :  Adam-Nicolas Pelletier
### Last modified On: 

import numpy as np
import itertools
import random
import os
import pandas as pd
from Bio import SeqIO
from StringIO import StringIO
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import sys
from outliers import smirnov_grubbs as grubbs
import scipy as sp 
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
#matplotlib.style.use('ggplot')


##Goal:  take raw luciferase readings from a text file and create an analysis report

########################################################################################################################################################
########################################################## USER INPUT ##################################################################################


########################################################################################################################################################
########################################################################################################################################################




x = int(raw_input("Input filenames in Terminal (1) or in Filenames.txt (2)?"))
dogrubbs = raw_input("Use Grubbs test to remove potential outliers? (y/n) : ")

if x == 1: ## Terminal USer input for filenames
	luciffile = raw_input("Results File (Luciferase Read):")
	


	luciflegendfile = raw_input("Legend File Name:")
	luciffile = raw_input("Dual Glo Luciferase Readings File:")
	outputfile2 = raw_input("Save Report as:")
	pdfreport = raw_input("Save PDF Report as:")
elif x ==2: # Filenames in textfile "Filenames.txt"
	filelistbare = []
	filelist = open("Filenames.txt", "rw+")
	filelistlist = filelist.readlines()
	for i in filelistlist:
		a = i
		a = a.replace("\n", "").replace("\r", "")
		filelistbare.append(a)
	luciflegendfile = filelistbare[0].split(":")[1]
	luciffile = filelistbare[1].split(":")[1]
	outputfile2 = filelistbare[2].split(":")[1]
	pdfreport= filelistbare[3].split(":")[1]
else: #Return error, abort script
	print str(x) + " is an invalid choice. Aborting."
	sys.exit()
	
filelist.close()



def luciferaseplate(readfile, legendfile):
	""" Takes a Renilla and firefly dual glo luciferase experiment and generates a report with computed values""" 
##### 
	return "Prout"


###Read legend file
legend = open(luciflegendfile, "rw+")
plan = []
plantemp = legend.readlines()
for i in plantemp: 
	a = i
	a = a.replace("\n", "").replace("\r", "")
	plan.append(a)


planlist = []
for j in plan[2:10]: #Remove perioheral wells
	a = j.replace("X", "0")


	planlist.append(a.split("\t")[1:])

plannp = np.array(planlist, dtype="int32")
plandict = {}


legtemp= plan[11:len(plan)]  
leg = []
for i in legtemp: 
	leg.append(i.split("\t"))


promoterdict = {}
ctrldict = {}
testdict = {}


for i in leg:
	try:
		if (i[3]) == "E":
			promoterdict[i[0]] = i[2]
		elif (i[3]) == "C":
			ctrldict[i[2]] = i[0]
		else:
			pass
		testdict[i[0]] = i[4]
	except IndexError:
		pass


legdict = {}
legdictrn = {}
expdict = {}
controldict = {}
try: 
	for j in leg:
		legdict[int(j[0])] = [j[1], j[2],j[3],j[4]]
		legdictrn[int(j[0])] = [j[1], j[2],j[3],j[4]]
		expdict[int(j[0])] = [j[2]]
		controldict[int(j[0])] = [j[3]]
except IndexError:
	pass

legend.close()



rawreads = open(luciffile, "rw+")
lucif = []
luciflist = []
fireflylist = []
renillalist = []
luciftemp = rawreads.readlines()
for i in luciftemp: 
	a = i
	a = a.replace("\n", "").replace("\r", "")
	lucif.append(a)
for j in lucif[5:13]:
	splitlist = j.split("\t")
	fireflylist.append(splitlist[2:14])
	renillalist.append(splitlist[17:29])


ffnp = np.array(fireflylist)
rnnp = np.array(renillalist)



unique,pos = np.unique (plannp,return_inverse=True)
counts = np.bincount(pos)

counts2 = np.delete(counts,0)
maxpos = counts2.argmax()

maxcount= counts2[maxpos]


ffheader = maxcount * ["Firefly_Repl"] + ["Firefly_Average"] + ["Firefly_St.Dev"]
rnheader = maxcount * ["Renilla_Repl"] + ["Renilla_Average"] + ["Renilla_St.Dev"]
ratioheader = maxcount * ["RLU"] + ["RLU_Average"] + ["RLU_St.Dev"]
fcheader = maxcount * ["Fold_Change"] + ["FC_Average"] + ["FC_St.Dev"]
header = ["Condition", "Promoter", "Exp_vs_Control", "Reference_Control"] + ffheader + rnheader + ratioheader + fcheader


uniquelist = list(unique)
uniquelist.pop(0)
ratiomeandict = {}

for i in uniquelist:
	
	if legdict[i][0] == "-":
		legdict[i] = ["-"] *maxcount
		legdictrn[i] = ["-"] *maxcount
	else:
		iid = []
		iid.append(i)
		ffmeanlist = []
		rnmeanlist = []
		ffstdlist = []
		rnstdlist = []
		ratiomean = []
		indexcount = 0
		fflist = maxcount* list("-")
		rnlist = maxcount * list("-")
		ratiolist = maxcount * list("-")
		foldlist = []
		
		for j in xrange(len(list(np.where(plannp == 1))[0])):
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
		
		legdict[i] = legdict[i] + fflist
		legdictrn[i] = legdictrn[i] + rnlist
	

legdf = pd.DataFrame.from_dict(legdict)
legdf = pd.DataFrame.transpose(legdf)
legdf = legdf.apply(lambda x: x.str.strip() if isinstance(x, str) else x).replace('', np.nan).replace("-", np.nan)
legdf = legdf.apply(pd.to_numeric,errors="ignore")

legdfrn = pd.DataFrame.from_dict(legdictrn)
legdfrn = pd.DataFrame.transpose(legdfrn)
legdfrn = legdfrn.apply(lambda x: x.str.strip() if isinstance(x, str) else x).replace('', np.nan).replace("-", np.nan)
legdfrn = legdfrn.apply(pd.to_numeric,errors="ignore")

headfflist = legdf.columns.values.tolist()



ffdict = {}
rndict = {}
sdlistc = [] #make a list of standard deviations for controls and for experimental values, that will be used to identify outliers with only 2 replicates. 
sdliste = []
sdlistcrn = [] #make a list of standard deviations for controls and for experimental values, that will be used to identify outliers with only 2 replicates. 
sdlistern = []

for i in xrange(len(legdf)): ###detect any outlier with either Firefly values or Renilla values, null both values if the case for that well, and average the remaining valid wells. 
	dfempty = pd.DataFrame.transpose(pd.DataFrame(np.nan, index=[i+1], columns=headfflist[4:]))
	dfempty2 = pd.DataFrame.transpose(pd.DataFrame(np.nan, index=[i+1], columns=headfflist[4:]))
	nature = legdf.iloc[int(i),headfflist[2]]
	data= legdf.iloc[int(i),headfflist[4:]]
	datarn= legdfrn.iloc[int(i),headfflist[4:]]
	data[data<0] = np.nan
	datarn[datarn<0] = np.nan

	if dogrubbs == "y":
		print "Removing outliers using Grubbs test ..."
		if len(list(data)) > 2:
			try:
				dfempty.update(grubbs.test(data, alpha=0.05))
				dfempty2.update(grubbs.test(datarn, alpha=0.05))
				z = dfempty.to_dict()
				y = dfempty2.to_dict()
				ffdict.update(z)
				rndict.update(y)
			except TypeError:
				dfempty.update(data)
				dfempty2.update(datarn)
				z = dfempty.to_dict()
				y = dfempty2.to_dict()
				ffdict.update(z)
				rndict.update(y)
	else:
		print 
		dfempty.update(data)
		dfempty2.update(datarn)
		z = dfempty.to_dict()
		y = dfempty2.to_dict()
		ffdict.update(z)
		rndict.update(y)

	if nature == "C":
		sdlistc.append(data.std())
		sdlistcrn.append(datarn.std())
	elif nature == "E":
		sdliste.append(data.std())
		sdlistern.append(datarn.std())


for i in xrange(len(legdf)): ###detect any outlier with either Firefly values or Renilla values, null both values if the case for that well, and average the remaining valid wells. 
	dfempty = pd.DataFrame.transpose(pd.DataFrame(np.nan, index=[i+1], columns=headfflist[4:]))
	dfempty2 = pd.DataFrame.transpose(pd.DataFrame(np.nan, index=[i+1], columns=headfflist[4:]))
	nature = legdf.iloc[int(i),headfflist[2]]
	data= legdf.iloc[int(i),headfflist[4:]]
	datarn= legdfrn.iloc[int(i),headfflist[4:]]
	if nature == "C":
		sdlist = sdlistc
		sdlistrn = sdlistern
	else:
		sdlist = sdliste
		sdlistrn = sdlistern


	if len(list(data)) == 2:
		if (list(data)[0]-list(data)[1]) > 2*max(sdlist):
			thresh = 1
		else:
			thresh = 0

		if thresh == 1:
			data[data >0] = np.nan
			print "Problematic replicates #%s Firefly" % i+1
		
	if len(list(data)) == 2:
		if (list(datarn)[0]-list(datarn)[1]) > 2*max(sdlistrn):
			thresh = 1
		else:
			thresh = 0

		if thresh == 1:
			datarn[datarn >0] = np.nan
			print "Problematic replicates #%s Renilla" % i+1

legdff = pd.DataFrame.from_dict(ffdict).transpose()
legdfrnf = pd.DataFrame.from_dict(rndict).transpose()
const = legdf.iloc[:,:4]

legdf = pd.concat([const,legdff], axis=1)
fmask = pd.isnull(legdf)

legdfrn = pd.concat([const,legdfrnf], axis=1)
rnmask = pd.isnull(legdfrn)

legdf = legdf.where(~rnmask, other=np.nan)
legdfrn = legdfrn.where(~fmask, other=np.nan)

legdf["Firefly_Average"] = legdff.mean(axis=1)
legdf["Firefly_St.Dev"] = legdff.std(axis=1)

legdfrn["Renilla_Average"] = legdfrnf.mean(axis=1)
legdfrn["Renilla_St.Dev"] = legdfrnf.std(axis=1)

legdf = pd.concat([legdf,legdfrn.iloc[:,4:]], axis=1)



lenlist = [x+1 for x in range(len(legdf))]



ratiodict = {}



dfempty3 = pd.DataFrame(np.nan, index=lenlist, columns=headfflist[4:])

for i in headfflist[4:]:
	#print legdff[i]
	
	dfempty3[i]= legdff[i]/legdfrnf[i] * 1000


xx = dfempty3.to_dict()
ratiodict.update(xx)


ratiodf = pd.DataFrame.from_dict(ratiodict)


ratiodf["RLU_Average"] = ratiodf.mean(axis=1)
ratiodf["RLU_St.Dev"] = ratiodf.std(axis=1)





ratiodf = pd.concat([legdf.iloc[:,0:4],ratiodf], axis=1)

legdf = pd.concat([legdf,ratiodf.iloc[:,4:]], axis=1)

controldf = ratiodf[ratiodf.iloc[:,2] == "C"]
rludict = pd.Series(controldf.RLU_Average.values, index=controldf.iloc[:,3]).to_dict()
ratiodf["RLU_CTRL"] = ratiodf.iloc[:,3].map(rludict)




fchangedf = pd.DataFrame(np.nan, index=lenlist, columns=headfflist[4:])  # set up a dataframe for fold changes
for i in headfflist[4:]:
	#print legdff[i]
	
	fchangedf[i] = (ratiodf[i] / ratiodf["RLU_CTRL"])


fchangedf["FC_Average"] = fchangedf.mean(axis=1)
fchangedf["FC_St.Dev"] = fchangedf.std(axis=1)


legdf = pd.concat([legdf,fchangedf], axis=1)
legdf.columns = header



legdf.to_csv(outputfile2,sep="\t")



pp = PdfPages(pdfreport)
fig, axes = plt.subplots(nrows=2, ncols=2)
fig.set_canvas(plt.gcf().canvas)
plt.suptitle(pdfreport)

ax1 = legdf.plot(kind="bar", x="Condition", y="Firefly_Average" , yerr="Firefly_St.Dev",ax=axes[0,0],legend=None, color= "#fff68f") 
ax1.set_xlabel("Condition")
ax1.set_ylabel("Firefly")
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('none')



ax2 = legdf.plot(kind="bar", x="Condition", y="Renilla_Average" , yerr="Renilla_St.Dev",ax=axes[1,0],legend=None, color= "#ffed7c")  
ax2.set_xlabel("Condition")
ax2.set_ylabel("Renilla")
ax2.ticklabel_format(style="sci", axis="y")
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('none')



ax3 = legdf.plot(kind="bar", x="Condition", y="RLU_Average" , yerr="RLU_St.Dev",ax=axes[0,1],legend=None, color= "#00ff00")
ax3.set_xlabel("Condition")
ax3.set_ylabel("RLU")
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('none')


ax4 = legdf.plot(kind="bar", x="Condition", y="FC_Average" , yerr="FC_St.Dev",ax=axes[1,1],legend=None,color= "#00ffff" )
ax4.set_xlabel("Condition")
ax4.set_ylabel("Fold Change")
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.yaxis.set_ticks_position('left')
ax4.xaxis.set_ticks_position('none')


matplotlib.rcParams.update({'font.size': 6})




plt.tight_layout()
fig.subplots_adjust(top=0.95)

pp.savefig()
pp.close()

# plt.show()


# fig.ff = legdf.plot(kind="bar", x="Condition", y="FF_Average" , yerr="FF_Std.Dev")
# fig.rn = legdf.plot(kind="bar", x="Condition", y="RN_Average" , yerr="RN_Std.Dev")
# fig.rlu = legdf.plot(kind="bar", x="Condition", y="RLU_Average" , yerr="RLU_Std.Dev")
# fig_fc = legdf.plot(kind="bar", x="Condition", y="FC_Average" , yerr="FC_Std.Dev")
# fig_fc.

# fig_fc.suptitle("boldfigure suptitle", fontsize=14, fontweight="bold")
# ax = fig.add_subplot(111)
# ax.set_title("Fold Change")

# ax.set_xlabel("Condition")
# ax.set_ylabel("Fold Change")

# ax.plot([2], [1], "o")

# ax.axis([0,10,0,10])



# controldf = legdf[legdf["Exp_vs_Control"] == "C"]



# rludict = pd.Series(controldf.RLU_Average.values, index=controldf.Reference_Control).to_dict()




# legdf["RLU_Ctrl"] = legdf["Reference_Control"].map(rludict)



# legdf["Fold_Change"] = (legdf["RLU_Average"] / legdf["Reference_Control"].map(rludict))




