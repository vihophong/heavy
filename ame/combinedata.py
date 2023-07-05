from __future__ import print_function
import struct
import numpy as np
import math
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


data_masstable = np.load("mass_1.mas20.npy",allow_pickle='TRUE')
data_Qbn = np.load("rct1.mas20.npy",allow_pickle='TRUE')
data_Sn = np.load("rct2_1.mas20.npy",allow_pickle='TRUE')


data_search = data_masstable
for i in range(len(data_masstable)):
	Qbn_entry = next((item for item in data_Qbn if (item["A"] == data_masstable[i]["A"] and item["Z"] == data_masstable[i]["Z"] and item["EL"] == data_masstable[i]["EL"])), None)
	if (Qbn_entry==None):
		print("Error1",data_masstable[i]["A"],data_masstable[i]["Z"])
	Sn_entry = next((item for item in data_Sn if (item["A"] == data_masstable[i]["A"] and item["Z"] == data_masstable[i]["Z"] and item["EL"] == data_masstable[i]["EL"])), None)
	if (Sn_entry==None):
		print("Error1",data_masstable[i]["A"],data_masstable[i]["Z"])
	data_masstable[i]['ZA'] = data_masstable[i]["Z"]*1000+data_masstable[i]["A"]
	
	data_masstable[i]['Sn'] = Sn_entry["Sn"]
	data_masstable[i]['D_Sn'] = Sn_entry["D_Sn"]
	data_masstable[i]['is_ex_Sn'] = Sn_entry["is_ex_Sn"]
	data_masstable[i]['Sp'] = Sn_entry["Sp"]
	data_masstable[i]['D_Sp'] = Sn_entry["D_Sp"]
	data_masstable[i]['is_ex_Sp'] = Sn_entry["is_ex_Sp"]
	data_masstable[i]['Q4b'] = Sn_entry["Q4b"]
	data_masstable[i]['D_Q4b'] = Sn_entry["D_Q4b"]
	data_masstable[i]['is_ex_Q4b'] = Sn_entry["is_ex_Q4b"]
	data_masstable[i]['Qda'] = Sn_entry["Qda"]
	data_masstable[i]['D_Qda'] = Sn_entry["D_Qda"]
	data_masstable[i]['is_ex_Qda'] = Sn_entry["is_ex_Qda"]
	data_masstable[i]['Qpa'] = Sn_entry["Qpa"]
	data_masstable[i]['D_Qda'] = Sn_entry["D_Qda"]
	data_masstable[i]['is_ex_Qda'] = Sn_entry["is_ex_Qda"]
	data_masstable[i]['Qna'] = Sn_entry["Qna"]
	data_masstable[i]['D_Qna'] = Sn_entry["D_Qna"]
	data_masstable[i]['is_ex_Qna'] = Sn_entry["is_ex_Qna"]

	data_masstable[i]['S2n'] = Qbn_entry["S2n"]
	data_masstable[i]['D_S2n'] = Qbn_entry["D_S2n"]
	data_masstable[i]['is_ex_S2n'] = Qbn_entry["is_ex_S2n"]
	data_masstable[i]['S2p'] = Qbn_entry["S2p"]
	data_masstable[i]['D_S2p'] = Qbn_entry["D_S2p"]
	data_masstable[i]['is_ex_S2p'] = Qbn_entry["is_ex_S2p"]
	data_masstable[i]['Qa'] = Qbn_entry["Qa"]
	data_masstable[i]['D_Qa'] = Qbn_entry["D_Qa"]
	data_masstable[i]['is_ex_Qa'] = Qbn_entry["is_ex_Qa"]
	data_masstable[i]['Q2b'] = Qbn_entry["Q2b"]
	data_masstable[i]['D_Q2b'] = Qbn_entry["D_Q2b"]
	data_masstable[i]['is_ex_Q2b'] = Qbn_entry["is_ex_Q2b"]
	data_masstable[i]['Qep'] = Qbn_entry["Qep"]
	data_masstable[i]['D_Qep'] = Qbn_entry["D_Qep"]
	data_masstable[i]['is_ex_Qep'] = Qbn_entry["is_ex_Qep"]
	data_masstable[i]['Qbn'] = Qbn_entry["Qbn"]
	data_masstable[i]['D_Qbn'] = Qbn_entry["D_Qbn"]
	data_masstable[i]['is_ex_Qbn'] = Qbn_entry["is_ex_Qbn"]

# Add In133 and In134 from Izzo21
	if (data_masstable[i]['EL']=='In' and data_masstable[i]['A']==133):
		data_masstable[i]['M'] = -57678.
		data_masstable[i]['DM'] = 41.
		data_masstable[i]['is_ex_M'] = 1
	if (data_masstable[i]['EL']=='In' and data_masstable[i]['A']==134):
		data_masstable[i]['M'] = -51855
		data_masstable[i]['DM'] = 44.
		data_masstable[i]['is_ex_M'] = 1

	if ((data_masstable[i]['EL']=='In' and data_masstable[i]['A']==133) or (data_masstable[i]['EL']=='In' and data_masstable[i]['A']==134)):
		S1n_entry = next((item for item in data_search if (item["Z"] == data_masstable[i]["Z"] and item["N"] == data_masstable[i]["N"]-1)), None)
		S2n_entry = next((item for item in data_search if (item["Z"] == data_masstable[i]["Z"] and item["N"] == data_masstable[i]["N"]-2)), None)
		Qbn_entry = next((item for item in data_search if (item["Z"] == data_masstable[i]["Z"]+1 and item["N"] == data_masstable[i]["N"]-2)), None)
		Qb_entry = next((item for item in data_search if (item["Z"] == data_masstable[i]["Z"]+1 and item["N"] == data_masstable[i]["N"]-1)), None)
		if (S1n_entry!=None):
			S1n_mass = -data_masstable[i]["M"] + S1n_entry["M"] + 8071.31806
			D_S1n_mass = np.sqrt(data_masstable[i]["DM"]**2+S1n_entry["DM"]**2)
		if (S2n_entry!=None):
			S2n_mass = -data_masstable[i]["M"] + S2n_entry["M"] + 8071.31806*2.
			D_S2n_mass = np.sqrt(data_masstable[i]["DM"]**2+S2n_entry["DM"]**2)
		if (Qb_entry!=None):
			Qb = data_masstable[i]["M"] - Qb_entry["M"]
			D_Qb = np.sqrt(data_masstable[i]["DM"]**2+Qb_entry["DM"]**2)
		if (Qbn_entry!=None):
			Qbn = data_masstable[i]["M"] - Qbn_entry["M"] - 8071.31806
			D_Qbn = np.sqrt(data_masstable[i]["DM"]**2+Qbn_entry["DM"]**2)
		print("Sn",S1n_mass,data_masstable[i]['Sn'])
		print("S2n",S2n_mass,data_masstable[i]['S2n'])
		print("Qb",Qb,data_masstable[i]['Qb'])
		print("Qbn",Qbn,data_masstable[i]['Qbn'])
		data_masstable[i]['Sn'] = S1n_mass
		data_masstable[i]['D_Sn'] = D_S1n_mass
		data_masstable[i]['is_ex_Sn'] = 1
		data_masstable[i]['S2n'] = S2n_mass
		data_masstable[i]['D_S2n'] = D_S2n_mass
		data_masstable[i]['is_ex_S2n'] = 1
		data_masstable[i]['Qb'] = Qb
		data_masstable[i]['DQb'] = D_Qb
		data_masstable[i]['is_ex_Qb'] = 1
		data_masstable[i]['Qbn'] = Qbn
		data_masstable[i]['D_Qbn'] = D_Qbn
		data_masstable[i]['is_ex_Qbn'] = 1
		## Sp?S2p? I'm tired...


np.save("all_ame20.npy",data_masstable)

def drawbox(N,Z,fcolor='None',ecolor='gray', falpha = 1):
	"""
	Draw box
	Parameters:
	   N ( int ): Neutron number
	   P ( int ): Proton number
	   ecolor ( str ): Color code
	"""
	rec = plt.Rectangle((N-0.5,Z-0.5),1,1,facecolor=fcolor,edgecolor=ecolor,alpha = falpha)
	#plt.text(N-0.4,Z-0.1,'$\mathregular{^{'+str(Z+N)+'}'+elements[Z]+'}$')
	return rec 

for i in range(len(data_masstable)):
	plt.gca().add_patch(drawbox(data_masstable[i]["N"],data_masstable[i]["Z"],fcolor='gray',ecolor='k',falpha = 0.5))
	if (data_masstable[i]["is_ex_M"]==1):
		plt.gca().add_patch(drawbox(data_masstable[i]["N"],data_masstable[i]["Z"],fcolor='r',ecolor='None',falpha = 0.5))
	
plt.xlabel('Neutron number, $N$')
plt.ylabel('Proton number, $Z$')
plt.xlim([0,250])
plt.ylim([0,136])
plt.show()