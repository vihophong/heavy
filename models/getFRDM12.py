"""This contains functions to manipulate reaclib v2 data file"""

import struct
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


elements={"h": 1, "he": 2, "li": 3, "be": 4, "b": 5, "c": 6, "n": 7, "o": 8, "f": 9, "ne": 10, "na": 11, "mg": 12, "al": 13, 
"si": 14, "p": 15, "s": 16, "cl": 17, "ar": 18, "k": 19, "ca": 20, "sc": 21, "ti": 22, "v": 23, "cr": 24, "mn": 25, "fe": 26,
 "co": 27, "ni": 28, "cu": 29, "zn": 30, "ga": 31, "ge": 32, "as": 33, "se": 34, "br": 35, "kr": 36, "rb": 37, "sr": 38, "y": 39,
  "zr": 40, "nb": 41, "mo": 42, "tc": 43, "ru": 44, "rh": 45, "pd": 46, "ag": 47, "cd": 48, "in": 49, "sn": 50, "sb": 51, "te": 52,
   "i": 53, "xe": 54, "cs": 55, "ba": 56, "la": 57, "ce": 58, "pr": 59, "nd": 60, "pm": 61, "sm": 62, "eu": 63, "gd": 64, "tb": 65,
    "dy": 66, "ho": 67, "er": 68, "tm": 69, "yb": 70, "lu": 71, "hf": 72, "ta": 73, "w": 74, "re": 75, "os": 76, "ir": 77, "pt": 78,
     "au": 79, "hg": 80, "tl": 81, "pb": 82, "bi": 83, "po": 84, "at": 85, "rn": 86, "fr": 87, "ra": 88, "ac": 89, "th": 90, "pa": 91,
      "u": 92, "np": 93, "pu": 94, "am": 95, "cm": 96, "bk": 97, "cf": 98, "es": 99, "fm": 100, "md": 101, "no": 102, "lr": 103, "rf": 104,
       "db": 105, "sg": 106, "bh": 107, "hs": 108, "mt": 109, "ds": 110, "rg": 111, "cn": 112, "nh": 113, "fl": 114, "mc": 115, "lv": 116, "ts": 117, "og": 118,
	   "119": 119,"120": 120,"121": 121,"122": 122,"123": 123,"124": 124,"125": 125,"126": 126,"127": 127,"128": 128,"129": 129,"130": 130,
	   "131": 131,"132": 132,"133": 133,"134": 134,"135": 135,"136": 136}

Zele = []
Zele.append("n")
for key in elements:
	Zele.append(key)

def getnamebyz(z):
	"""
	Get element name by atomic number Z
	
	Parameters:
	   z ( int ): Atomic number Z
	"""
	return Zele[z]

def getZ(input):
	"""
	Get atomic number Z by element name
	
	Parameters:
	   input ( str ): Element name
	"""
	if (input==""):
		return -999999
	else:
		sep=re.split('(\d+)',input)
		if len(sep)==1:
			if sep[0]=="n":
				return int(0)
			elif (sep[0]=="p" or sep[0]=="d" or sep[0]=="t"):
				return int(1)			
			else:
				print("Something wrong! ",input)
		else:
			return int(elements[sep[0]])

def getA(input):
	"""
	Get mass number A by element name
	
	Parameters:
	   input ( str ): Element name
	"""
	if (input==""):
		return -999999
	else:
		sep=re.split('(\d+)',input)
		if len(sep)==1:
			if sep[0]=="n":
				return 1
			elif sep[0]=="p":
				return 1
			elif sep[0]=="d":
				return 2
			elif sep[0]=="t":
				return 3
			else:
				print("Something wrong! ",input)
		else:
			return int(sep[1])

def load_txt(infile):
	"""
	Load frdm table and write it to an array of dictionaries
	
	Parameters:
	   infile ( str ): File path-name
	"""
	n_lines = sum(1 for line in open(infile))
	file1 = open(infile);
	count = 0
	datafrdm=[]
	while True:
		count+=1
		line1 = file1.readline()
		if (not line1):
			break
		if (line1[0]=="#"):
			continue
		line1 = line1[0:len(line1)-1]
		line1 = bytes(line1,encoding='utf8')
		(Z,N,A,e2,e3,e4,e6,b2,b3,b4,b6,Esp,Emic,Ebind,Mth,Mexp,sexp,EFLmic,MFLth) = struct.unpack("5s5s5s10s10s10s10s10s10s10s10s10s10s10s10s10s10s10s10s",line1)
		(Z,N,A,e2,e3,e4,e6,b2,b3,b4,b6,Esp,Emic,Ebind,Mth,Mexp,sexp,EFLmic,MFLth) = map(lambda x: x.decode('utf-8').strip(),(Z,N,A,e2,e3,e4,e6,b2,b3,b4,b6,Esp,Emic,Ebind,Mth,Mexp,sexp,EFLmic,MFLth))
		(Z,N,A) = map(lambda x: int(x),(Z,N,A))
		(e2,e3,e4,e6,b2,b3,b4,b6,Esp,Emic,Ebind,Mth,Mexp,sexp,EFLmic,MFLth) = map(lambda x: float(x) if len(x)>0 == 0 else -999999,(e2,e3,e4,e6,b2,b3,b4,b6,Esp,Emic,Ebind,Mth,Mexp,sexp,EFLmic,MFLth))
		datafrdm.append({"ZA": Z*1000+N+Z,"N": N,"Z": Z,"A": N+Z,"EL": getnamebyz(Z),"Ebind":Ebind,"Mth":Mth})
	return datafrdm

Mn_MeV = 8.07131806
Mp_MeV = 7.288971064
#M = N*8.07131806+Z*7.288971064-Ebind
def getdriplines():
	datafrdm = load_txt('ADNDT-FRDM2012-TABLE.dat')
	data_bound = []
	data_bound_Qbn = []
	for i in range(len(datafrdm)):
		S1n_entry = next((item for item in datafrdm if (item["Z"] == datafrdm[i]["Z"] and item["N"] == datafrdm[i]["N"]-1)), None)
		S2n_entry = next((item for item in datafrdm if (item["Z"] == datafrdm[i]["Z"] and item["N"] == datafrdm[i]["N"]-2)), None)
		Qbn_entry = next((item for item in datafrdm if (item["Z"] == datafrdm[i]["Z"]+1 and item["N"] == datafrdm[i]["N"]-2)), None)
		Qb_entry = next((item for item in datafrdm if (item["Z"] == datafrdm[i]["Z"]+1 and item["N"] == datafrdm[i]["N"]-1)), None)
		S1n = -999999
		S1n_mass = -999999
		S2n = -999999
		Qbn = -999999
		Qb = -999999
		if (S1n_entry!=None):
			S1n = datafrdm[i]["Ebind"]-S1n_entry["Ebind"]
			S1n_mass = -datafrdm[i]["Mth"] + S1n_entry["Mth"] + 8.07131806			
		if (S2n_entry!=None):
			S2n = datafrdm[i]["Ebind"]-S2n_entry["Ebind"]
		if (Qb_entry!=None):
			Qb = datafrdm[i]["Mth"] - Qb_entry["Mth"]
		if (Qbn_entry!=None):
			Qbn = datafrdm[i]["Mth"] - Qbn_entry["Mth"] - 8.07131806
		S1p_entry = next((item for item in datafrdm if (item["Z"] == datafrdm[i]["Z"]-1 and item["N"] == datafrdm[i]["N"])), None)
		S2p_entry = next((item for item in datafrdm if (item["Z"] == datafrdm[i]["Z"]-2 and item["N"] == datafrdm[i]["N"])), None)
		S1p = -999999
		S2p = -999999
		if (S1p_entry!=None):
			S1p = datafrdm[i]["Ebind"]-S1p_entry["Ebind"]
		if (S2p_entry!=None):
			S2p = datafrdm[i]["Ebind"]-S2p_entry["Ebind"]
		isbound = False
		if (S1n>0 and S2n>0 and S1p>0 and S2p>0):
			isbound = True
		
		data_bound.append({"ZA": datafrdm[i]["ZA"],"N": datafrdm[i]["N"],"Z": datafrdm[i]["Z"],"A": datafrdm[i]["A"],"EL": datafrdm[i]["EL"],"Ebind":datafrdm[i]["Ebind"],"mass":datafrdm[i]["Mth"],"S1n":S1n,"S2n":S2n,"S1p":S1p,"S2p":S2p,"Qb":Qb,"Qbn":Qbn,"isbound":isbound})
	return data_bound
data_bound = getdriplines()
np.save("data_frdm12.npy",data_bound)

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

for i in range(len(data_bound)):
	if (data_bound[i]["isbound"]):
		plt.gca().add_patch(drawbox(data_bound[i]["N"],data_bound[i]["Z"],fcolor='gray',ecolor='None',falpha = 0.5))
	if (data_bound[i]["Z"]==50):
		print(data_bound[i]["N"],data_bound[i]["N"]*8.07131806+data_bound[i]["Z"]*7.288971064-data_bound[i]["Ebind"],data_bound[i]["mass"],data_bound[i]["Ebind"])

plt.xlabel('Neutron number, $N$')
plt.ylabel('Proton number, $Z$')
plt.xlim([0,250])
plt.ylim([0,136])
plt.show()


