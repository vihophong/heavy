"""This contains functions to manipulate reaclib v2 data file"""

import struct
import numpy as np
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

time_factor = {'s':1., 'y':31536000., 'ms': 0.001, 'd' : 86400., 'ky' : 31536000000, 'm' : 60., 'h': 3600.}

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
		return -8888
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
		return -9999
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
	Load nubase file and write it to an array of dictionaries
	
	Parameters:
	   infile ( str ): File path-name
	"""
	n_lines = sum(1 for line in open(infile))
	file1 = open(infile);
	count = 0
	nubase=[]
	nubase_stable=[]
	nubase_bminus=[]
	nubase_bplus=[]
	nubase_alpha=[]
	while True:
		count+=1
		line1 = file1.readline()
		if (not line1):
			break
		if (line1[0]=="#"):
			continue
		line1 = line1[0:len(line1)-1]
		if len(line1)<220:
			for i in range(220-len(line1)):
				line1 = line1+" "
		line1 = bytes(line1,encoding='utf8')

		(A,Zi,Ael,s_type,Mass,dMass,Exc,dExc,Orig,Isom_Unc,Isom_Inv,T12,T12_unit,dT12,Jpi,Ensdf_year,Discov_year,BR) = struct.unpack("3s1x4s3x5s1s1x13s11s12s11s2s1s1s9s2s1x7s14s2s10x4s90s12x",line1)
		(A,Zi,Ael,s_type,Mass,dMass,Exc,dExc,Orig,Isom_Unc,Isom_Inv,T12,T12_unit,dT12,Jpi,Ensdf_year,Discov_year,BR) = map(lambda x: x.decode('utf-8').strip(),(A,Zi,Ael,s_type,Mass,dMass,Exc,dExc,Orig,Isom_Unc,Isom_Inv,T12,T12_unit,dT12,Jpi,Ensdf_year,Discov_year,BR))
		
		# Process data
		A = int(A)
		Z = int(Zi[0:3])
		N = int(A)-int(Zi[0:3])
		is_gs = False
		if (Zi[3:4]=="0"):
			is_gs = True

		# Save data
		if (T12=="stbl" or T12_unit=="Zy" or T12_unit=="My" or T12_unit=="Ey" or T12_unit=="Gy" or T12_unit=="Yy" or T12_unit=="Py" or T12_unit=="Ty" or T12_unit=="My"):
			nubase_stable.append({"A":A,"Z":Z,"N":N})
		if (len(BR)>2):
			if ((BR[0:2] == "B+" or BR[0:2] == "IT") and is_gs and T12_unit!="Zy" and T12_unit!="My" and T12_unit!="Ey" and T12_unit!="Gy" and T12_unit!="Yy" and T12_unit!="Py" and T12_unit!="Ty" and T12_unit!="My"):
				#if (T12[-1]!="#"):
				nubase_bplus.append({"A":A,"Z":Z,"N":N})
		if (len(BR)>1):
			if (BR[0:1] == "A" and T12_unit!="" and dT12!="" and is_gs and T12_unit!="Zy" and T12_unit!="My" and T12_unit!="Ey" and T12_unit!="Gy" and T12_unit!="Yy" and T12_unit!="Py" and T12_unit!="Ty" and T12_unit!="My"):
				if (T12[-1]!="#"):
					nubase_alpha.append({"A":A,"Z":Z,"N":N})

		# for Beta minus
		if (len(BR)>2):
			if ((BR[0:2] == "B-" or BR[0:2] == "EC") and T12_unit!="" and dT12!="" and is_gs and T12_unit!="Zy" and T12_unit!="My" and T12_unit!="Ey" and T12_unit!="Gy" and T12_unit!="Yy" and T12_unit!="Py" and T12_unit!="Ty"):
				if (T12[-1]!="#"):
					time_f = time_factor[T12_unit]
					T12 = time_f * float(T12)
					dT12 = time_f * float(T12)
					P1n = 0.
					dP1n = 0.
					P2n = 0.
					dP2n = 0.
					if (BR.find("B-n")!=-1):
						if (BR.find("B-n ?")!=-1 or BR.find("B-n=?")!=-1 or BR.find("B-n= ?")!=-1):
							P1n = -9999.
						else:
							val = BR.split(";")
							if (val[1][0:4]=="B-n<" or val[1][0:4]=="B-n~" or val[1][0:4]=="B-n>"):
								if (val[1][0:4]=="B-n<"):
									P1n = -9999.
								if (val[1][0:4]=="B-n>"):
									P1n = -9999.
								if (val[1][0:4]=="B-n~"):
									valP1n = val[1][4:].split()
									P1n = float(valP1n[0])
									if (len(valP1n)>1):
										dP1n = float(valP1n[1])
							else:
								valP1n = val[1][4:].split()
								P1n = float(valP1n[0])
								if (len(valP1n)>1):
									dP1n = float(valP1n[1])
					if (BR.find("B-2n")!=-1):
						if (BR.find("B-2n ?")!=-1 or BR.find("B-2n=?")!=-1 or BR.find("B-2n= ?")!=-1):
							P2n = -9999.
						else:
							val = BR.split(";")
							if (val[2][0:5]=="B-2n<" or val[2][0:5]=="B-2n~" or val[2][0:5]=="B-2n>"):
								if (val[2][0:5]=="B-2n<"):
									P2n = -9999.
								if (val[2][0:5]=="B-2n>"):
									P2n = -9999.
								if (val[2][0:5]=="B-2n~"):
									valP2n = val[2][5:].split()
									P2n = float(valP2n[0])
									if (len(valP2n)>1):
										dP2n = float(valP2n[1])
							else:
								valP2n = val[2][5:].split()
								P2n = float(valP2n[0])
								if (len(valP2n)>1):
									dP2n = float(valP2n[1])
					nubase_bminus.append({"A":A,"Z":Z,"N":N, "T12":T12, "dT12":dT12, "P1n": P1n, "dP1n": dP1n, "P2n": P2n, "dP2n": dP2n})
	np.save("nubase_stable.npy",nubase_stable)
	np.save("nubase_bminus.npy",nubase_bminus)
	np.save("nubase_bplus.npy",nubase_bplus)
	np.save("nubase_alpha.npy",nubase_alpha)

load_txt('nubase_3.mas20.txt')

def drawbox(N,Z,fcolor='None',ecolor='gray', falpha = 1,linewidth=1):
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

def plot_nubase():
	magic_num = [2, 8, 20, 28, 50, 82, 126]
	for i in magic_num:
		plt.axhline(y=i+0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axhline(y=i-0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axvline(x=i+0.5,color='b',linestyle='--',linewidth=0.2)
		plt.axvline(x=i-0.5,color='b',linestyle='--',linewidth=0.2)

	nubase_stable = np.load("nubase_stable.npy",allow_pickle='TRUE')
	for i in range(len(nubase_stable)):
		plt.gca().add_patch(drawbox(nubase_stable[i]["N"],nubase_stable[i]["Z"],fcolor='k',ecolor='None',falpha = 1))

	nubase_bminus = np.load("nubase_bminus.npy",allow_pickle='TRUE')
	for i in range(len(nubase_bminus)):
		plt.gca().add_patch(drawbox(nubase_bminus[i]["N"],nubase_bminus[i]["Z"],fcolor='g',ecolor='k',falpha = 1,linewidth=0.001))
	nubase_bplus = np.load("nubase_bplus.npy",allow_pickle='TRUE')
	for i in range(len(nubase_bplus)):
		plt.gca().add_patch(drawbox(nubase_bplus[i]["N"],nubase_bplus[i]["Z"],fcolor='r',ecolor='k',falpha = 1,linewidth=0.001))
	nubase_alpha = np.load("nubase_alpha.npy",allow_pickle='TRUE')
	for i in range(len(nubase_alpha)):
		plt.gca().add_patch(drawbox(nubase_alpha[i]["N"],nubase_alpha[i]["Z"],fcolor='y',ecolor='k',falpha = 1,linewidth=0.001))
	
	plt.xlabel('Neutron number, $N$')
	plt.ylabel('Proton number, $Z$')
	plt.xlim([9.5,200])
	plt.ylim([9.5,116])
	
plot_nubase()
plt.show()
