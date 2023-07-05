from __future__ import print_function
import struct
import numpy as np
import math
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


data_masstable = np.load("all_ame20.npy",allow_pickle='TRUE')

k=0
for i in range(len(data_masstable)):
    notFoundFlag = False
    Qb_entry = next((item for item in data_masstable if (item["A"] == data_masstable[i]["A"] and item["Z"] == data_masstable[i]["Z"])), None)
    if (Qb_entry==None):
        notFoundFlag = True
        # print("Error Qb_entry",data_masstable[i]["A"],data_masstable[i]["Z"])
    S2n_daugter_entry = next((item for item in data_masstable if (item["A"] == data_masstable[i]["A"] and item["Z"] == data_masstable[i]["Z"]+1)), None)
    if (S2n_daugter_entry==None):
        notFoundFlag = True
        # print("Error S2n_daugter_entry",data_masstable[i]["A"],data_masstable[i]["Z"])
    is_ex_Qb2n = 0
    Qb2n = -9999.
    D_Qb2n = -9999.
    if (not notFoundFlag):
        k+=1
        # print(data_masstable[i]["A"],data_masstable[i]["Z"])
        Qb2n = Qb_entry["Qb"]-S2n_daugter_entry["S2n"]
        D_Qb2n = math.sqrt(Qb_entry["DQb"]*Qb_entry["DQb"]+S2n_daugter_entry["D_S2n"]*S2n_daugter_entry["D_S2n"])
        if (Qb_entry["is_ex_Qb"]==1 and S2n_daugter_entry["is_ex_S2n"]==1):
            is_ex_Qb2n = 1
    data_masstable[i]['Qb2n'] = Qb2n
    data_masstable[i]['D_Qb2n'] = D_Qb2n
    data_masstable[i]['is_ex_Qb2n'] = is_ex_Qb2n


np.save("all_ame20_wQb2n.npy",data_masstable)