
import time
import sys
import numpy as np
# import zmq
import random
import datetime
from array import array

from ROOT import TTree, TFile


if __name__ == '__main__':
    #print all content
    thefile = open(sys.argv[1],"r")
    start_time_info_logfile = open(sys.argv[2]+".log","w")
    flag_read = False
    ch_flag = np.zeros(10,dtype=bool)
    with open(sys.argv[2]+".csv",'w') as csvfile:
        csvfile.write('rawdata,ch,edge,t,sweeps,tag,data_lost\n')
        while True:
            line1 = thefile.readline()
            if (not line1):
                break
            if (not flag_read):
                for i in range(10):
                    if(ch_flag[i]):
                        #REPORT-FILE from 10/31/2023 12:47:25.227  written 10/31/2023 12:47:25
                        line1spl = line1.split()
                        sec = int(line1spl[3].split(":")[2].split(".")[0])
                        sec_remainder  = "0."+line1spl[3].split(":")[2].split(".")[1]
                        epochTimeStr = datetime.datetime(int(line1spl[2].split("/")[2]),
                        int(line1spl[2].split("/")[0]),int(line1spl[2].split("/")[1]),
                        int(line1spl[3].split(":")[0]),int(line1spl[3].split(":")[1]),sec).strftime("%s")
                        epochTime = float(epochTimeStr)+float(sec_remainder)
                        line_to_write = "CHN "+str(i)+" "+line1[:-1]+" "+str(epochTime)+"\n"
                        start_time_info_logfile.writelines(line_to_write)
                        ch_flag[i]=False
                    if (line1==("[CHN%d]\n"%i)):
                        ch_flag[i]=True
            line1 = line1.strip()
            if (flag_read):
                fulldata = int(line1.strip(),base=16)
                ch = fulldata&0x7;
                edge = (fulldata&0x8)>>3;
                t = (fulldata&0xFFFFFFF0)>>4;
                sweeps = (fulldata&0xFFFF00000000)>>32;
                tag = (fulldata&0x7FFF000000000000)>>48;
                data_lost = (fulldata&0x8000000000000000)>>63;
                csvdata = "%s,%d,%d,%d,%d,%d,%d\n" % (line1.strip(),ch,edge,t,sweeps,tag,data_lost)
                csvfile.write(csvdata)
                # yield line1
            if (line1 == "[DATA]"):
                flag_read = True
        start_time_info_logfile.close()

