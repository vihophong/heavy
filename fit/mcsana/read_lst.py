#!/home/daq/venv/bin/python

import time
import sys
import numpy as np
# import zmq
import random
import datetime
from array import array

from ROOT import TTree, TFile

def follow(thefile):
    #print all content
    
    start_time_info_logfile = open(sys.argv[2]+".log","w")
    flag_read = False
    ch_flag = np.zeros(10,dtype=bool)
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
            # time.sleep(0.001)
            yield line1
        if (line1 == "[DATA]"):
            flag_read = True
    print("eof")
    start_time_info_logfile.close()
    # thefile.seek(0,2)
    # while True:
    #     line = thefile.readline()
    #     if not line:
    #         time.sleep(0.1)
    #         continue
    #     yield line


if __name__ == '__main__':
    logfile = open(sys.argv[1],"r")
    nsent = 0
    loglines = follow(logfile)

    output_file = TFile.Open(sys.argv[2], 'recreate')
    # some_float = array('f', [0.])
    ch = array('i', [0])
    edge = array('i', [0])
    t = array('i', [0])
    sweeps = array('i', [0])
    tag = array('i', [0])
    data_lost = array('i', [0])
    tree = TTree('tree', 'tree')
    tree.Branch('ch', ch, 'ch/I')
    tree.Branch('edge', edge, 'edge/I')
    tree.Branch('t', t, 't/I')
    tree.Branch('sweeps', sweeps, 'sweeps/I')
    tree.Branch('tag', tag, 'tag/I')
    tree.Branch('data_lost', data_lost, 'data_lost/I')
    
    for line in loglines:
        fulldata = int(line.strip(),base=16)
        ch[0] = fulldata&0x7;
        edge[0] = (fulldata&0x8)>>3;
        t[0] = (fulldata&0xFFFFFFF0)>>4;
        sweeps[0] = (fulldata&0xFFFF00000000)>>32;
        tag[0] = (fulldata&0x7FFF000000000000)>>48;
        data_lost[0] = (fulldata&0x8000000000000000)>>63;
        # print(fulldata,ch,edge,t,sweeps,tag,data_lost)
        nsent+=1
        tree.Fill()
    tree.Write()
    output_file.Close()

