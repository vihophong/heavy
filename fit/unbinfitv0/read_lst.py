#!/home/daq/venv/bin/python

import time
import sys
import zmq
import random
from array import array

from ROOT import TTree, TFile

def follow(thefile):
    #print all content
    flag_read = False
    while True:
        line1 = thefile.readline()
        if (not line1):
            break
        line1 = line1.strip()
        if (flag_read):
            # time.sleep(0.001)
            yield line1
        if (line1 == "[DATA]"):
            flag_read = True
    print("eof")
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

