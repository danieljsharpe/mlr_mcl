'''
Read min.data & ts.data files and write .txt input file for Infomap

Arguments:
-Calculate rate constant: (1) simple PE TST (2) read from file "ts_weights.dat"
-min.data file name (e.g. "min.data") -ts.data file name (e.g. "ts.data")
-output file name (e.g. "ktn.txt")

Daniel J. Sharpe
'''

import numpy as np
import sys

class Data_2_txt_infomap(object):

    k_B = 1.9872036E-03 # Boltzmann constant / kcal K^{-1} mol^{-1}
    h = 1.5836687E-37 # Planck's constant / kcal s

    def __init__(self):
        pass

    def read_mindata(self,fname):

        self.min_energies = []
        with open(fname,"r") as mdf:
            for min_idx, line in enumerate(mdf.readlines()):
                line = line.split()
                min_en = float(line[0])
                self.min_energies.append(min_en)

    def read_ts_weights(self):

        self.ts_weights = []
        with open("ts_weights.dat","r") as tswf:
            for line in tswf.readlines():
                self.ts_weights.append(float(line))

    def read_tsdata(self,fname,T,kappa,k_opt):

        connections = []
        if k_opt==1: prefactor = kappa*Data_2_txt_infomap.k_B*T/Data_2_txt_infomap.h
        if k_opt==2: wts_count = 0
        with open(fname,"r") as tsdf:
            for ts_idx, line in enumerate(tsdf.readlines()):
                ts_cost_uv, ts_cost_vu = 0., 0.
                line = line.split()
                ts_en = float(line[0])
                min1_idx = int(line[3])
                min2_idx = int(line[4])
                if min1_idx == min2_idx: continue # skip TSs that don't connect 2 minima
                if k_opt==1: # Arrhenius / basic TST-type law for u->v elementary transition. k / s^{-1} mol^{-1}
                    ts_cost_uv = prefactor*np.exp(-(ts_en-self.min_energies[min1_idx-1])/(Data_2_txt_infomap.k_B*T))
                    ts_cost_vu = prefactor*np.exp(-(ts_en-self.min_energies[min2_idx-1])/(Data_2_txt_infomap.k_B*T))
                elif k_opt==2: # read from file containing logs of rate constants
                    ts_cost_uv = np.exp(self.ts_weights[wts_count])
                    ts_cost_vu = np.exp(self.ts_weights[wts_count+1])
                    wts_count += 2
                else:
                    quit("not provided a valid value for k_opt")
                connections.append([min1_idx,min2_idx,ts_cost_uv])
                connections.append([min2_idx,min1_idx,ts_cost_vu]) # recall: bidirectional edges
        self.connections = sorted(connections,key=lambda x: (x[0],x[1]))

    def check_for_rpt_conns(self):

        blacklist_ts_idx = []
        for i in range(len(self.connections)-1):
           if self.connections[i][0] == self.connections[i+1][0] and \
               self.connections[i][1] == self.connections[i+1][1]:
               if self.connections[i][2] < self.connections[i+1][2]:
                   blacklist_ts_idx.append(i+1)
               else:
                   blacklist_ts_idx.append(i)
        n_ts_del = 0
        for ts_idx in blacklist_ts_idx:
            del self.connections[ts_idx-n_ts_del]
            n_ts_del += 1

    def write_txt_out(self,txt_out_fname):

        with open(txt_out_fname,"w") as txtf:
            for connection in self.connections:
                txtf.write(str(connection[0])+" "+str(connection[1])+" "+ \
                           str(connection[2])+"\n")

if __name__=="__main__":

    k_opt = int(sys.argv[1])
    mindata_fname = sys.argv[2]
    tsdata_fname = sys.argv[3]
    txt_out_fname = sys.argv[4]
    T = 298. # temperature / K (NOT reduced units)
    kappa = 1. # dimensionless transmission coeff for TST rate const expression

    parser1 = Data_2_txt_infomap()

    if k_opt==1: parser1.read_mindata(mindata_fname)
    elif k_opt==2: parser1.read_ts_weights()
    parser1.read_tsdata(tsdata_fname,T,kappa,k_opt)
    parser1.check_for_rpt_conns()
    parser1.write_txt_out(txt_out_fname)
