'''
Read min.data & ts.data files and write .txt input file for Infomap

Daniel J. Sharpe
'''

import numpy as np

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

    def read_tsdata(self,fname,T,kappa):

        connections = []
        prefactor = kappa*Data_2_txt_infomap.k_B*T/Data_2_txt_infomap.h
        with open(fname,"r") as tsdf:
            for ts_idx, line in enumerate(tsdf.readlines()):
                line = line.split()
                ts_en = float(line[0])
                min1_idx = int(line[3])
                min2_idx = int(line[4])
                if min1_idx == min2_idx: continue # skip TSs that don't connect 2 minima
                # Arrhenius / basic TST-type law for u->v elementary transition. k / s^{-1} mol^{-1}
                ts_cost_uv = prefactor*np.exp(-(ts_en-self.min_energies[min1_idx-1])/(Data_2_txt_infomap.k_B*T))
                ts_cost_vu = prefactor*np.exp(-(ts_en-self.min_energies[min2_idx-1])/(Data_2_txt_infomap.k_B*T))
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

    mindata_fname = "min.data.regrouped"
    tsdata_fname = "ts.data.regrouped"
    txt_out_fname = "ktn.txt"
    T = 298. # temperature / K (NOT reduced units)
    kappa = 1. # dimensionless transmission coeff for TST rate const expression

    parser1 = Data_2_txt_infomap()

    parser1.read_mindata(mindata_fname)
    parser1.read_tsdata(tsdata_fname,T,kappa)
    parser1.check_for_rpt_conns()
    parser1.write_txt_out(txt_out_fname)
