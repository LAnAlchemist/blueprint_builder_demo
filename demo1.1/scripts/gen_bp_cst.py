#!/usr/bin/env python
import os
import sys
import copy
#Gyu Rie Lee, IPD, UW, 2020.
#An example python script to generate a blueprint and constraint file for a beta barrel with given topology L1E9L3E13L2E10L3E13L2E8L5E11L2E8L3E10L1
#This script is not written for general cases. It relies on pre-defined requirements for de novo beta barrel construction.
#Assignments on the positions of beta-bulges, kinks, and registry shifts were adapted from the reference;
#[Dou and Vorobieva et al., De novo design of a fluorescence-activating β-barrel, Nature, 2018.]


def bp_format(resNo, aa, ss, remodel=True):
    fmt = '%d  %s  %-2s  %s'
    remodel_stat = '.'
    if (remodel):
        resNo = 0
        remodel_stat = 'R'
    return fmt%(resNo,aa,ss,remodel_stat)

def bb_hbond_constraint(don_resNo, acc_resNo, dist_cst_parm=[3.0,0.5], ang_cst_parm=[3.1,0.3]):
    dist_cst_fmt = 'AtomPair N %d O %d HARMONIC %3.1f %3.1f'
    ang_cst_fmt_one = 'Angle N %d H %d O %d CIRCULARHARMONIC %3.1f %3.1f'
    ang_cst_fmt_two = 'Angle H %d O %d C %d CIRCULARHARMONIC %3.1f %3.1f'
    #
    cst_s = []
    cst_s.append(dist_cst_fmt%(don_resNo,acc_resNo,dist_cst_parm[0],dist_cst_parm[1]))
    cst_s.append(ang_cst_fmt_one%(don_resNo,don_resNo,acc_resNo,ang_cst_parm[0],ang_cst_parm[1]))
    cst_s.append(ang_cst_fmt_two%(don_resNo,acc_resNo,acc_resNo,ang_cst_parm[0],ang_cst_parm[1]))
    return cst_s

class BluePrintConstructor:
    def __init__(self, in_str, in_gly_kink, out_bp_fn, out_cst_fn=None):
        self.in_topo_str = in_str.strip()
        self.gly_kink_pos = in_gly_kink.strip()
        self.out_bp_fn = out_bp_fn.strip()
        self.out_cst_fn = None
        if out_cst_fn != None:
            self.out_cst_fn = out_cst_fn.strip()
        #
        self.Nterm_helix_len = 4
        #
        self.SSseg_keys = []
    def parse_SS_segment_lens(self):
        E_splits = self.in_topo_str.split('E')
        self.SSseg_to_len = {}
        for i_E,E_split in enumerate(E_splits):
            L_splits = E_split.split('L')
            L_id = i_E+1
            L_len = int(L_splits[-1])
            self.SSseg_to_len['L%d'%L_id] = L_len
            if (L_splits[0] != ''):
                E_id = copy.deepcopy(i_E)
                E_len = int(L_splits[0])
                self.SSseg_to_len['E%d'%E_id] = E_len
                last_E_id = E_id
        #
        #Order of SS segments is important for bp writing
        for i in range(1,last_E_id+1):
            self.SSseg_keys.append('L%d'%i)
            self.SSseg_keys.append('E%d'%i)
        self.SSseg_keys.append('L%d'%(i+1))
        return
    def parse_SS_segment_gly_kink_pos(self):
        E_kinks = self.gly_kink_pos.split(';')
        self.E_to_gly_kinks = {}
        for E_kink in E_kinks:
            tmp = E_kink.split('.')
            i_pos = [int(x)-1 for x in tmp[-1].split(',')]
            self.E_to_gly_kinks[tmp[0]] = i_pos
        return
    def Nterm_helix_bp(self):
        resNo = 0
        cont = []
        cont.append(bp_format(1,'V','L',remodel=False))
        resNo += 1
        cont.append(bp_format(0,'V','L'))
        resNo += 1
        for i in range(self.Nterm_helix_len):
            cont.append(bp_format(0,'V','H'))
            resNo += 1
        cont.append(bp_format(0,'V','L'))
        resNo += 1
        return cont,resNo
    def append_SS_segment_bps(self,cont=[],resNo=0):
        #Beta turn types
        L_len_to_abego = {1:' ',2:'AA',3:'AAG',4:'AAAG'}
        #
        self.SSseg_to_res_group = {}
        resNo_group = []
        for SSseg in self.SSseg_keys:
            ss_len = self.SSseg_to_len[SSseg]
            ss_state = SSseg[0]
            ss_id = SSseg[1]
            #
            tmp_seq = ['V' for i in range(ss_len)]
            if ss_state == 'L':
                abego_s = L_len_to_abego[ss_len]
            elif ss_state == 'E':
                tmp_abego = ['B' for i in range(ss_len)]
                #Even number indexed-E has beta bulge
                if (int(ss_id)%2) == 0 and int(ss_id) < 8:
                    tmp_abego[-2] = 'A'
                #Replace V with G for assigned glycine kinks
                if SSseg in self.E_to_gly_kinks.keys():
                    i_kinks = self.E_to_gly_kinks[SSseg]
                    for i_kink in i_kinks:
                        tmp_seq[i_kink] = 'G'
                        tmp_abego[i_kink] = 'E'
                abego_s = ''.join(tmp_abego)
            #
            for i_res in range(ss_len):
                resNo += 1
                resNo_group.append((resNo,abego_s[i_res]))
                ss_abego = '%s%s'%(ss_state,abego_s[i_res])
                aa = tmp_seq[i_res]
                cont.append(bp_format(0,aa,ss_abego.strip()))
            self.SSseg_to_res_group[SSseg] = resNo_group
            resNo_group = []
        #
        fout = open('%s'%self.out_bp_fn,'w')
        fout.write('%s\n'%('\n'.join(cont)))
        fout.close()
        return
    def make_blueprint_file(self):
        cont = self.Nterm_helix_bp()
        self.parse_SS_segment_lens()
        self.parse_SS_segment_gly_kink_pos()
        cont,resNo = self.Nterm_helix_bp()
        self.append_SS_segment_bps(cont,resNo)
        return
    def get_bb_hb_pairs_and_write_csts(self):
        if len(self.SSseg_keys) == 0:
            return False
        cst_s = []
        #E1-E2, E3-E4, E5-E6, E7-E8 bb hbonds
        #E2-E3, E4-E5, E6-E7 bb hbonds
        for i_E in [1,3,5,7]:
            E_odd = 'E%d'%i_E
            E_odd_res_group = self.SSseg_to_res_group[E_odd]
            E_odd_len = len(E_odd_res_group)
            E_even = 'E%d'%(i_E+1)
            E_even_res_group = self.SSseg_to_res_group[E_even]
            #
            for i_res in range(1,E_odd_len+1,2):
                i_resNo = E_odd_res_group[-i_res][0]
                j_res = i_res - 1
                j_resNo = E_even_res_group[j_res][0]
                cst_s.extend(bb_hbond_constraint(i_resNo,j_resNo))
                cst_s.extend(bb_hbond_constraint(j_resNo,i_resNo))
            #
            if i_E == 1:
                continue
            E_even = 'E%d'%(i_E-1)
            E_even_res_group = self.SSseg_to_res_group[E_even]
            #
            for i_res in range(0,E_odd_len,2):
                i_resNo = E_odd_res_group[i_res][0]
                j_res = -(i_res + 2)
                j_resNo_info = E_even_res_group[j_res]
                j_resNo = j_resNo_info[0]
                #if bulge
                if j_resNo_info[1] == 'A':
                    cst_s.extend(bb_hbond_constraint(j_resNo,i_resNo))
                    cst_s.extend(bb_hbond_constraint(j_resNo+1,i_resNo))
                    cst_s.extend(bb_hbond_constraint(i_resNo,j_resNo+1))
                else:
                    cst_s.extend(bb_hbond_constraint(j_resNo,i_resNo))
                    cst_s.extend(bb_hbond_constraint(i_resNo,j_resNo))
        #E8-E1 bb hbonds
        E_last_res_group = self.SSseg_to_res_group['E8']
        E_first_res_group = self.SSseg_to_res_group['E1']
        #
        E_first_len = len(E_first_res_group)
        for i_res in range(2,E_first_len,2):
            j_res = i_res + 1
            j_resNo = E_last_res_group[j_res][0]
            i_resNo = E_first_res_group[-i_res][0]
            cst_s.extend(bb_hbond_constraint(j_resNo,i_resNo))
            if i_res == 8:
                continue
            cst_s.extend(bb_hbond_constraint(j_resNo,j_resNo))
        #
        #Add trp corner constraints
        arg_resNo = j_resNo - 1
        cst_s.append('#For Tryptophan corner formation')
        cst_s.append('Dihedral N 8 CA 8 C 8 N 9 CIRCULARHARMONIC 2.35 0.25')
        cst_s.append('Dihedral N 8 CA 8 C 8 N 9 CIRCULARHARMONIC 2.35 0.25')
        cst_s.append('Dihedral N 7 CA 7 C 7 N 8 CIRCULARHARMONIC 5.75 0.25')
        cst_s.append('Dihedral C 6 N 7 CA 7 C 7 CIRCULARHARMONIC 4.90 0.25')
        cst_s.append('AtomPair CA 11 O 8 BOUNDED 6.9 7.5 0.5')
        cst_s.append('AtomPair CA 7 CA %d BOUNDED 8.5 10.0 0.5'%arg_resNo)
        cst_s.append('AtomPair O 6 CA %d HARMONIC 8.5 0.5'%arg_resNo)
        #
        tip_resNo_s = [self.SSseg_to_res_group['L3'][0][0],\
                       self.SSseg_to_res_group['L5'][0][0],\
                       self.SSseg_to_res_group['L7'][0][0]]
        cst_s.append('AtomPair CA 7 CA %d BOUNDED 9.5 10.5 0.5'%tip_resNo_s[1])
        cst_s.append('AtomPair CA 7 CA %d BOUNDED 7.0 8.0 0.5'%tip_resNo_s[2])
        cst_s.append('AtomPair CA 6 CA %d HARMONIC 12.0 0.5'%tip_resNo_s[0])
        cst_s.append('AtomPair CA 6 CA %d HARMONIC 9.5 0.5'%tip_resNo_s[1])
        cst_s.append('AtomPair CA 6 CA %d HARMONIC 9.5 0.5'%tip_resNo_s[2])
        #
        fout = open('%s'%self.out_cst_fn,'w')
        fout.write('%s\n'%('\n'.join(cst_s)))
        fout.close()
        return True
    def process(self):
        self.make_blueprint_file()
        status = self.get_bb_hb_pairs_and_write_csts()
        if not status:
            print ("Error: SSsegment info not set up and can't write cst file.")
        return


def main():
    if len(sys.argv) < 4:
        print ('USAGE:python gen_bp_cst.py [input_topology_string] [glycine_kink_pos_string] [output_bp_fn] [(optional)output_cst_fn]')
        return
    #
    topo_string = sys.argv[1]
    gly_kink_pos = sys.argv[2]
    out_bp_fn = sys.argv[3]
    out_cst_fn = None
    if (len(sys.argv) > 4):
        out_cst_fn = sys.argv[4]
    #
    bp_construct = BluePrintConstructor(topo_string, gly_kink_pos,\
                                        out_bp_fn, out_cst_fn)
    bp_construct.process()
    return

if __name__ == '__main__':
    main()
