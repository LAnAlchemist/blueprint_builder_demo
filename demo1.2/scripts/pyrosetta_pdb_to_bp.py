#!/usr/bin/env python
'''
created by Gyu Rie Lee, IPD, UW, 2020.

- This is an example python script is written for report 'De novo protein design using the blueprint builder in Rosetta',
	to generate a blueprint for the supplied pdb.
- USAGE:
	>>python /path_to_script/pyrosetta_pdb_to_bp your_pdb.pdb /path_to_your_pdb/you_pdb.pdb
'''
import os
import sys
import pyrosetta as py
cwd = os.getcwd()

class BluePrintResidue:
    def __init__(self,ss,abego='',remodel=True,aa_seq='V'):
        self.ss = ss
        #
        self.abego = ''
        if (abego != ''):
            self.abego = abego
        #
        self.remodel = 'R'
        if not (remodel):
            self.remodel = '.'
        #
        self.aa_seq = 'V'
        if (aa_seq != 'V'):
            self.aa_seq = aa_seq
        return
    def write_bp_line(self,resNo_in=0):
        fmt = '%-5d%-5s%-5s%-5s'
        resNo = 0
        if (self.remodel == '.'):
            resNo = resNo_in
        ss_def = '%s%s'%(self.ss,self.abego)
        bp_info = (resNo,self.aa_seq,ss_def,self.remodel)
        return fmt%bp_info

class PoseToSSinfo:
    def __init__(self,pose_in):
        self.pose_in = pose_in
    def read_phipsi(self,resNo):
        self.tnr = self.pose_in.total_residue()
        phi = self.pose_in.phi(resNo)
        psi = self.pose_in.psi(resNo)
        return (phi,psi)
            
class PdbToBlueprint:
    def __init__(self,pdbname):
        self.pdbname = pdbname.strip()
        self.pdb_header = self.pdbname.split('/')[-1].split('.')[0]
        self.bp_fn = '%s/%s_bp'%(cwd,self.pdb_header)
        #
        self.stub_length = 1
        self.keep_seq = ['C','P','G']  #conserve C, P G residues
    def read_pose(self):
        self.pose = py.pose_from_pdb('%s'%self.pdbname)
        self.seq = self.pose.sequence()
        return
    def read_ss(self):
        dssp = py.rosetta.protocols.moves.DsspMover()
        if self.pose != None:
            dssp.apply(self.pose) #use dssp to check ss of each pose
        else:
            return False
        self.ss = self.pose.secstruct()
        return
    def read_abego(self):
        abego = py.rosetta.core.sequence.get_abego(self.pose)
        self.abego = list(abego)
        return
    def read_for_bp(self):
        self.read_pose()
        self.read_ss()
        self.read_abego()
        return
    def make_bp(self):
        self.read_for_bp()
        nres = len(self.ss)
        bp_s = []
        for i_res in range(nres):
            ss = self.ss[i_res]
            abego = self.abego[i_res]
            if (abego == '-'):
                abego = ''
            remodel = False
            if (self.seq[i_res] in self.keep_seq):
                aa_seq = self.seq[i_res]
            else:
                aa_seq = 'V'	#mutate all residues to Val
            #
            #
            resNo = i_res + 1
            if (resNo <= self.stub_length):
                remodel = False
            else:
                remodel = True
            #
            bp_line = BluePrintResidue(ss,abego,remodel,aa_seq)
            bp_s.append(bp_line.write_bp_line(resNo))
        #
        fout = open('%s'%self.bp_fn,'w')
        fout.write('%s\n'%('\n'.join(bp_s)))
        fout.close()
        return

def main():
    py.init()
    pdbname = sys.argv[1]
    pdb_bp = PdbToBlueprint(pdbname)
    pdb_bp.make_bp()
    return

if __name__ == '__main__':
    main()
