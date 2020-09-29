#!/usr/bin/env python
import os
import sys
import copy
'''
Gyu Rie Lee, IPD, UW, 2020
- An example script to generate resfile and xml file for sequence design on blueprint-built beta-barrel backbones, written for 'In Silico Design of de novo Protein Using Blueprints Pipeline'
- The design rules and specifications are all adapted from [Dou and Vorobieva et al., De novo design of a fluorescence-activating β-barrel, Nature, 2018.]
- Usage: >>python gen_resfile_xml_from_bp.py [input_bp_fn] [output_resfile_fn] [output_xml_fn]
'''
def resfile_format(resNo,seq_list,aa_spec='PIKAA'):
    fmt = '%d A %s %s'%(resNo,aa_spec,seq_list)
    return fmt

class BluePrintResidue:
    def __init__(self,line,resNo=0):
        self.line = line.strip()
        x = self.line.split()
        self.resNo = resNo
        self.aa = x[1]
        self.ss_abego = x[2].strip()
        self.ss_state = self.ss_abego[0]
        self.abego = ''
        if len(self.ss_abego) == 2:
            self.abego = x[2][1]
            
class BluePrintToSeqDesignPrep:
    def __init__(self,bp_fn,out_resfile_fn,out_xml_fn):
        self.bp_fn = bp_fn.strip()
        self.resfile_fn = out_resfile_fn.strip()
        self.xml_fn = out_xml_fn.strip()
    def parse_blueprint(self):
        self.ss_state_seq = []
        resNo = 0
        bp_res_s = {}
        with open('%s'%self.bp_fn) as fp:
            for line in fp:
                resNo += 1
                bp_res = BluePrintResidue(line,resNo)
                bp_res_s[bp_res.resNo] = bp_res
        tnr = copy.deepcopy(resNo)
        #
        E_id = 0
        L_id = 0
        E_block = []
        self.Eid_to_blocks = {}
        self.E_kinks = []
        L_block = []
        self.Lid_to_blocks = {}
        for resNo in range(1,tnr+1):
            curr_bp_res = bp_res_s[resNo]
            self.ss_state_seq.append(curr_bp_res.ss_state)
            #Read strands
            if curr_bp_res.ss_state == 'E':
                if curr_bp_res.abego == 'E' and curr_bp_res.aa == 'G':
                    self.E_kinks.append(resNo)
                #
                #First strand or start of new strand
                if resNo > 1 and bp_res_s[resNo-1].ss_state != 'E':
                    if E_id > 0:
                        self.Eid_to_blocks[E_id] = E_block
                    E_id += 1
                    E_block = [curr_bp_res.resNo]
                elif resNo > 1 and bp_res_s[resNo-1].ss_state == 'E':
                    E_block.append(curr_bp_res.resNo)
            #
            #Turns in between strands
            elif E_id > 0 and curr_bp_res.ss_state == 'L':
                if resNo > 1 and bp_res_s[resNo-1].ss_state != 'L':
                    if L_id > 0:
                        self.Lid_to_blocks[L_id] = L_block
                    L_id += 1
                    L_block = [curr_bp_res.resNo]
                elif resNo > 1 and bp_res_s[resNo-1].ss_state == 'L':
                    L_block.append(curr_bp_res.resNo)
        #
        self.Eid_to_blocks[E_id] = E_block
        self.Lid_to_blocks[L_id] = L_block
        return
    def gen_resfile(self):
        n_E_ids = max(self.Eid_to_blocks.keys())
        n_L_ids = max(self.Lid_to_blocks.keys())
        #
        resNo_used_and_G = []
        cont = []
        cont.append('ALLAA')
        cont.append('start')
        #Trp corner
        start_first_strand_resNo = self.Eid_to_blocks[1][0]
        end_last_strand_resNo = self.Eid_to_blocks[n_E_ids][-1]
        cont.append(resfile_format(start_first_strand_resNo-1,'P'))
        cont.append(resfile_format(start_first_strand_resNo,'G'))
        cont.append(resfile_format(start_first_strand_resNo+2,'W'))
        cont.append(resfile_format(end_last_strand_resNo-1,'R'))
        resNo_used_and_G.append(start_first_strand_resNo)
        self.trp_resNo = start_first_strand_resNo+2
        #
        #Turns
        for L_id in range(1,n_L_ids+1):
            L_block = self.Lid_to_blocks[L_id]
            len_L = len(L_block)
            #
            i = L_block[0]
            if len_L == 2:
                cont.append(resfile_format(i-1,'ST'))
                cont.append(resfile_format(i,'P'))
                cont.append(resfile_format(i+1,'DEHTY'))
            elif len_L == 3:
                cont.append(resfile_format(i-1,'ELRST'))
                cont.append(resfile_format(i,'AEKS'))
                cont.append(resfile_format(i+1,'D'))
                cont.append(resfile_format(i+2,'G'))
                resNo_used_and_G.append(i+2)
            elif len_L == 4:
                cont.append(resfile_format(i-1,'DN'))
                cont.append(resfile_format(i,'P'))
                cont.append(resfile_format(i+3,'G'))
                resNo_used_and_G.append(i+3)
        #
        #glycine kinks
        for resNo in self.E_kinks:
            if resNo in resNo_used_and_G:
                continue
            cont.append(resfile_format(resNo,'G'))
        #
        fout = open('%s'%self.resfile_fn,'w')
        fout.write('%s\n'%('\n'.join(cont)))
        fout.close()
        return
    def modify_xml_file(self):
        resfile_resNo_s = []
        with open('%s'%self.resfile_fn) as fp:
            for line in fp:
                if line.startswith('ALLAA') or line.startswith('start'):
                    continue
                x = line.strip().split()
                resfile_resNo_s.append(int(x[0]))
        resfile_resNo_s.sort()
        #
        resfile_res_str = ','.join('%d'%i for i in resfile_resNo_s)
        ss_seq_str = ''.join(self.ss_state_seq)
        #
        cont = open('././scripts/template_bb_design.xml').read()
        ch = cont.replace('{resfile_res_list}','%s'%resfile_res_str)
        ch = ch.replace('{secstruc}','%s'%ss_seq_str)
        ch = ch.replace('{TrpResNo}','%d'%self.trp_resNo)
        #
        fout = open('%s'%self.xml_fn,'w')
        fout.write(ch)
        fout.close()
        return
    def process(self):
        self.parse_blueprint()
        self.gen_resfile()
        self.modify_xml_file()
        return

def main():
    if len(sys.argv) < 3:
        print ('USAGE:python gen_resfile_xml_from_bp.py [input_bp_fn] [output_resfile_fn] [output_xml_fn]')
        return
    #
    bp_fn = sys.argv[1]
    out_resfile_fn = sys.argv[2]
    out_xml_fn = sys.argv[3]
    #
    bp_to_resfile_xml = BluePrintToSeqDesignPrep(bp_fn,out_resfile_fn,out_xml_fn)
    bp_to_resfile_xml.process()
    return

if __name__ == '__main__':
    main()

                
                
                    
                    
                
