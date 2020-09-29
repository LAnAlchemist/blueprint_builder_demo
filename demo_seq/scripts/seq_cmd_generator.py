#!/usr/bin/env python3

import os, sys, glob
'''
- written by Linna An, IPD, UW 2020/08/27 for 'In Silico Design of de novo Protein Using Blueprints Pipeline'
- Function: this script is used to generate sequence design commands for selected backbones.
- Usage: 
	>>python ./scripts/seq_cmd_generator.py [path_to_your_rosetta] [path_to_your_xml] [path_to_your_resfile] [path_to_your_input_backbone_folder]
'''

def GenSeqDesignCmd(ROSETTA,xml,resfile,pdb_dir):

	cmd_list = []

	for pdb in glob.glob(pdb_dir+'/*.pdb'):
		cmd = '{} -s {} -parser:protocol {} -parser:script_vars resfile={} cst_trp=./cst_trp -out:path:pdb ./seq_result -out:path:score ./seq_result -out:suffix _seq -nstruct 1'.\
		format(ROSETTA,pdb,xml,resfile)
		cmd_list.append(cmd)
	fp = open('./seq_design_cmd','w')
	fp.write('\n'.join(cmd_list))
	fp.close()

def main():

	if len(sys.argv) != 5:
		print ('USAGE:>>python ./scripts/seq_cmd_generator.py [path_to_your_rosetta] [path_to_your_xml] [path_to_your_resfile] [path_to_your_input_backbone_folder]')
		return

	ROSETTA = sys.argv[-4]
	xml = sys.argv[-3]
	resfile = sys.argv[-2]
	pdb_dir = sys.argv[-1]
	GenSeqDesignCmd(ROSETTA,xml,resfile,pdb_dir)

	return

if __name__ == '__main__':
	main()