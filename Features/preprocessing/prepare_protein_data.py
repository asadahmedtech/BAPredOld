from Bio import *
import os

#PDB Parser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.NACCESS import run_naccess, process_rsa_data
import pickle 

#Labels for files
"""Secondary Structure in Millers Format
	H - Alpha Helix (4-12)            = 1
	B - Isolated Beta Bridge residue  = 2	
	E - Strand                        = 3
	G - 3-10 helix                    = 4
	I - Pi Helix                      = 5
	T - Turn                          = 6
	S - bend                          = 7
	- - None                          = 0
"""
SS_Labels = {'H' : 1, 'B' : 2, 'E' : 3, 'G' : 4, 'I' : 5, 'T' : 6, 'S' : 7, '-' : 0}
"""Relactive Solvent Accessiblity (RSA)
	Threshold = 25
	Exposed (> Threshold) = 1
	Burried (<= Threshold) = 0
"""
RSA_Threshold = 25

def parse_PSSM(path, file):
	pssm = {}
	with open(os.path.join(path, file), 'r') as f:
		lines = f.readlines()
		# lines = [i.split() if(len(i.split()) == 44) for i in lines]
		lines_new = []
		for i in lines:
			i = i.split()
			if(len(i) == 44):
				lines_new.append(i)
		lines = [i[:22] for i in lines]
	
	for i in lines:
		scores = i[2:]
		scores = [int(temp_i) for temp_i in scores]
		pssm[i[0]+'_'+i[1]] = scores
	return pssm
def calc_features(PATH, PDB_id, OUTPATH):

	#Loading the files
	parser = PDBParser(PERMISSIVE = 1)

	filename = os.path.join(PATH, PDB_id + ".pdb")
	structure = parser.get_structure(PDB_id, filename)
	model = structure[0]

	#DSSP Analysis for SS, PHI, PSI
	dssp = DSSP(model, filename)

	#NACCESS Analysis for SASA
	rsa, asa = run_naccess(model, filename)
	rsa = process_rsa_data(rsa)

	#Feature mapping to each atomic coordinate
	dssp_present, dssp_not_present = 0, 0
	
	feature = dict() #The feature dictionary

	for model in structure:
		for chain in model:
			for residue in chain:
				for atom in residue:
					
					print(atom.get_full_id())
					ID = (atom.get_full_id()[2], atom.get_full_id()[3])

					if(ID in list(dssp.keys())):
						if(rsa[ID]["all_atoms_abs"] > Threshold):
							rsa_label = 1
						else:
							rsa_label = 0
						feat = (SS_Labels[dssp[ID][2]], dssp[ID][4]/360, dssp[ID][5]/360, rsa_label)
						feature[tuple(atom.get_coord())] = feat

						print(ID, atom.get_coord(), feat)
						dssp_present += 1

					else:
						print("==> ID not present : ", atom.get_full_id())
						dssp_not_present += 1

	#Printing the Stats 
	print("==> STATS : PDBID : %s , DSSP PRESENT : %s , DSSP NOT PRESENT : %s"%(PDB_id, dssp_present, dssp_not_present))

	#Saving the feature to each PDB file
	with open(os.path.join(OUTPATH, PDB_id + ".dat"), "wb+") as f:
		pickle.dump(feature, f)
		print("==> Dump completed")


if __name__ == '__main__':
	input_dir = '/home/binnu/Asad/dataset/new_db/pocket_pdb/'
	output_dir = "/home/binnu/Asad/dataset/new_db/pocket_pdb_featurized/"
	files = os.listdir(input_dir)

	for file in files:
		if(file.endswith(".pdb")):
			if('2W97' in file.upper()):
				print("==> Converting file : ", file)
				calc_features(input_dir, file[:-4], output_dir)
