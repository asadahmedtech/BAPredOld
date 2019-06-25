import numpy as np
import pandas as pd
import h5py

import pybel
from extract_features import Featurizer

import os


def input_file(path):
    """Check if input file exists."""

    path = os.path.abspath(path)
    if not os.path.exists(path):
        raise IOError('File %s does not exist.' % path)
    return path


def output_file(path):
    """Check if output file can be created."""

    path = os.path.abspath(path)
    dirname = os.path.dirname(path)

    if not os.access(dirname, os.W_OK):
        raise IOError('File %s cannot be created (check your permissions).'
                      % path)
    return path


def string_bool(s):
    s = s.lower()
    if s in ['true', 't', '1', 'yes', 'y']:
        return True
    elif s in ['false', 'f', '0', 'no', 'n']:
        return False
    else:
        raise IOError('%s cannot be interpreted as a boolean' % s)

def create_features(pocket_dir, ligand_dir, ID, datafile, affinities, pocket_format = "mol2", ligand_format = "mol2"):
    ligand_files = os.listdir(ligand_dir)
    ligand_list = []
    for file in ligand_files:
        if(ID in file):
            ligand_list.append(file)
    del ligand_files


    for file in ligand_list:
        pocket = next(pybel.readfile(pocket_format, os.path.join(pocket_dir, ID.split("_")[0].lower() + "_pocket.%s"%(pocket_format))))
        ligand = next(pybel.readfile(ligand_format, os.path.join(ligand_dir, file)))
            
        pocket_coords, pocket_features = featurizer.get_features(pocket, ID[:4].lower() ,molcode=-1)
        ligand_coords, ligand_features = featurizer.get_features(ligand, ID[:4].lower() ,molcode=1)
            
        centroid = ligand_coords.mean(axis=0)
        ligand_coords -= centroid
        pocket_coords -= centroid
        data = np.concatenate(
            (np.concatenate((ligand_coords, pocket_coords)),
            np.concatenate((ligand_features, pocket_features))),
            axis=1,
        )

        dataset = datafile.create_dataset(file[:-5  ], data=data, shape=data.shape, dtype='float32', compression='lzf')
        dataset.attrs['affinity'] = affinities.loc[ID]
    print("===> File dumped : ", file)

    # for file in ligand_list:        
    #         with h5py.File(os.path.join(output_dir, file, "w")) as f:

    #             pocket = next(pybel.readfile(pocket_format, os.path.join(pocket_dir, ID.split("_")[0].lower() + "_pocket.%s"%(pocket_format))))
    #             ligand = next(pybel.readfile(ligand_format, os.path.join(ligand_dir, file))
                
    #             pocket_coords, pocket_features = featurizer.get_features(pocket, molcode=-1)
    #             ligand_coords, ligand_features = featurizer.get_features(ligand, molcode=1)
                
    #             centroid = ligand_coords.mean(axis=0)
    #             ligand_coords -= centroid
    #             pocket_coords -= centroid

    #             data = np.concatenate(
    #                 (np.concatenate((ligand_coords, pocket_coords)),
    #                 np.concatenate((ligand_features, pocket_features))),
    #                 axis=1,
    #             )

    #             dataset = f.create_dataset(name, data=data, shape=data.shape,
    #                                            dtype='float32', compression='lzf')

    #             if affinities is not None:
    #                 dataset.attrs['affinity'] = affinities.loc[ID]
    #         print("===> File dumped : ", file)

if __name__ == '__main__':
    pocket_dir = '/home/binnu/Asad/dataset/new_db/pocket_mol2/'
    ligand_dir = '/home/binnu/Asad/dataset/new_db/ligand_mol2_filtered/'
    output_dir = "/home/binnu/Asad/dataset/new_db/complex_feature/"
    protein_feature_path = '/home/binnu/Asad/dataset/new_db/protein_pdb_featurized/'

    #TODO Scrap Affinity value from web :: Status : Done
    affinities = '/home/binnu/Asad/dataset/new_db/affinity.csv'

    featurizer = Featurizer()
    charge_idx = featurizer.FEATURE_NAMES.index('partialcharge')
    global featurizer, charge_idx


    
    # files = os.listdirs(pocket_dir)

    if(affinities is not None):
        affinities = pd.read_csv(affinities)
        if 'affinity' not in affinities.columns:
            raise ValueError('There is no `affinity` column in the table')
        elif 'name' not in affinities.columns:
            raise ValueError('There is no `name` column in the table')
        affinities = affinities.set_index('name')['affinity']
    else:
        affinities = None


    # for file in files:
    #     if(file.endswith(".pdb")):
    #         print("==> Creating Feature file : ", file)
    #         calc_features(pocket_dir, ligand_dir, file[:-4], output_dir)


    pdb_ligand_ID = '/home/binnu/Asad/dataset/new_db/PDB_ligands_ID_final_NR.txt'
    with open(pdb_ligand_ID, "r") as f:
        pdb_ligand_ID = f.readlines()
    for i in range(len(pdb_ligand_ID)):
        pdb_ligand_ID[i] = pdb_ligand_ID[i][:-1]
    iterr = 1

    segmentation_fault = ['2QJY_SMA']

    file_done = []
    with open('file_done.txt', 'r') as f:
        file_done = f.readlines()
        # file_done = [i[:-1] for i in file_done]
    once = True
    for ID in pdb_ligand_ID:
        if(once):
            datafile = h5py.File(os.path.join(output_dir, 'data.hdf'), "w")
            once = False
        else:
            datafile = h5py.File(os.path.join(output_dir, 'data.hdf'), "a")

        if(ID not in segmentation_fault and ID not in file_done):
            print("==> Creating Feature file : ", ID, iterr)
            create_features(pocket_dir, ligand_dir, ID, datafile, affinities)
            file_done.append(ID)
            with open('file_done.txt', 'w') as f:
                f.writelines(file_done)
            iterr += 1

        datafile.close()

    # create_features(pocket_dir, ligand_dir, '2QM7_JDP',datafile, affinities)
    