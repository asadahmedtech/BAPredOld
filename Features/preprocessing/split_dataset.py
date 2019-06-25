from sklearn.utils import shuffle
import h5py

import os

import argparse

parser = argparse.ArgumentParser(description='Split dataset int training,'
                                 ' validation and test sets.')
parser.add_argument('--input_path', '-i', default='/home/binnu/Asad/dataset/new_db/complex_feature',
                    help='directory with pdbbind dataset')
parser.add_argument('--output_path', '-o', default='/home/binnu/Asad/dataset/new_db/complex_feature',
                    help='directory to store output files')
parser.add_argument('--size_val', '-s', type=int, default=200,
                    help='number of samples in the validation set')
parser.add_argument('--size_test', '-t', type=int, default=500,
                    help='number of samples in the validation set')
args = parser.parse_args()

# create files with the training and validation sets
with h5py.File('%s/training_set.hdf' % args.output_path, 'w') as g, \
     h5py.File('%s/validation_set.hdf' % args.output_path, 'w') as h, \
     h5py.File('%s/test_set.hdf' % args.output_path, 'w') as i:
    with h5py.File('%s/data.hdf' % args.input_path, 'r') as f:
        data_shuffled = shuffle(list(f.keys()), random_state=123)
        for pdb_id in data_shuffled[:args.size_val]:
            ds = h.create_dataset(pdb_id, data=f[pdb_id])
            ds.attrs['affinity'] = f[pdb_id].attrs['affinity']
        for pdb_id in data_shuffled[args.size_val:args.size_test+args.size_val]:
            ds = g.create_dataset(pdb_id, data=f[pdb_id])
            ds.attrs['affinity'] = f[pdb_id].attrs['affinity']
        for pdb_id in data_shuffled[args.size_test+args.size_val:]:
            ds = i.create_dataset(pdb_id, data=f[pdb_id])
            ds.attrs['affinity'] = f[pdb_id].attrs['affinity']