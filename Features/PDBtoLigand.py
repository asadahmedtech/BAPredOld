import argparse
import os

parser = argparse.ArgumentParser(
    description='Extracts Ligand information from PDB files',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

# parser.add_argument('-f', required=True, type=str, nargs='+',
#                     help='PDB File')
# parser.add_argument('--ligand_format', type=str, default='mol2',
#                     help='file format for the ligand,'
#                          ' must be supported by openbabel')
# parser.add_argument('--output', '-o', default='./complexes.hdf',
#                     type=str,
#                     help='name for the file with the prepared structures')
# parser.add_argument('--water', '-w', default=False, type=bool,
#                     help='whether to take HOH molecules')
parser.add_argument('--verbose', '-v', default=False, type=bool,
                    help='whether to print messages')

args = parser.parse_args()

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

def ExtractDatafromPDB(path):
    if(args.verbose):
        print("==> Parsing PDB file")

    with open(path, "r") as temp_file :
        data = temp_file.readlines()
    formatted_data = [i.split() for i in data]

    hetatm, hetatmcoords, conect, hoh, hohcoords = [], [], [], [], []

    for line_number in range(len(formatted_data)):
        i = formatted_data[line_number]
        if(i[0] in ['CONECT', 'HETATM', 'HEADER']): 
            if(i[0] == 'HETATM' and i[3] != 'HOH'):
                hetatm.append(data[line_number])
                hetatmcoords.append(i[1])
            if(i[0] == 'CONECT'):
                conect.append(data[line_number])
            if(i[0] == 'HETATM' and i[3] =='HOH'):
                hoh.append(data[line_number])
                hohcoords.append(i[1])
            if(i[0] == 'HEADER'):
                header = data[line_number]

    del data, formatted_data

    data = {"HETATM" : hetatm, "HETATMCoords" : hetatmcoords, 
            "HOH" : hoh, "HOHCoords" : hohcoords,
            "CONECT" : conect, "HEADER" : header}
    return data

def createLigandPDB(path, outname):

    data = ExtractDatafromPDB(path)

    if(args.verbose):
        print("==> Creating PDB of Ligand")

    PDBout = []
    PDBout.append(data["HEADER"])
    PDBout += data["HETATM"]

    for i in data["CONECT"]:
        temp = i.split()
        temp.pop(0)
        FLAG = 0
        for j in temp : 
            if(j not in data["HETATMCoords"] or j in data["HOHCoords"]):
                FLAG = 1
        if(not FLAG):
            PDBout.append(i)

    with open(os.path.join("out/", outname), "w") as f:
        f.writelines(PDBout)
        f.flush()
    return

def convertBulk(path):
    files = os.listdir(path)
    for file in files:
        if(file.endswith(".pdb")):
            print("==> Converting ", file)
            createLigandPDB(os.path.join(path, file), file[:-4] + "_ligand.pdb")

if __name__ == "__main__":
    convertBulk("in/")