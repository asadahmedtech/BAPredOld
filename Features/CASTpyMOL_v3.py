## PyMol plugin for the CASTp server (http://sts.bioe.uic.edu/castp)
## CASTpyMOL v3.0
## (C) 2018 Wei Tian
##
## The plugin was suggested by Sebastien Moretti
## Andrew Binkowski implemented v1.0,
## and Joe Dundas maintained v2.0.
## Wei Tian reworked on v2.0 and implemented v3.0. 

from Tkinter import *
from pymol import cmd
from tkFileDialog import *
def __init__(self):

    self.menuBar.addcascademenu('Plugin', 'MyPlugin', 'CASTp file selection',
                                label='CASTp pocket loader'
                                )
   
    self.menuBar.addmenuitem('MyPlugin', 'command',
                        'Remote PDB',
                        label='CASTp by PDB/job ID',
                        command = lambda s=self : RemotePDB(s) )
    
    self.menuBar.addmenuitem('MyPlugin', 'command',
                        'Local PDB',
                        label='CASTp from local files',
                        command = lambda s=self : LocalPDB(s) )

    self.menuBar.addmenuitem('MyPlugin', 'command', 
                        'FeedbackForm', 
                        label='Feedback', 
                        command = lambda s=self : Feedback(s) )

ver = 3.1


class Feedback:
    def __init__(self,app):
        import tkMessageBox
        tkMessageBox.showinfo('Feedback','Your feedback is important to us!\nPlease email your feedback about CASTp and/or this plugin to wtian7@uic.edu')        

#######################################################################################################
# Get pocket information from CASTp web server database.                                              #
#######################################################################################################
class RemotePDB:

    def __init__(self,app):
        if not _valid(ver,app):
            return
        
        import tkSimpleDialog
        import tkMessageBox
        import urllib
        import os

        pdburl = 'http://sts.bioe.uic.edu/castp/data/{subdatdir}/{moldir}/{pjid}.pdb'
        pocurl = 'http://sts.bioe.uic.edu/castp/data/{subdatdir}/{moldir}/tmp/{pjid}.poc'
        
        with open(urllib.urlretrieve('http://sts.bioe.uic.edu/castp/plugin/path')[0]) as f:
            pdburl = f.readline().strip()
            pocurl = f.readline().strip()
        
        pjid = tkSimpleDialog.askstring('PDB Loader','Enter the PDB or Job ID\n', parent=app.root)
        
        # user click cancel or close the dialog
        if pjid is None:
            return

        pjid = pjid.strip().lower()
        
        haserror = False
        
        # Get Pocket information from CASTp web server!  Size 4 : full structure file
        # Size 5 : a single chain structure file
        if len(pjid) in [4,5,6]:#TODO
            subdatdir = 'pdb'
            moldir = pjid[1:3]+'/'+pjid
        elif len(pjid) == 15:
            subdatdir = 'tmppdb'
            moldir = pjid 
        else:
            emessage = pjid + ' does not appear to be a PDB code'   
            haserror = True

        path = pdburl.format(subdatdir=subdatdir, moldir=moldir, pjid=pjid)
        pocpath = pocurl.format(subdatdir=subdatdir, moldir=moldir, pjid=pjid)
        
        # Try to retrieve the files if there are no previous errors.
        if not haserror:
            pdbfile = urllib.urlretrieve(path)[0]
            pocfile = urllib.urlretrieve(pocpath)[0]
            if(os.path.getsize(pdbfile) < 400 or os.path.getsize(pocfile) < 400):
                emessage = pjid + ' is not in the CASTp database'
                haserror = True
        

        if haserror:
            tkMessageBox.showerror('Sorry', emessage, parent=app.root)
            return           

        _showpoc(pjid,pdbfile,pocfile)





#######################################################################################################
# Load pocket information from files on the local machine                                             #
#######################################################################################################
class LocalPDB:
    
    def __init__(self,app):
        #if not _valid(ver):
        #    return

        import tkFileDialog
        import os
        pdbfile = tkFileDialog.askopenfilename(parent=app.root, title='Open the structure file\n'+
            'Within the same directory you must have the corresponding .poc file generated by the CASTp webserver.',
            filetypes = (("structure file","*.pdb"),("pocket file","*.poc")))
        
        # if user click canel or close the window
        if pdbfile == '':
            return 

        filedir = os.path.dirname(pdbfile)
        pdbid = pdbfile.split(os.sep)[-1]
        pdbid = pdbid[:pdbid.rfind('.')]

        pocfile = filedir + os.sep + pdbid + '.poc'
        _showpoc(pdbid, pdbfile, pocfile)




def _showpoc(pjid,pdbfile,pocfile):
    # Load the file
    cmd.load(pdbfile, pjid)
    with open(pocfile) as f:
        lines = f.readlines()
    
    pocdict = {}
    for line in lines:
        if not line.strip():
            continue
        atmid = line[6:11].strip()
        pocid = int(line[67:70].strip())
        if pocid not in pocdict:
            pocdict[pocid] = []
        pocdict[pocid].append(atmid)
    
    for pocid in sorted(pocdict.keys()):
        atomlist = pocdict[pocid]
        subsize = 20
        prefix = ''
        while len(atomlist)>0:
            sublist = atomlist[:min(subsize,len(atomlist))]
            atomlist = atomlist[min(subsize,len(atomlist)):]
            cmd.do("select {selid}, {prefix} id {sele},1,1".format(selid='Poc'+str(pocid).zfill(4), sele='+'.join(sublist), prefix=prefix))
            prefix = selid='Poc'+str(pocid).zfill(4) + ' or '


def _valid(currver,app):
    import urllib
    import tkMessageBox
    
    with open(urllib.urlretrieve('http://sts.bioe.uic.edu/castp/plugin/version')[0]) as f:
        oldestver = f.readline()
    if currver < float(oldestver):
        tkMessageBox.showerror('Sorry', 'This plugin is out-of-date\nPlease download the latest plugin from the CASTp server.', parent=app.root)
        return False
    return True
