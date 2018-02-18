from Grid import Grid
from Molecule import Molecule
from BindPocFinder import BindPocFinder
import myTools

class PocketStream:
  """ Perform pocket analysis on a sequence of molecules."""

  def __init__(self,pdbfile,outputDir,includeH,minNeighbours=9,dobThreshold=9,minClust=50):

    self.pdbfile = pdbfile
    self.outputDir = outputDir
    self.includeH=includeH
    self.minNeighbours=minNeighbours
    self.dobThreshold=dobThreshold
    self.minClust=minClust

    # Set grid width and margin parameters
    self.gridmargin = 0.8
    self.gridwidth = 0.8

    self.grid = None
    self.poEx_out = [('%15s%15s%15s%15s%15s%20s\n' % (
    'System', 'Pocket id', 'Pocket volume', 'Total volume', 'BDS pocket', 'Tot. BDS volume'))]

    self.pdb_name=self.pdbfile[:-4]
    self.directory = myTools.createSubDirect (self.outputDir, self.pdb_name)
    self.mol = Molecule(self.pdbfile, self.pdb_name, self.includeH)
    self.bindpocfinder = BindPocFinder (self.mol, self.pdb_name, None, MINNEIGH=minNeighbours,DOB=dobThreshold, MINCLUST=minClust, grd=self.grid)
    self.bindpocfinder.findpocket()

    if self.bindpocfinder.clusters:
      fmt = self.directory + '/' + self.pdb_name + '_pocket%i.pdb'
      for ncluster, cluster in enumerate (self.bindpocfinder.clusters):
        fname = fmt % (ncluster + 1)
        self.write_coord_to_pdb (fname, cluster)

  def write_coord_to_pdb(self,name, coords):
    file = open(name, 'w')
    atom_serial_number = 0
    for ind in coords:
      atom_serial_number = atom_serial_number + 1
      alt_loc = ' '
      res_name = 'COO'
      chain_id = ' '
      res_num = 0
      code = ' '
      occupancy = 0.0
      temp_factor = 0.0
      seg_id = 'XXXX'
      ele_symbol = 'XX'
      charge = 'XX'
      atom_name = 'XX'
      print >> file, ('ATOM  %5i %4s%s%3s %s%4i%s   %8.3f%8.3f%8.3f%6.2f%6.2f    %4s%2s%2s') % \
        (atom_serial_number, atom_name, alt_loc, res_name,chain_id,res_num,code,ind[0],ind[1],ind[2], \
        occupancy, temp_factor, seg_id, ele_symbol, charge)
    file.close()
    self.mol.get_res_name_ident(coords)