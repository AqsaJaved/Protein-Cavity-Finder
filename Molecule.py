import re
import numpy, random

class Molecule:
    def __init__(self, pdbfile, pdbcode=None, includeH=True):
        """ Initialise molecule with pdb filename and code (set to pdbfile if missing).

        The molecule is read from pdbfile. """
        self.initialized = False
        self.load (pdbfile, pdbcode, includeH=includeH)

    def load(self, pdbfile, pdbcode=None, includeH=True):
        """ Load molecule from PDB file pdbfile. pdbcode is set to pdbfile if missing """

        self.pdbfile = pdbfile
        if pdbcode is None:
            self.pdbcode = pdbfile
        else:
            self.pdbcode = pdbcode

        # Does not handle file-not-found exceptions: this is done up-front
        f = open (pdbfile, "r")
        lines = f.readlines ()
        f.close ()

        self.atomcoords = []
        self.atmnames = []
        self.atmsymbols = []
        self.resnames = []
        self.resnum = []
        self.atmkeys = []
        self._residueInfo = dict ()
        self._residuenames = dict ()
        self._atomInfo = []
        count = -1
        reH = re.compile ('H')
        for line in lines:
            line = line.strip ()

            # Filter for ATOM
            if not ((line[0:4] == 'ATOM')): continue
            coords = [float (line[30:38]), float (line[38:46]), float (line[46:54])]
            name = line[12:16]
            symbol = line[13:14]
            resname = line[17:20]
            resID = int (line[22:26])

            # Account for PDB files in which the element symbol is shifted from column 14
            # VMD writes such PDB files, ptraj does not
            # Fully compliant PDB files should also have the element in columns 77-78
            # Option "nowrap" to ptraj's "trajout" command may well control this behaviour
            if ((symbol != 'H') and reH.match (name)): symbol = 'H'
            if not includeH and (symbol == 'H'): continue
            count = count + 1

            self.atomcoords.append (coords)
            self.atmnames.append (name)
            self.atmsymbols.append (symbol)
            self.atmkeys.append (0)

            self.resnames.append (resname)
            self.resnum.append (resID)

            self._atomInfo.append (dict (id=count, coords=coords, name=name, symbol=symbol, residue=resID, residue_name=resname,key=0))
            if not self._residueInfo.has_key (resID):
                self._residueInfo[resID] = dict (id=resID, atomID=[], name=resname)
                self._residuenames[resID] = dict (id=resID, name=resname)
            self._residueInfo[resID]['atomID'].append (count)

        self.nAtoms = len (self.atmnames)
        self.nCoords = 3 * self.nAtoms
        self.framebytes = (self.nCoords) * 8 + (
        self.nCoords / 10 + 1)  # Numeric fields + EOL characters (in crd format)
        if self.nCoords % 10 == 0: self.framebytes -= 1  # Special case if ncoords exactly divisible by 10
        self.moltype = None
        self.initialized = True

    def numAtoms(self):
        """ Returns number of atoms in molecule """
        return self.nAtoms

    def numCoords(self):
        """ Returns number of coordinates in molecule """
        return self.nCoords

    def numResidues(self):
        """ Returns number of residues in molecule """
        resIDs = set (self.resnum)
        return len (resIDs)

    def coordinates(self):
        for atom in self._atomInfo:
            yield atom["coords"]

    def atoms(self):
        for atom in self._atomInfo:
            yield atom

    def getResidueInformation(self, resIDs=None, atomIDs=None):
        """ Return residue information for a list of residue IDs or a list of atom IDs

        Information returned is dictionary with resID, resName, and list of atomIDs
        """
        if resIDs is None:
            resIDs = set ()
        else:
            resIDs = set (resIDs)

        if atomIDs is not None:
            for i in atomIDs:
                resIDs.add (self._atomInfo[i]["residue"])
            return self.getResidueInformation (resIDs=resIDs)

        resIDs = list (resIDs)
        resIDs.sort ()
        str=''
        for res in resIDs:
            str=str+self._residueInfo[res]["name"]
        print self._shortenResidue (str)
        str = ''
        return dict ((resID, self._residueInfo[resID]) for resID in resIDs)

    def getAtomInformation(self, resIDs=None, atomIDs=None):
        """ Return atom information for a list of residue IDs or a list of atom IDs

        Information returned is dictionary with atomID, name, symbol, coords, key and resID
        """
        if atomIDs is None:
            atomIDs = set ()
        else:
            atomIDs = set (atomIDs)

        if resIDs is not None:
            for i in resIDs:
                atomIDs.update (self._residueInfo[i]['atomID'])
            return self.getAtomInformation (atomIDs=atomIDs)

        atomIDs = list (atomIDs)
        atomIDs.sort ()
        return dict ((atomID, self._atomInfo[atomID]) for atomID in atomIDs)

    def __str__(self):
        s = ["PDB file %s" % self.pdbfile]
        s.append ("Number of atoms: %d" % self.numAtoms ())
        s.append ("Number of residues: %d" % self.numResidues ())
        return "\n\t".join (s) + "\n"

    def __del__(self):
        pass

    def calc_geom_center(self):
        fX = 0.0
        fY = 0.0
        fZ = 0.0

        for coord in self.atomcoords:
            fX = fX + coord[0]
            fY = fY + coord[1]
            fZ = fZ + coord[2]
        fX = fX / self.nAtoms
        fY = fY / self.nAtoms
        fZ = fZ / self.nAtoms
        center = (fX, fY, fZ)
        return center

    def calc_main_axis(self):
        """ Calcs sym. (3,3)-matrix of main moment
        ( sum(y_i**2+z_i**2)  sum(-x_i*y_i)     sum(-x_i*z_i)     )
        ( sum(-x_i*y_i)     sum(x_i**2+z_i**2)  sum(-y_i*z_i)     )
        ( sum(-x_i*z_i)     sum(-y_i*z_i)     sum(x_i**2+y_i**2)  )
        and determines its eigen vectors.
        """
        #Clarify why the above step has been done
        c0, c1, c2 = self.calc_geom_center ()
        M = numpy.zeros ((3, 3), dtype=float)
        M = [[0] * 3, [0] * 3, [0] * 3]
        for x in self.atomcoords:
            xi = x[0] - c0
            yi = x[1] - c1
            zi = x[2] - c2
            M[0][0] = M[0][0] + xi * xi
            M[0][1] = M[0][1] + xi * yi
            M[0][2] = M[0][2] + xi * zi
            M[1][1] = M[1][1] + yi * yi
            M[1][2] = M[1][2] + yi * zi
            M[2][2] = M[2][2] + zi * zi
        M[1][0] = M[0][1]
        M[2][0] = M[0][2]
        M[2][1] = M[1][2]
        M = numpy.array (M)
        d = sum (numpy.diag (M))
        M = -M
        M[0, 0] = M[0, 0] + d
        M[1, 1] = M[1, 1] + d
        M[2, 2] = M[2, 2] + d

        eigenVals, eigenVecs = numpy.linalg.eig (M)
        eigenVecs = eigenVecs.transpose ()
        return eigenVecs

    def get_transf_matrix(self):

        # This function gets out the transformation matrix which is needed to
        # transform the main axis of the protein to the x-, y-, z- axis.We use
        # this kind of transformation when we want to put our protein within a
        # grid which has minimum volume.

        main_axis = self.calc_main_axis ()

        #    aux_mat = zeros((3,3), Float)
        #    for i in (0,1,2):
        #       for j in (0,1,2):
        #        aux_mat[i][j] = main_axis[j][i]
        #
        #    rot_mat = zeros((3,3), Float)
        #    rot_mat = inverse(aux_mat)

        #    return rot_mat
        return main_axis

    def transform_molecule(self, transf_mat):
        """ Apply transformation matrix to coordinates of molecule and return transformed coordinates.
        The molecule is unchanged. """
        #Understand how the transformation is working
        transf_molecule = []
        for ind in self.atomcoords:
            transf_coord = [0.0, 0.0, 0.0]
            for i in (0, 1, 2):
                transf_coord[0] = transf_coord[0] + transf_mat[0][i] * ind[i]
                transf_coord[1] = transf_coord[1] + transf_mat[1][i] * ind[i]
                transf_coord[2] = transf_coord[2] + transf_mat[2][i] * ind[i]
            transf_molecule.append (transf_coord)
        return transf_molecule

        # Function will be used for PPIAnalyzer by C. Pfleger
        # returns for a given coordinate set (atoms) a list of residue names and residue identifier

    def get_res_name_ident(self, list_of_coord):  # pragma: no cover
        list_of_resi = ['', '']  # residue name / residue identifier
        for i in range (len (self.atomcoords)):
            coord = self.atomcoords[i]
            if coord in list_of_coord: #problemo
                print "Bin ich hier?"
                akt_res_name = self.resnames[i]  # !# Possibly out one indent ..
                akt_res_iden = self.resnum[i]
                help = [akt_res_name, akt_res_iden]
                if help not in list_of_resi:
                    list_of_resi.extend (help)
            help = None
        return list_of_resi

    def write_transf_molecule(self, filename, transf_coord):  # pragma: no cover
        file = open (filename, 'w')
        atom_number = 0
        for ind in transf_coord:
            alt_loc = ' '
            res_name = self.resnames[atom_number]
            chain_id = ' '
            res_num = self.resnum[atom_number]
            icode = ' '
            occupancy = 0.0
            temp_factor = 0.0
            seg_id = 'XXXX'
            ele_symbol = 'XX'
            charge = 'XX'
            atom_name = self.atmnames[atom_number]
            atom_number = atom_number + 1
            print self.id
            print >> file, ('ATOM  %5i %4s%s%3s %s%4i%s   %8.3f%8.3f%8.3f%6.2f%6.2f    %4s%2s%2s') % (
            atom_number, atom_name, alt_loc, res_name, chain_id, res_num, icode, ind[0], ind[1], ind[2], occupancy,
            temp_factor, seg_id, ele_symbol, charge)
        file.close ()

    def _shortenResidue(self, str):
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        if len (str) % 3 != 0:
            raise ValueError ('Input length should be a multiple of three')

        y = ''
        for i in range (len (str) / 3):
            y += d[str[3 * i:3 * i + 3]]
        return y