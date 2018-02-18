from Grid import Grid
import numpy
import sys


########################################################################

class BindPocFinder:
    def __init__(self, molecule, pdbcode, bdsData=None, MINNEIGH=9, DOB=9, MINCLUST=50,grd=None):
        # Input parameters
        self.mol = molecule
        #self.mol is object from Molecule class. Nothing significant in it can be removed later
        self.pdbcode = pdbcode
        self.bds = bdsData
        #Clarify about bds Data
        self.DOB = DOB
        self.MINNEIGH = MINNEIGH
        self.MINCLUST = MINCLUST
        self.grid = grd

        # Preset grid parameters
        self.gridmargin = 0.8  # These two parameters are only used for --mode="sep" -- otherwise set in PocketStream
        self.gridwidth = 0.8  # Could modify findpocket_sep to conform to other findpocket routines and delete these parameters
        self.minlength = 10

        self.bds_transf_Coords = []
        self.bds_atmnames = []
        self.bds_atmsymbols = []
        self.bds_atmkeys = []
        self.bds_resnames = []
        self.bds_resnum = []

    def findpocket(self):
        self._transform_molecule ()
        self._buildgrid ()
        print 'Grid initialized:', self.grid.initialized
        if (self.grid.initialized):
            self._get_buried_cubes ()
            self._identifycluster (self.invTransMatrix)
            self._identify_interface ()
    def _transform_molecule(self):
        """ Transform the molecule according to main axis.

        The idea is to have small coordinates for the grid, however this only rotates the molecule and doesn't
        center it at (0,0,0). We store the transformation matrix, so that we can rotate the grid points back to
        the frame of the protein.
        """
        """ Transform the molecule in the way for minimizing the grid that will be built afterwards.
        It consists of transforming the whole molecule in the same way as we have to do for 
        transforming the main axis of the molecule into the x-, y-, z- axis."""

        self.transf_Matrix = self.mol.get_transf_matrix () #Got the eigen vectors as main axis for the protein
        self.transf_Coords = self.mol.transform_molecule (self.transf_Matrix)
        self.invTransMatrix = numpy.linalg.inv (self.transf_Matrix) #Compute the (multiplicative) inverse of the transformed matrix.

    def _buildgrid(self):
        """ Build up grid """
        self.grid = Grid ()
        self.grid.initializeWithCoordinates (self.transf_Coords, self.gridmargin, self.gridwidth)
        self.grid.assignAtoms (self.mol, altCoordinates=self.transf_Coords)

    def _get_buried_cubes(self):
        self.ind_buried_cubes = self.grid.get_indices_of_buried_cubes (self.minlength, self.DOB, self.MINNEIGH)

    def _identifycluster(self, invTransMatrix):
        """ Cluster the buried cubes and transform their coordinates back to the original coordinate system """
        neighs = []
        self.clusters = []
        self.tr_clusters = []
        self.cube_clusters = []
        for cluster in self.grid.clusterBuriedCubes ():
            if len (cluster) <= self.MINCLUST: continue

            # Store coordinates
            cluster = self.grid.indexToCube (cluster)
            self.cube_clusters.append (cluster)
            coord_cluster = self.grid.cube_to_coord (cluster)
            self.tr_clusters.append (coord_cluster)
            self.clusters.append (self._transf_coord (invTransMatrix, coord_cluster))

            # Get indices of cluster points found
            ind_grid_tr_cluster = self.grid.ind_grid_vector (cluster)
            neighs.append ([len (self.grid.get_neighs (ind, ind_grid_tr_cluster)) for ind in cluster])

            #    self.tr_clusters = []
            #    self.clusters = []
            #
            #    neighs = []
            #
            #    ind_grid_bur = self.grid.ind_grid_vector(self.ind_buried_cubes)
            #
            #    # dictionary mark used to keep track of cubes already clustered
            #    mark = {}
            #    for i in self.ind_buried_cubes: mark[i] = 0
            #    for ind in self.ind_buried_cubes:
            #      if mark[ind] != 0: continue
            #      tr_cluster = []
            #      self.grid.dfs(ind, ind_grid_bur, mark, tr_cluster)
            #
            #      # we only collect cluster of a given size
            #      if len(tr_cluster)<=self.MINCLUST: continue
            #
            #      # Store coordinates
            #      coord_cluster = self.grid.cube_to_coord(tr_cluster)
            #      self.tr_clusters.append(coord_cluster)
            #      self.clusters.append(self._transf_coord(invTransMatrix, coord_cluster))
            #
            #      # Get indices of cluster points found
            #      ind_grid_tr_cluster= self.grid.ind_grid_vector(tr_cluster)
            #      neighs.append([len(self.grid.get_neighs(ind, ind_grid_tr_cluster)) for ind in tr_cluster])

        # Sort the clusters by size
        self.clusters.sort (lambda x, y: len (y) - len (x))
        self.tr_clusters.sort (lambda x, y: len (y) - len (x))
        self.cube_clusters.sort (lambda x, y: len (y) - len (x))
        neighs.sort (lambda x, y: len (y) - len (x))

        # Now get the envelope of the cluster by removing internal points
        self.envelopes = []
        for i, cl in enumerate (self.clusters):
            n = neighs[i]
            self.envelopes.append ([gp for j, gp in enumerate (cl) if n[j] < 26])

        self.tr_envelopes = []
        for i, cl in enumerate (self.tr_clusters):
            n = neighs[i]
            self.tr_envelopes.append ([gp for j, gp in enumerate (cl) if n[j] < 26])

    def _transf_coord(self, transf_mat, coords):
        transf_coords = []
        for ind in coords:
            transf_coord = [0.0, 0.0, 0.0]
            for i in (0, 1, 2):
                transf_coord[0] = transf_coord[0] + transf_mat[0][i] * ind[i]
                transf_coord[1] = transf_coord[1] + transf_mat[1][i] * ind[i]
                transf_coord[2] = transf_coord[2] + transf_mat[2][i] * ind[i]
            transf_coords.append (transf_coord)
        return transf_coords

    def _identify_interface(self):
        """ Determine the atoms lining the identified pockets. For efficiency, iterate only over the envelopes """
        self.interfaceAtoms = []

        for env in self.tr_envelopes:
            atomIDs = set ()
            ind_env = self.grid.coord_to_cube (env)
            # print 'Anzahl interface-cubes: ', len(ind_env)
            mark = numpy.zeros ((self.grid.ngrid[0], self.grid.ngrid[1], self.grid.ngrid[2]), numpy.int)
            for ind in ind_env:
                indx = ind[0]
                indy = ind[1]
                indz = ind[2]
                # Look at all neighbors
                for x in range (-3, 3):
                    newx = indx + x
                    if newx < 0: continue
                    if newx >= self.grid.ngrid[0]: continue
                    for y in range (-3, 3):
                        newy = indy + y
                        if newy < 0: continue
                        if newy >= self.grid.ngrid[1]: continue
                        for z in range (-3, 3):
                            # Avoid myself
                            if (x == 0 and y == 0 and z == 0): continue
                            newz = indz + z
                            if newz < 0: continue
                            if newz >= self.grid.ngrid[2]: continue
                            if mark[newx][newy][newz] == 1: continue  # Already treated this grid point

                            mark[newx][newy][newz] = 1
                            if not self.grid.cubes_occupied[newx][newy][newz]: continue

                            atomIDs.update ([atom["id"] for atom in self.grid.cubes_atoms[newx][newy][newz]])

            atomIDs = list (atomIDs)
            atomIDs.sort ()
            self.interfaceAtoms.append (self.mol.getAtomInformation (atomIDs=atomIDs))
            self.mol.getResidueInformation(atomIDs=atomIDs)