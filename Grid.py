import math

class GridError(Exception):
  """ Custom exception to be used in Grid class """
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)


class Grid:
  def __init__(self):
    """ Constructor of grid object. Need to initialise grid before use """
    self.initialized = False
    self.atomsAssigned = False
    self.buriedCubesDetected = False
    self.buriedCubesClustered = False
    self.OVinitialized = False

  def initializeWithCoordinates(self, coordinates, margin=None, gridspacing=None):
    """ Initialize with a list of coordinates """
    if len (coordinates) == 0:
      raise ValueError, "You must provide at least one point to initialize grid"
    min = list (coordinates[0])
    max = list (coordinates[0])
    for c in coordinates:
      if c[0] < min[0]: min[0] = c[0]
      if c[1] < min[1]: min[1] = c[1]
      if c[2] < min[2]: min[2] = c[2]
      if c[0] > max[0]: max[0] = c[0]
      if c[1] > max[1]: max[1] = c[1]
      if c[2] > max[2]: max[2] = c[2]
    width = [max[i] - min[i] for i in range (3)]
    center = [0.5 * (max[i] + min[i]) for i in range (3)]
    self.initializeWithCenter (center, width, margin=margin, gridspacing=gridspacing)

  def initializeWithCenter(self, center, width, margin=None, gridspacing=None):
    """ Initialize with center point and width of grid along coordinate axis.

    The center point will be a grid point and an equal number of grid points
    positioned in each direction. """
    self.initialized = False
    if margin is None: margin = 1.5
    if gridspacing is None: gridspacing = 1.0

    if margin < 0:
      raise ValueError, "Margin used for grid initialisation must be greater equal 0"
    if gridspacing <= 0:
      raise ValueError, "Grid spacing must be larger than 0"

    # Determine size of grid and minimum and maximum coordinates
    ngrid = [int (math.ceil ((0.5 * width[i] + margin) / gridspacing)) for i in range (3)]
    self.min = tuple (center[i] - ngrid[i] * gridspacing for i in range (3))
    self.max = tuple (center[i] + ngrid[i] * gridspacing for i in range (3))
    self.ngrid = tuple (2 * ngrid[i] + 1 for i in range (3))
    self.margin = margin
    self.gridspacing = gridspacing

    if max (self.ngrid) > 200:
      raise GridError, "Grid too large (%d, %d, %d). Increase the grid spacing. Current value %s" % (
      self.ngrid[0], self.ngrid[1], self.ngrid[2], str (self.gridspacing))

    # Build up the grid.
    # First we build the data structures for cubes and cubes_occupied
    # Cubes contains the indices of the atoms for each grid point
    # cubes_occupied is a flag indicating if a grid point is empty or not
    # cubes_atoms contains a list of the atoms for each grid point
    self.cubes = [None] * self.ngrid[0]
    self.cubes_occupied = [None] * self.ngrid[0]
    self.cubes_buried = [None] * self.ngrid[0]
    for x in range (self.ngrid[0]):
      self.cubes[x] = [None] * self.ngrid[1]
      self.cubes_occupied[x] = [None] * self.ngrid[1]
      self.cubes_buried[x] = [None] * self.ngrid[1]
      for y in range (self.ngrid[1]):
        # We cannot use [[]]*self.steps[2] here
        self.cubes[x][y] = [[] for i in range (self.ngrid[2])]
        self.cubes_occupied[x][y] = [False] * self.ngrid[2]
        self.cubes_buried[x][y] = [False] * self.ngrid[2]

    # and then we create copies for the remaining arrays
    import cPickle as pickle
    c = pickle.dumps (self.cubes)
    self.cubes_atoms = pickle.loads (c)

    self.initialized = True
    self.buriedCubesDetected = False
    self.buriedCubesClustered = False
    self.atomsAssigned = False

  def assignAtoms(self, molecule, altCoordinates=None):
      """ Reset grid and process molecule atoms """
      # Assign atoms to cubes
      # Loop over all coordinates
      # Make radius dependent on atomname (see Bondi radii)
      # Decreased all radii by 0.1A (HG 23/08/05)
      # IC - Radii appear to be vdW radii, and are consistent with those in A. Bondi, J. Phys. Chem., 1964, 68, 441.
      # IC - Added F=1.37 for fluorine but should also be OK-ish for Fe3+
      if not self.initialized: raise GridError, "Grid not initialized."
      self.resetGrid ()
      atomRadii = dict (H=1.1, C=1.6, N=1.45, O=1.42, S=1.7, P=1.7, F=1.37)
      rDefault = 1.6
      stahlTol = 0.8

      for i, atom in enumerate (molecule.atoms ()):
        atomsymbol = atom["symbol"]
        if altCoordinates is not None:
          coords = altCoordinates[i]
        else:
          coords = atom["coords"]

        radius = atomRadii.get (atomsymbol, rDefault)

        # IC - Symmetrised around coords[i] using ceil/floor
        #      Represents the atom as a group of cubes, rather than a sphere ..
        #      imin = [int((coords[i] - radius - self.min[i]) / self.gridspacing) for i in range(3)]
        #      imax = [int((coords[i] + radius - self.min[i]) / self.gridspacing) for i in range(3)]
        imin = [int (round ((coords[i] - radius - self.min[i]) / self.gridspacing)) for i in range (3)]
        imax = [int (round ((coords[i] + radius - self.min[i]) / self.gridspacing)) for i in range (3)]
        # Make sure we stay inside grid
        imin = [max (0, imin[i]) for i in range (3)]
        imax = [min (self.ngrid[i], imax[i] + 1) for i in range (3)]
        for xind in range (imin[0], imax[0]):
          for yind in range (imin[1], imax[1]):
            for zind in range (imin[2], imax[2]):
              self.cubes[xind][yind][zind].append (atom["id"])
              self.cubes_occupied[xind][yind][zind] = True
              self.cubes_atoms[xind][yind][zind].append (atom)
      self.atomsAssigned = True
      self.buriedCubesDetected = False
      self.buriedCubesClustered = False

  def get_indices_of_cubes(self):
    """ Returns a list of indices to access the grid points.
    Implemented as generator """
    if not self.initialized: raise GridError, "Grid not initialized."
    for ix in range(self.ngrid[0]):
      for iy in range(self.ngrid[1]):
        for iz in range(self.ngrid[2]):
          yield (ix,iy,iz)

  def resetGrid(self):
    """ Reset grid data structures to reflect an empty grid """
    for x,y,z in self.get_indices_of_cubes(): self.cubes[x][y][z] = []
    for x,y,z in self.get_indices_of_cubes(): self.cubes_occupied[x][y][z] = False
    for x,y,z in self.get_indices_of_cubes(): self.cubes_buried[x][y][z] = False
    for x,y,z in self.get_indices_of_cubes(): self.cubes_atoms[x][y][z] = []
    self.buriedCubesDetected = False
    self.buriedCubesClustered = False
    self.atomsAssigned = False

  def get_indices_of_buried_cubes(self, minLength=10, dobThreshold=9, minNeighbours=9):
    """ Returns a list of indices to access buried grid points. """
    if not self.initialized: raise GridError, "Grid not initialized."
    if not self.atomsAssigned: raise GridError, "No atoms assigned to Grid"
    if not self.buriedCubesDetected:
      self.detectBuriedCubes(minLength,dobThreshold,minNeighbours)
    return [(ix,iy,iz) for ix,iy,iz in self.get_indices_of_cubes()
                           if self.cubes_buried[ix][iy][iz]]

  def detectBuriedCubes(self, minLength=10, dobThreshold=9, minNeighbours=9):
    """ First, identifies indices of empty cubes that have more than DOB non-empty neighbours.
        The resulting list of cubes from this first step is reduced by ensuring that each of these cubes
        has at least MINNEIGH neighbours.
    """
    # Reset information about empty cubes
    for x,y,z in self.get_indices_of_cubes(): self.cubes_buried[x][y][z] = False

    buried_cubes = []
    for x,y,z in self.get_indices_of_empty_cubes():
      buriedness = self.count_buriedness(minLength, x,y,z)
      if buriedness > dobThreshold:
        buried_cubes.append( (x,y,z) )

    ind_buried_cubes = []
    ind_grid_bur = self.ind_grid_vector(buried_cubes)
    for ind in buried_cubes:
      neighs = self.get_neighs(ind,ind_grid_bur)
      if (len(neighs) > minNeighbours):
        self.cubes_buried[ind[0]][ind[1]][ind[2]] = True
        ind_buried_cubes.append(ind)

    self.buriedCubesDetected = True
    self.buriedCubesClustered = False
    return ind_buried_cubes

  def get_indices_of_empty_cubes(self):
    """ Returns a list of indices to access empty grid points.
    Implemented as generator """
    if not self.initialized: raise GridError, "Grid not initialized."
    for ix,iy,iz in self.get_indices_of_cubes():
      if not self.cubes_occupied[ix][iy][iz]: yield (ix,iy,iz)

  def get_neighs(self, ind_of_cube,ind_grid_vector):
    """ Get list of neighbours for an index (ind_of_cube) """
    neighs = []
    indx = ind_of_cube[0]
    indy = ind_of_cube[1]
    indz = ind_of_cube[2]
    # Look at all neighbors avoiding out of bounds
    for x in (-1, 0, 1):
      newx = indx + x
      if(newx < 0 or newx >= self.ngrid[0]): continue
      for y in (-1, 0, 1):
        newy = indy + y
        if(newy < 0 or newy >= self.ngrid[1]): continue
        for z in (-1, 0, 1):
          # Avoid myself
          if(x == 0 and y == 0 and z == 0): continue
          newz = indz + z
          if(newz < 0 or newz >= self.ngrid[2]): continue
          # Look if neighbor is in ind_of_cubes
          index = newx*self.ngrid[2]*self.ngrid[1] + newy*self.ngrid[2] + newz
          if ind_grid_vector[index] != -1:
            neighs.append( (newx, newy, newz) )
    return neighs

  def ind_grid_vector(self, ind_cubes):
    """ This function maps a list of indices onto a vector that contains -1 or the index position for indices
        in the list. """
    vector_index = [-1]*(self.ngrid[0]*self.ngrid[1]*self.ngrid[2])
    for ind in ind_cubes:
      index = ind[0]*self.ngrid[2]*self.ngrid[1] + ind[1]*self.ngrid[2] + ind[2]
      vector_index[index] = index
    return vector_index

  def count_buriedness(self, minLength, xind, yind, zind):
    """ Determines number of occupied grid points along the grid axis and the diagonal

    For a grid point, directions checked are
      (+00)(-00) (0+0)(0-0) (00+)(00-) six along grid axis
      (+++)(---) (++-)(--+) (+--)(-++) (+-+)(-+-) eight along diagonals

    Maximum values are 4 for a corner grid point, 6 for a point on the grid edge, and 10 for
    a grid point inside.
    IC - Isn't 14 the maximum for a grid point inside?
    minLength determines the number of grid points to check away from current
    grid point (xind, yind, zind).

    PG - This code is optimized for speed. Test correctness of code by comparing with
    count_buriedness2!
    """
    minSteps = minLength / self.gridspacing
    maxSteps = int (minSteps) + 2
    px = xind + maxSteps
    if px > self.ngrid[0]: px = self.ngrid[0]
    py = yind + maxSteps
    if py > self.ngrid[1]: py = self.ngrid[1]
    pz = zind + maxSteps
    if pz > self.ngrid[2]: pz = self.ngrid[2]
    mx = xind - maxSteps
    if mx < -1: mx = -1
    my = yind - maxSteps
    if my < -1: my = -1
    mz = zind - maxSteps
    if mz < -1: mz = -1

    cntMain = 0
    # X-coords (up & down,side directions)
    for n in xrange (xind + 1, px):
      if self.cubes_occupied[n][yind][zind]:
        cntMain = cntMain + 1
        break
    # print "+00",cntMain
    for n in xrange (xind - 1, mx, -1):
      if self.cubes_occupied[n][yind][zind]:
        cntMain = cntMain + 1
        break
    # print "-00",cntMain

    # Y-coords (up & down, side directions)
    for n in xrange (yind + 1, py):
      if self.cubes_occupied[xind][n][zind]:
        cntMain = cntMain + 1
        break
    # print "0+0",cntMain
    for n in xrange (yind - 1, my, -1):
      if self.cubes_occupied[xind][n][zind]:
        cntMain = cntMain + 1
        break
    # print "0-0",cntMain

    # Z-coords (up & down, side directions)
    for n in xrange (zind + 1, pz):
      if self.cubes_occupied[xind][yind][n]:
        cntMain = cntMain + 1
        break
    # print "00+",cntMain
    for n in xrange (zind - 1, mz, -1):
      if self.cubes_occupied[xind][yind][n]:
        cntMain = cntMain + 1
        break
    # print "00-",cntMain

    ## for all 8 possible diagonals

    # (+1,+1,+1)
    newx = xind + 1
    newy = yind + 1
    newz = zind + 1
    while (newx < px) and (newy < py) and (newz < pz):
      if self.cubes_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx + 1
      newy = newy + 1
      newz = newz + 1
    # print "+++",cntMain

    # (-1,-1,-1)
    newx = xind - 1
    newy = yind - 1
    newz = zind - 1
    while (newx > mx) and (newy > my) and (newz > mz):
      if self.cubes_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx - 1
      newy = newy - 1
      newz = newz - 1
    # print "---",cntMain

    # (+1,+1,-1)
    newx = xind + 1
    newy = yind + 1
    newz = zind - 1
    while (newx < px) and (newy < py) and (newz > mz):
      if self.cubes_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx + 1
      newy = newy + 1
      newz = newz - 1

    # (-1,-1,+1)
    newx = xind - 1
    newy = yind - 1
    newz = zind + 1
    while (newz < pz) and (newx > mx) and (newy > my):
      if self.cubes_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx - 1
      newy = newy - 1
      newz = newz + 1

    # (+1,-1,-1)
    newx = xind + 1
    newy = yind - 1
    newz = zind - 1
    while (newx < px) and (newy > my) and (newz > mz):
      if self.cubes_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx + 1
      newy = newy - 1
      newz = newz - 1
    # print "+--",cntMain

    # (-1,+1,+1)
    newx = xind - 1
    newy = yind + 1
    newz = zind + 1
    while (newx > mx) and (newy < py) and (newz < pz):
      if self.cubes_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx - 1
      newy = newy + 1
      newz = newz + 1
    # print "-++",cntMain

    # (+1,-1,+1)
    newx = xind + 1
    newy = yind - 1
    newz = zind + 1
    while (newx < px) and (newy > my) and (newz < pz):
      if self.cubes_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx + 1
      newy = newy - 1
      newz = newz + 1
    # print "+-+",cntMain

    # (-1,+1,-1)
    newx = xind - 1
    newy = yind + 1
    newz = zind - 1
    while (newx > mx) and (newy < py) and (newz > mz):
      if self.cubes_occupied[newx][newy][newz]:
        cntMain = cntMain + 1
        break
      newx = newx - 1
      newy = newy + 1
      newz = newz - 1
    # print "-+-",cntMain
    return cntMain

  def clusterBuriedCubes(self):
    if not self.initialized: raise GridError, "Grid not initialized."
    if not self.atomsAssigned: raise GridError, "No atoms assigned to Grid"
    if not self.buriedCubesDetected: raise GridError, "Buried cubes not detected. Call detectBuriedCubes first."

    indBuriedCubes = self.get_indices_of_buried_cubes ()
    clusters = []
    mark = set (i for i in indBuriedCubes)
    ind_grid_bur = self.ind_grid_vector (indBuriedCubes)
    for ind in indBuriedCubes:
      if not ind in mark: continue
      tr_cluster = self._dfs (ind, ind_grid_bur, mark)
      clusters.append ([i for i in self.ind_grid_vector (tr_cluster) if i > -1])
    self.clusters = clusters
    self.buriedCubesClustered = True
    return clusters

  def _dfs(self, index, ind_of_cubes, mark):
    """ This is used to cluster indices.
    DFS stands for depth first search

    IC - Isn't this actually a breadth-first search?
    """
    cl = set ()
    cl.add (index)
    mark.remove (index)
    toCheck = set (self.get_neighs (index, ind_of_cubes))
    while len (toCheck) > 0:
      n = set ()
      for i in toCheck:
        mark.remove (i)
        cl.add (i)
        n.update (self.get_neighs (i, ind_of_cubes))
      n.difference_update (cl)

      toCheck = n.difference (cl)
    return list (cl)

  def indexToCube(self, indexList):
    """ This function is roughly the inverse of ind_grid_vector .. """
    cubes = []
    n21 = self.ngrid[2]*self.ngrid[1]
    n2 = self.ngrid[2]
    for idx in indexList:
      i0 = int(math.floor(idx/n21))
      idx = idx-i0*n21
      i1 = int(math.floor(idx/n2))
      i2 = int(idx-i1*n2)
      cubes.append( (i0,i1,i2) )
    return cubes

  def cube_to_coord(self, ind_of_cubes):
    """ Convert indices to coordinates """
    return [[self.gridspacing * (ind[i] + 0.5) + self.min[i] for i in (0, 1, 2)]
            for ind in ind_of_cubes]

  def coord_to_cube(self, coords):
    """ convert coordinates to indices """
    return [[int (((coord[i] - self.min[i]) / self.gridspacing)) for i in (0, 1, 2)]
            for coord in coords]