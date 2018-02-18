import sys
import os
import re
import time
from optparse import OptionParser
from PocketStream import PocketStream
from MoleculeSource import MoleculeSource

def parseCommandLineOptions():

    parser = OptionParser()
    parser.add_option('-i', '--input', dest='pdbfile', help='Pdb file')
    #parser.add_option ('-d', '--degree', dest='degree_buriedness', help='Degree of buriedness')
    #parser.add_option ('-m', '--minneigh', dest='min_neighbour', help='Minimum number of buried neighbours')
    #parser.add_option ('-c', '--minclustsize', dest='min_clust_size', help='Minimum cluster size')

    (options, args) = parser.parse_args()

    if options.pdbfile is None:
        options.pdbfile = raw_input('Enter PDB file:')
    #if options.degree_buriedness is None:
    #    options.degree_buriedness = raw_input('Enter Degree of Buriedness:')
    #if options.min_neighbour is None:
    #    options.min_neighbour = raw_input('Enter minimum no of buried neighbours:')
    #if options.min_clust_size is None:
    #    options.min_clust_size = raw_input('Enter minimum cluster size:')

    message = 'The pdb file is ' + options.pdbfile + '.'

    print message
    return options
class CavFinder:
    # Hard coding parameters as of now
    def __init__(self, options):
        self.pdb = options.pdbfile
        self.out_directory= os.getcwd()
        self.bind_region = None
        self.min_neighbour=15
        self.degree_buriedness=11
        self.min_clust_size=50
        #self.min_neighbour= options.min_neighbour
        #self.degree_buriedness = options.degree_buriedness
        #self.min_clust_size = options.min_clust_size
        self.ncl = 0 # ncl for no of clusters
        self.includeH= False # ignoring Hydrogen for now

        # Check for validity of pdb file
        if not os.path.isfile (self.pdb):
            print "File not found in current directory: "
            sys.exit(0)

        self.pdbPattern = re.compile ("^[^#].*pdb$")
        if not self.pdbPattern.search(self.pdb):
            print "Enter a file ending with .pdb: "
            sys.exit (0)

        # Perform pocket analysis
        self.molsource = MoleculeSource (pdbfile=self.pdb, tempDir=self.out_directory,
                                         includeH=self.includeH)
        self.pokstream = PocketStream (pdbfile=self.pdb,
                                       outputDir=self.out_directory,
                                       includeH=self.includeH,
                                       minNeighbours=self.min_neighbour,
                                       dobThreshold=self.degree_buriedness,
                                       minClust=self.min_clust_size,
                                       )
        

def main(argv):
    t0 = time.time ()
    if sys.platform != 'win32': os.system ('clear')

    options = parseCommandLineOptions ()
    CavFinder (options)

    tdiff = time.time () - t0
    print 'Overall time: %8.2f s' % tdiff

    sys.exit (0)
###########################################

if __name__ == '__main__':
    main (sys.argv)