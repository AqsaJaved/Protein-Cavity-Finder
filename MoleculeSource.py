import os,sys,re, shutil

class MoleculeSource:
    def __init__(self, pdbfile=None, tempDir=".", includeH=False):
        """ First, process the input files to create a dictionary of file locations
            for setting up the various MoleculeStreams.

            Second, initialise a list of MoleculeStream objects and a pointer to
            indicate which is currently in use.

            Third, set up a single projMask list covering all MoleculeStreams."""

        self.pdbfile = pdbfile
        self.tempDir = tempDir
        self.includeH = includeH