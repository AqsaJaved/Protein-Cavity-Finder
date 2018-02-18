import os
import sys
from re import *

def createSubDirect(root,sub):
  sub2 = os.path.basename(sub)
  newdir = root + '/' + sub2
  if not os.path.isdir(newdir):
    try:
      os.mkdir(newdir)
    except:
      print 'ERROR:\tModul: myTools.py\tfunction: createSubDirect()'
      print '\tNot able to create directory [' + root + '/' + sub+']'
  return root + '/' + sub2