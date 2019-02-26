#!/usr/bin/env python3
import sys
from pylab import *
close('all')


def Treatement(arg):
    path = '/media/gabriel/MyPassport/theseGab/TraitementData/sources'
    sys.path.append(path)
    from Treatement import Treatement
    return Treatement(arg)


t = Treatement("tmaintien")
t.autoTreatement()
