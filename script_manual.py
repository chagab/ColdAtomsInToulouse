#!/usr/bin/env python3
import sys
from pylab import *
close('all')


def CAT(arg):
    path = '/media/gabriel/MyPassport/theseGab/TraitementData/sources'
    sys.path.append(path)
    from CAT import CAT
    return CAT(arg)


t = CAT("nbrePeriodes")
t.setNoisyArea(noise=(480, 680))
t.setAreaOfInterest(angle=4, area=(100, 1292, 550, 730))
t.computeAllOD()
t.computeAllProfile()
ordres = [[85 + i * 70, 155 + i * 70] for i in range(15)]
t.plotAllODAtOnce(ordres)
t.computeEvolutionOfOrder()
t.plotAllOrderAtOnce()
t.sortPopulationBySign()
t.plotPopulationBySign()
t.sortTunnelPopulation()
t.plotTunnelPopulation()
