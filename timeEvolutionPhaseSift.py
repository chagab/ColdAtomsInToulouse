# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 09:42:45 2019

@author: Gabriel Chatelain and Maxime Arnal
"""

import pylab as py
# py.ion()

py.rc('font',**{'family':'serif','serif':['Palatino Linotype'],'size': 14})
#py.rc('text', usetex=True)


class timeEvolutionPhaseShift():
    """docstring for timeEvolution."""

    #Planck constant
    h=6.62607015*1e-34
    # frequency linked to E_L in Hz
    nu_L = 8111.26
    #E_L
    E_L=h*nu_L
    # number of bands
    M = 21
    # number of sites
    Ns = 30
    # length step
    dx = 1/100.0
    ordreMax = 25
    ordreMaxPlot = 4
    # position range
    x = py.arange(-Ns / 2, Ns / 2 + dx, dx)
    # momentum range
    k = py.arange(-ordreMax, ordreMax + 1, 1)
    # quasimomentum range
    q = py.arange(-0.5, 0.5 + 0.001, 0.001)
    # middle index to find q=0
    middleQ = int((len(q)-1)/2 +1)
    # fourier coefficients indexes
    m = py.arange(int(-(M - 1)/2), int((M - 1)/2) + 1)

    def __init__(self, s0=10, angle=25, Texp=40e-6, nstep=40):
        # optical lattice depth
        self.s0 = s0
        # shift angle
        self.angle = angle
        # experimental ime of the experiment
        self.Texp = Texp
        # time step
        self.dt = self.Texp / nstep
        # time range
        self.t = py.arange(0, self.Texp + self.dt / 2, self.dt)[py.newaxis, :]

        self.H2 = py.diag((self.M - 1) * [-self.s0/4], 1)
        self.H3 = py.diag((self.M - 1) * [-self.s0/4], -1)

    def solveCentralEquation(self, q=0):
        self.H1 = py.diag((q + self.m)**2)
        H = self.H1 + self.H2 + self.H3
        self.H = H
        # compute the eigenfunctions and eigenenergies
        self.Eq, self.Cq = py.eigh(H)
#        self.Eq = (self.Eq + self.s / 2)


    def computeBandStructure(self):
        bandsIndex = range(self.M)
        bands = {}
        for i in bandsIndex:
            bands[i] = py.zeros(len(self.q))

        for i, q in enumerate(self.q):
            self.solveCentralEquation(q)
            for j in bandsIndex:
                bands[j][i] = self.Eq[j]
        self.bands = bands

        freqOscill=(bands[2][self.middleQ]-bands[0][self.middleQ])*self.nu_L/2
        freqOscill=round(freqOscill,2)

        print('s=' + str(self.s0))
        print('Ligne 2 phonons, frequence oscillations = ' + str(freqOscill) + ' Hz.')


    def computeEigenFunctions(self):
        self.solveCentralEquation()
        # generate matrix with M x-rows
        XX = (py.ones([len(self.x), self.M]).T @ py.diag(self.x)).T
        argument = 1j * XX @ py.diag(2 * py.pi * self.m)
        self.argument = argument
        eigenfuncs=[]
        for l in range(self.M):
            eigenf = py.sum(py.exp(argument) @ py.diag(
                self.Cq[:, l]), axis=1).T
            if l==0:
                eigenfuncs = eigenf
            else:
                eigenfuncs=py.vstack((eigenfuncs,eigenf))
        self.eigenfuncs = eigenfuncs


    def computeTimeEvolution(self):
        #step 1: calculate eigenfunctions
        self.computeEigenFunctions()
        #step 2: phase got by fourier coefficients due to translation along x-axis
        phase = py.diag(py.exp(-1j * self.angle * 2 * py.pi / 180 * self.m))
        phaseFactor = []
        for l in range(self.M):
            # Dot product between Bloch function l and shifted one
            phaseFactor.append(self.Cq[:, l] @ phase @ self.Cq[:, 0])
        #step 3: omegas are
        self.omegas = (self.Eq * 2 * py.pi * self.nu_L)[:, py.newaxis]
        #step 4: phase at each time
        self.timePhase = py.exp(-1j * self.t.T @ self.omegas.T)
        #step 5: calculate time evolution
        self.timeEvolution = self.timePhase \
            @ py.diag(phaseFactor) \
            @ self.eigenfuncs
        #step 6: space density
        self.density = abs(self.timeEvolution)**2

    def computeMomentumEvolution(self):
        self.computeTimeEvolution()
        self.momentumPhase = py.exp(-1j * 2 * py.pi * py.diag(self.x)
                                    @ py.ones([len(self.x), len(self.k)])
                                    @ py.diag(self.k))
        self.momentumEvolution = self.timeEvolution @ self.momentumPhase
        self.momentumDensity = abs(self.momentumEvolution)**2
        # normalisation of momentum density
        self.momentumDensity = self.momentumDensity/ \
            py.sum(self.momentumDensity,axis=1)[:,None]

    def plotTimeEvolution(self):
        fig = py.figure()
        j1 = 1400
        j2 = 1600
        for j in range(py.shape(self.t)[1]):
            py.plot(self.x[j1: j2], self.density[j, j1:j2], 'b')
            py.ylim([0, 0.2])
            py.pause(0.01)
            fig.clear()
        py.show()

    def plotMomentumEvolution(self):
        self.computeMomentumEvolution()
        py.figure()
        py.imshow(self.momentumDensity, cmap='jet', aspect='auto')
        py.colorbar()
        py.xlabel('k/kL')
        py.ylabel('temps (µs)')
        py.show()

    def plotOrderEvolution(self):
        self.computeMomentumEvolution()
        numbLines=int(2*self.ordreMaxPlot/3)+1
        py.figure()
        for order in py.arange(-self.ordreMaxPlot, self.ordreMaxPlot+1, 1):
            middleK = int(len(self.k) / 2)
            self.middleK = middleK
            ax = py.subplot(numbLines, 3, order + self.ordreMaxPlot + 1)
            py.plot(self.t.T / 1e-6,
                    self.momentumDensity[:, middleK+order], 'b')
            ax.set_xlim([0,None])
            ax.set_ylim([0,None])
            py.title('Ordre ' + str(order))
            py.xlabel('Temps (µs)')
            py.ylabel('Pop')
        py.show()
        py.tight_layout()

    def plotBandStructure(self, NbandsToPlot=8, typeOfModulation='phase',
                          freq=0, deltaFreq=100, color=False):
        self.computeBandStructure()
        freq=freq/self.nu_L;
        deltaFreq=deltaFreq/self.nu_L;

        fig, ax = py.subplots()

        if color:
            py.title('Depth: ' + str(self.s0)
                     + '\n Freq: ' + str(round(freq*self.nu_L/1000,2)) + ' kHz'
                     + '\n Modulation: ' + typeOfModulation)
            uncertainty = deltaFreq
            bandsIndex = range(NbandsToPlot)
            for i in bandsIndex:
                bandColor = self.chooseColor(
                    self.bands[i], i + 1, typeOfModulation)
                ax.plot(self.q, self.bands[i], bandColor)
                ax.annotate(str(i + 1), (self.q[self.middleQ + int(len(self.q)/4)],
                            self.bands[i][self.middleQ + int(len(self.q)/4)]))
                if freq != 0:
                    for j in bandsIndex:
                        if j > i:
                            transitions = self.bands[j] - self.bands[i]
                            for k, transition in enumerate(transitions):
                                if(transition > (freq - uncertainty)
                                   and transition < (freq + uncertainty)):
                                    ax.arrow(self.q[k], self.bands[i][k], 0,
                                             self.bands[j][k] - self.bands[i][k],
                                             head_width=0.05, head_length=0.1, fc='k',
                                             ec='k')
        else:
           py.title('Depth: ' + str(self.s0) )
           bandsIndex = range(NbandsToPlot)
           for i in bandsIndex:
                ax.plot(self.q, self.bands[i])
                ax.annotate(str(i + 1), (self.q[self.middleQ + int(len(self.q)/4)],
                            self.bands[i][self.middleQ + int(len(self.q)/4)]))

        ax.set_xlabel('Quasimomentum (units of k_L/2)')
        ax.set_ylabel('Energy (units of E_L)')
        ax.set_xlim([-0.5,0.5])
        py.show()

    def chooseColor(self, band, bandNumber, typeOfModulation):
        if max(band) < 0:
            if typeOfModulation == 'amplitude':
                if bandNumber % 2 == 0:
                    return '--g'
                else:
                    return 'g'
            if typeOfModulation == 'phase':
                if bandNumber % 2 == 1 and bandNumber != 1:
                    return '--g'
                else:
                    return 'g'
        else:
            if typeOfModulation == 'amplitude':
                if bandNumber % 2 == 0:
                    return '--r'
                else:
                    return 'r'
            if typeOfModulation == 'phase':
                if bandNumber % 2 == 1 and bandNumber != 1:
                    return '--r'
                else:
                    return 'r'


if __name__ is '__main__':

    py.close('all')
    t = timeEvolutionPhaseShift(s0=21, angle=45, nstep=1000, Texp=40e-6)
#    t.computeBlochFunctions()
#    t.computeEigenFunctions()
#
#    t.computeTimeEvolution()
#    t.computeMomentumEvolution()
#
#    t.computeBandStructure()
    t.plotBandStructure(color=False)
#
    t.plotMomentumEvolution()
    t.plotOrderEvolution()
#    t.plotTimeEvolution()