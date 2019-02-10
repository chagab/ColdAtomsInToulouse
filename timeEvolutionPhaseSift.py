from pylab import *
ion()


class timeEvolutionPhaseShift():
    """docstring for timeEvolution."""

    # frequency linked to E_L in Hz
    nu_L = 8111.26
    # number of bands
    M = 21
    Q = 2
    # number of sites
    Ns = 30
    # length step
    dx = 0.01
    # momentum step
    dk = 0.01
    tmax = 1
    ordreMax = 10
    ordreMaxPlot = 3
    # position range
    x = arange(-Ns / 2, Ns / 2 + dx, dx)
    # momentum range
    k = arange(-ordreMax, ordreMax + dk, dk)
    # quasimomentum range
    q = arange(-0.5, 0.5 + 0.1, 0.1)
    # fourier coefficients indexes
    m = arange(int(-(M - 1) / 2), int((M - 1) / 2) + 1)

    def __init__(self, s0=10, angle=25, Texp=40e-6, nstep=40):
        # optical lattice depth
        self.s0 = s0
        self.s = -self.s0 * 4
        self.Vq = -self.s / 4
        # shift angle
        self.angle = angle
        # experimental ime of the experiment
        self.Texp = Texp
        # time step
        self.dt = 1 / nstep
        # time range
        self.t = arange(0, self.tmax + self.dt / 2, self.dt)[newaxis, :]

    def computeBandStructure(self):
        bandsIndex = range(self.M)
        bands = {}
        for i in bandsIndex:
            bands[i] = zeros((len(self.q)))

        H2 = diag((self.M - 1) * [self.Vq], 1)
        H3 = diag((self.M - 1) * [self.Vq], -1)

        for i, q in enumerate(self.q):
            H1 = diag((q - self.m)**2)
            H = H1 + H2 + H3
            Eq, Cq = eigh(H)
            Eq = (Eq + self.s / 2) * self.nu_L / 4
            for j in bandsIndex:
                bands[j][i] = Eq[j]
        self.bands = bands

    def computeBlochFunctions(self):
        # compute the Bloch functions at the center of the Brillouin zone
        q = 0
        # construct the Bloch Hamiltonian
        H1 = diag((q - self.m)**2)
        H2 = diag((self.M - 1) * [self.Vq], 1)
        H3 = diag((self.M - 1) * [self.Vq], -1)
        H = H1 + H2 + H3
        # compute the eigenfunctions and eigenenergies
        Eq, Cq = eigh(H)
        # normalize eigenenergies with respect to nu_L
        Eq = (Eq + self.s / 2) * self.nu_L / 4
        self.Eq = Eq
        self.Cq = Cq
        # the omegas are
        self.omegas = (Eq * 2 * pi * self.Texp)[:, newaxis]

    def computeEigenFunctions(self):
        # generate matrix with M x-rows
        XX = (ones([len(self.x), self.M]).T @ diag(self.x)).T
        argument = 1j * XX @ diag(2 * pi * self.m)
        eigenfuncs = []
        for l in range(self.M):
            eigenf = sum(exp(argument) @ diag(self.Cq[:, l]), axis=1).T
            eigenf /= sqrt(self.Ns)
            eigenfuncs.append(eigenf)
        self.eigenfuncs = eigenfuncs

    def computePhaseFactor(self):
        # phase got by fourier coefficients due to translation along x-axis
        phase = diag(exp(-1j * self.angle * 2 * pi / 180 * self.m))
        phaseFactor = []
        for l in range(self.M):
            # Dot product between Bloch function l and shifted one
            phaseFactor.append(self.Cq[:, l].T @ phase @ self.Cq[:, 1])
        self.phaseFactor = phaseFactor

    def computeTimeEvolution(self):
        self.computePhaseFactor()
        self.timePhase = exp(-1j * self.t.T @ self.omegas.T)
        self.timeEvolution = self.timePhase \
            @ diag(self.phaseFactor) \
            @ self.eigenfuncs
        self.density = abs(self.timeEvolution)**2

    def computeMomentumEvolution(self):
        self.momentumPhase = exp(-1j * 2 * pi * diag(self.x)
                                 @ ones([len(self.x), len(self.k)])
                                 @ diag(self.k))
        self.momentumEvolution = self.timeEvolution @ self.momentumPhase
        # normalisation of momentum
        self.momentumEvolution *= self.dx / sqrt(2 * pi)
        self.momentumDensity = abs(self.momentumEvolution)**2
        X = array(sum(self.momentumDensity, axis=1)**(-1))[:, newaxis]
        # normalisation of momentum density
        self.momentumDensity = self.momentumDensity * \
            concatenate(len(self.k) * [X], axis=1) / self.dk

    def plotTimeEvolution(self):
        fig = figure()
        j1 = 1400
        j2 = 1600
        for j in range(shape(self.t)[1]):
            plot(self.x[j1: j2], self.density[j, j1:j2], 'b')
            ylim([0, 0.2])
            pause(0.01)
            fig.clear()
        show()

    def plotMomentumEvolution(self):
        figure()
        imshow(self.momentumDensity, cmap='jet', aspect='auto')
        colorbar()
        xlabel('k/kL')
        ylabel('temps (µs)')
        show()

    def plotOrderEvolution(self):
        figure()
        for order in arange(-self.ordreMaxPlot, self.ordreMaxPlot, 1):
            middleK = floor(len(self.k) / 2) + 1
            subplot(2, 4, order + self.ordreMaxPlot + 1)
            slice = int(middleK + fix(order / self.ordreMax * middleK))
            plot(self.t.T / 1e-6 * self.Texp,
                 self.momentumDensity[:, slice], 'b')
            title('Ordre ' + str(order))
            xlabel('time (µs)')
            ylabel('density')
        tight_layout()
        show()

    def plotBandStructure(self,
                          NbandsToPlot=8,
                          typeOfModulation='phase',
                          freq=0,
                          deltaFreq=5e2):
        fig, ax = subplots()
        title('depth: ' + str(self.s0) + '\n freq: ' + str(freq) +
              ' kHz' + '\n modulation: ' + typeOfModulation)
        uncertainty = deltaFreq
        bandsIndex = range(NbandsToPlot)
        l = len(self.q)
        for i in bandsIndex:
            bandColor = self.chooseColor(
                self.bands[i], i + 1, typeOfModulation)
            ax.plot(self.q, self.bands[i], bandColor)
            ax.annotate(str(i + 1), (self.q[l - 7], self.bands[i][l - 7]))
            for j in bandsIndex:
                if j > i:
                    transitions = self.bands[j] - self.bands[i]
                    for transition, k in zip(transitions, range(len(transitions))):
                        if(transition > (freq - uncertainty) and transition < (freq + uncertainty)):
                            ax.arrow(self.q[k], self.bands[i][k], 0, self.bands[j][k] - self.bands[i][k], head_width=0.05,
                                     head_length=0.1, fc='k', ec='k')
        show()

    def chooseColor(slef, band, bandNumber, typeOfModulation):
        if max(band) < 0:
            if typeOfModulation is 'amplitude':
                if bandNumber % 2 is 0:
                    return '--g'
                else:
                    return 'g'
            if typeOfModulation is 'phase':
                if bandNumber % 2 is 1 and bandNumber is not 1:
                    return '--g'
                else:
                    return 'g'
        else:
            if typeOfModulation is 'amplitude':
                if bandNumber % 2 is 0:
                    return '--r'
                else:
                    return 'r'
            if typeOfModulation is 'phase':
                if bandNumber % 2 is 1 and bandNumber is not 1:
                    return '--r'
                else:
                    return 'r'


close('all')
t = timeEvolutionPhaseShift(s0=10, angle=90, nstep=54, Texp=16.3e-6)
t.computeBlochFunctions()
t.computeEigenFunctions()

t.computeTimeEvolution()
t.computeMomentumEvolution()

t.computeBandStructure()
t.plotBandStructure()

t.plotMomentumEvolution()
t.plotOrderEvolution()
# t.plotTimeEvolution()
