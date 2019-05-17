from Treatement import Treatement
import pylab as py

py.ion()
py.rc('text', usetex=True)
py.rc('font', family='serif')


class CAT(Treatement):
    """docstring for ."""

    left = []
    right = []
    zero = []
    tunnel = []
    nonTunnel = []

    def sortPopulationByMomentumSign(self):
        """
        Argument :
            - None.
        Return :
            - None.

        Fills in two array that correspond to the population with negative
        momentum (left array) and positive momentum (right array). The sign is
        given by the position of the excited cloud with respect to the central
        cloud.
        """
        n = len(self.variableArray)
        self.left, self.right, self.zero = (py.zeros(n),
                                            py.zeros(n), py.zeros(n))
        # the iteration is done on each orders where the keys "k" of the dict
        # "self.orders" give the momentum's sign  of the corresponding
        # order.
        for k, v in self.orders.items():
            if k < 0:
                self.left += v
            elif k > 0:
                self.right += v
            elif k == 0:
                self.zero += v

    def plotPopulationByMomentumSign(self, save=True):
        """
        Argument :
            - save : boolean. If True, the resulting figure is saved. Default
            behavior is True.
        Return :
            - None.

        Plot the evolution of the population of atoms with respect to it's
        momentum sign.
        """
        self.sortPopulationByMomentumSign()
        fig = py.figure()
        py.title("Population of signed momentum")
        py.xlabel(self.variable)
        py.ylabel("Normalized population")
        py.plot(self.variableArray, self.left, 'b+--')
        py.plot(self.variableArray, self.right, 'r+--')
        py.plot(self.variableArray, self.zero, 'y+--')
        py.legend(['$p<0$', '$p>0$', '$p=0$'])
        py.grid()
        py.show()
        if save is True:
            fig.savefig("populationByMomentumSign.svg", formt='svg')

    def sortTunnelPopulation(self):
        """
        Argument :
            - None.
        Return :
            - None.

        Fills in two array that correspond to the population that has tunneled
        during the experiment or not. To determine wether a cloud of atoms has
        tunneled or not, the following reasonning is applyied :
            - On an odd number of periods, non-tunneled have a positive
            momentum after the phase space rotation.
            - On an even number of periods, non-tunneled have a negative
            momentum after the phase space rotation.
        """
        self.sortPopulationByMomentumSign()
        self.tunnel = []
        self.nonTunnel = []
        for i, v in enumerate(self.variableArray):
            if v % 2 == 0:
                self.tunnel.append(self.right[i])
                self.nonTunnel.append(self.left[i])
            if v % 2 == 1:
                self.tunnel.append(self.left[i])
                self.nonTunnel.append(self.right[i])
        self.tunnel = py.array(self.tunnel)
        self.nonTunnel = py.array(self.nonTunnel)

    def plotTunnelPopulation(self, save=True):
        """
        Argument :
            - save : boolean. If True, the resulting figure is saved. Default
            behavior is True.
        Return :
            - None.

        Plot the evolution of the tunneled and non-tunneled population.
        """
        self.sortTunnelPopulation()
        fig = py.figure()
        py.title("Population oscillations")
        py.xlabel(self.variable)
        py.ylabel("Normalized population")
        py.plot(self.variableArray, self.tunnel, 'r+--')
        py.plot(self.variableArray, self.nonTunnel, 'b+--')
        py.plot(self.variableArray, self.zero, 'y+--')
        py.legend(['Tunneled population', 'Non-tunneled population',
                   'zeroth-order population'])
        py.grid()
        py.show()
        if save is True:
            fig.savefig("tunnelPopulation.svg", formt='svg')

    def computeMomentumSpreading(self):
        spreading = []
        meanPosMomentum = []
        for i in range(len(self.variableArray)):
            profile = self.profileArray[i] / max(self.profileArray[i])
            maxLocation = []
            maxVal = []
            for x0x1 in self.coords:
                x0 = x0x1[0]
                x1 = x0x1[1]
                maxLocation.append(x0 + py.argmax(profile[x0:x1]))
                maxVal.append(max(profile[x0:x1]))

            maxLocation = py.array(maxLocation)
            maxVal = py.array(maxVal)
            zeroIndex = int(len(maxLocation) / 2)
            momemtum = (py.array(range(len(profile))) -
                        maxLocation[zeroIndex]) / self.h_d_ratio
            maxLocationMomentum = (py.array(maxLocation)
                                   - maxLocation[zeroIndex]) / self.h_d_ratio
            meanPosMomentum.append(
                py.sum(maxLocationMomentum * maxVal) / py.sum(maxVal))
            spreading.append(
                py.sum(maxVal * (maxLocationMomentum - meanPosMomentum[-1])**2))
        spreading = py.array(spreading)
        spreading = spreading / max(spreading)
        self.spreading = spreading
        self.meanPosMomentum = meanPosMomentum

    def plotMeanMomentum(self):
        py.figure()
        py.ylabel('p*')
        py.xlabel(self.variable)
        py.plot(self.variableArray, self.meanPosMomentum, 'b+--')
        py.grid()
        py.show()

    def plotMomentumSpreading(self):
        py.figure()
        py.ylabel('Spreading (arb. units)')
        py.xlabel(self.variable)
        py.plot(self.variableArray, self.spreading, 'b+--')
        py.grid()
        py.show()
