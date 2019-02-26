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

    def sortPopulationBySign(self):
        N = len(self.variableArray)
        self.left, self.right, self.zero = (py.zeros(N),
                                            py.zeros(N), py.zeros(N))
        for k, v in self.orders.items():
            if k < 0:
                self.left += v
            elif k > 0:
                self.right += v
            elif k == 0:
                self.zero += v

    def plotPopulationBySign(self, save=True):
        self.sortPopulationBySign()
        fig = py.figure()
        py.title("Population of signed momentum")
        py.xlabel(self.variable)
        py.ylabel("Normalized density")
        py.xticks(range(len(self.variableArray)), self.variableArray)
        py.plot(self.variableArray, self.left, 'b+--')
        py.plot(self.variableArray, self.right, 'r+--')
        py.plot(self.variableArray, self.zero, 'y+--')
        py.legend(['negative momentum', 'positive momentum', 'zero momentum'])
        py.show()
        if save is True:
            fig.savefig("populationBySign.svg", formt='svg')

    def sortTunnelPopulation(self):
        self.sortPopulationBySign()
        self.tunnel = []
        self.nonTunnel = []
        for i, v in enumerate(self.variableArray):
            if v % 2 == 0:
                self.tunnel.append(self.right[i])
                self.nonTunnel.append(self.left[i])
            if v % 2 == 1:
                self.tunnel.append(self.left[i])
                self.nonTunnel.append(self.right[i])

    def plotTunnelPopulation(self, save=True):
        self.sortTunnelPopulation()
        fig = py.figure()
        py.title("Population oscillations")
        py.xlabel(self.variable)
        py.ylabel("Normalized density")
        py.plot(self.variableArray, self.tunnel, 'b+--')
        py.plot(self.variableArray, self.nonTunnel, 'r+--')
        py.legend(['Tunneled population', 'Non tunneled population'])
        py.show()
        if save is True:
            fig.savefig("tunnelPopulation.svg", formt='svg')
