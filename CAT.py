from Treatement_v4 import Treatement
from pylab import *
from os import listdir
from re import search
from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit


ion()
rc('text', usetex=True)
rc('font', family='serif')


class CAT(Treatement):
    """docstring for ."""

    left = []
    right = []
    zero = []
    tunnel = []
    nonTunnel = []

    def sortPopulationBySign(self):
        N = len(self.variableArray)
        self.left, self.right, self.zero = zeros(N), zeros(N), zeros(N)
        for k, v in self.orders.items():
            if k < 0:
                self.left += v
            elif k > 0:
                self.right += v
            elif k == 0:
                self.zero += v

    def plotPopulationBySign(self, save=True):
        fig = figure()
        title("Population of signed momentum")
        xlabel(self.variable)
        ylabel("Normalized density")
        xticks(range(len(self.variableArray)), self.variableArray)
        plot(self.left, 'b+--')
        plot(self.right, 'r+--')
        plot(self.zero, 'y+--')
        legend(['negative momentum', 'positive momentum', 'zero momentum'])
        show()
        if save is True:
            fig.savefig("populationBySign.svg", formt='svg')

    def sortTunnelPopulation(self):
        self.tunnel = []
        self.nonTunnel = []
        for v, i in zip(self.variableArray, range(len(self.variableArray))):
            if v % 2 == 0:
                self.tunnel.append(self.right[i])
                self.nonTunnel.append(self.left[i])
            if v % 2 == 1:
                self.tunnel.append(self.left[i])
                self.nonTunnel.append(self.right[i])

    def plotTunnelPopulation(self, save=True):
        fig = figure()
        title("Population oscillations")
        xlabel(self.variable)
        ylabel("Normalized density")
        xticks(range(len(self.variableArray)), self.variableArray)
        plot(self.tunnel, 'b+--')
        plot(self.nonTunnel, 'r+--')
        legend(['Tunneled population', 'Non tunneled population'])
        show()
        if save is True:
            fig.savefig("tunnelPopulation.svg", formt='svg')
