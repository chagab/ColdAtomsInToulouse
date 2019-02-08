from Treatement_v7 import Treatement
from pylab import *
from scipy.optimize import curve_fit


class Oscillations(Treatement):
    """docstring for Oscillations."""

    guesses = []
    coeffs = []
    matcov = []

    def func(self, t, amp, f, phi, offset):
        # function to fit
        return amp * cos(2 * pi * f * t + phi) + offset

    def fitOscillations(self, harmonics=None):
        t = self.variableArray
        dt = t[1] - t[0]
        data = self.orders[0]

        amp_guess = (data.max() - data.min()) / 2
        weights = abs(fft(data))**2
        index = weights.argmax() if harmonics is None else harmonics
        f_guess = fftfreq(t[-1] - t[0], dt)[index]
        offset_guess = data.mean()
        phi_guess = arccos(data[0] - offset_guess)
        guesses = [amp_guess, f_guess, phi_guess, offset_guess]
        self.guesses = guesses

        coeffs, matcov = curve_fit(self.func, t, data, self.guesses)
        self.coeffs = coeffs
        self.matcov = matcov

    def plotFit(self, plotGuess=False):
        figure()
        t = self.variableArray
        data = self.orders[0]
        ti, tf, dt = t[0], t[-1], t[1] - t[0]
        dt_precise = dt / 1000
        t_precise = arange(ti, tf, dt_precise)
        title(r"""$A\cos(2\pi ft + \phi) + o$
        $A={}, f={}, \phi={}, o={}$""".format(*self.coeffs))
        plot(t, data, 'b+--')
        plot(t_precise, self.func(t_precise, *self.coeffs), 'r')
        if plotGuess:
            plot(t_precise, self.func(t_precise, *self.guesses), 'g--')
        show()
