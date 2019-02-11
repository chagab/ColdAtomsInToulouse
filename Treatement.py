#######################################################################
# Program writen by Gabriel Chatelain at laboratory LCAR in Toulouse #
#######################################################################

import pylab as py
from os import listdir
from re import search
from scipy.ndimage.interpolation import rotate


py.ion()
py.rc('text', usetex=True)
py.rc('font', family='serif')


class Treatement(object):
    """docstring for Treatement."""
    # The path that points to the images folder.
    path = ''
    # The parameter that is varried during the experiment.
    variable = ''
    # The array that contains the value of the varied parameter.
    variableArray = []
    # The array that contains all complet OD images.
    ODArray = []
    # The array that contains all croped OD images.
    ODArrayCroped = []
    # The array that contains the y-axis integrated OD profile.
    profileArray = []
    # The array that contains the coordinate of each order of diffraction.
    coords = []
    # The dictionnary that contains the evolution of the order during the
    # experiment.
    orders = {}
    # The tuple that contains the coordinate of the noisy area of
    # the OD imaages.
    noise = ()
    # The array that contains the noisy area (area without the atoms) used to
    # comute the value of the background.
    noiseArray = []

    def __init__(self, variable, path="../Images/"):
        """
        Argument :
            - variable : string. The parameter that is varied during the
            experiment.
            - path : string. The path from this script to the images folder.
        Return :
            - None
        """
        self.variable = variable
        self.path = path

    def setAngle(self, angle):
        self.angle = angle

    def setArea(self, area):
        self.area = area

    def setOrdes(self, coords):
        self.coords = coords

    def calcOD(self, wiAt, noAt, back):
        """
        Argument :
            - wiAt : 2D array. Image of the atoms.
            - noAt : 2D array. Image without the atoms but with the laser on.
            - back : 2D array. Image without the atoms and without the laser
        Return :
            - OD : 2D array. Optical density (OD) image.

        Compute an OD image from three images. The OD is deduced using the Beer
        Lambert law : I = I_0 exp(-OD). In order to have a clean signal, one
        can remove the background noise "b" wich is given by the "back" image.
        Thus, OD is given by
            OD = -log((I - b)/(I_0 - b))
        """
        # remove the background
        AT = wiAt - back
        NO = noAt - back
        # change the values that can cause some problems with the log and
        # division
        AT[AT <= 1] = 1
        NO[NO <= 1] = 1
        # compute the optical density
        OD = -py.log(AT / NO)
        return OD

    def computeAllOD(self, area=(0, 1392, 0, 1040), angle=0):
        """
        Argument:
            - area : tuple. AOI in the form (x0, x1, y0, y1).
            - angle : integer. Angle of rotation of the OD image.
        Return:
            - None

        Generate an array of optical denisty images from
         a folder of tiff images.

         _Note_ : the "path" attribute has to be set before calling this
         method. """
        # unpack the area argument
        try:
            area = self.area
            angle = self.angle
        except AttributeError:
            pass
        x0, x1, y0, y1 = area
        variables = []
        wiAt, noAt, back = {}, {}, {}
        wiAtOffset, noAtOffset, backOffset = {}, {}, {}
        print("Begin scan")
        for f in listdir(self.path):
            # get the variable that is varaied
            i = int(  # convert the result to an integer
                # extract the value of the variable with regular expression
                search(self.variable + "([0-9p]+)", f)
                .group()
                .replace(self.variable, "")
                .replace("p", ".")
            )
            variables.append(i)
            # set the coordinate of where the noise is
            noiseRow = [i for i in range(
                1040)if i < self.noise[0] or i > self.noise[1]]
            path = self.path + f
            # crop, reshape the image and dispatch it in the correct list
            if f.endswith("_With.tif"):
                wiAt[i] = py.array(py.imread(path))
                wiAtOffset[i] = py.array(py.imread(path))[noiseRow, x0:x1]
            if f.endswith("_NoAt.tif"):
                noAt[i] = py.array(py.imread(path))
                noAtOffset[i] = py.array(py.imread(path))[noiseRow, x0:x1]
            if f.endswith("_Back.tif"):
                back[i] = py.array(py.imread(path))
                backOffset[i] = py.array(py.imread(path))[noiseRow, x0:x1]
        # filter the doublons from the variable array
        self.variableArray = list(sorted(set(variables)))
        print("Scan done")
        # In case of any changes of self.angle and self.area, the
        # self.ODArrayCroped is reset to an empty list
        self.ODArrayCroped = []
        # from the lists of images, create a list of optical density image
        print("Begin computing")
        for i in self.variableArray:
            OD = self.calcOD(wiAt[i], noAt[i], back[i])
            offset = self.calcOD(wiAtOffset[i], noAtOffset[i], backOffset[i])
            self.noiseArray.append(offset)
            OD -= offset.mean()
            OD[OD < 0] = 0
            self.ODArray.append(OD)
            ODcroped = rotate(OD, angle)[y0:y1, x0:x1]
            self.ODArrayCroped.append(ODcroped)
        print("Computing done")

    def computeFirstOD(self, area=(0, 1392, 0, 1040), angle=0):
        """
        Argument:
            - area : tuple. AOI in the form (x0, x1, y0, y1).
            - angle : integer. Angle of rotation of the OD image.
        Return:
            - 2D array.

        Generate the first optical denisty image from
         a folder of tiff images.

         _Note_ : the "path" attribute has to be set before calling this
         method. """
        wiAt, noAt, back = [], [], []
        files = listdir(self.path)[0:4]
        x0, x1, y0, y1 = area
        for f in files:
            if f.endswith("_With.tif"):
                wiAt = rotate(py.array(py.imread(self.path + f)),
                              angle)[y0:y1, x0:x1]
            if f.endswith("_NoAt.tif"):
                noAt = rotate(py.array(py.imread(self.path + f)),
                              angle)[y0:y1, x0:x1]
            if f.endswith("_Back.tif"):
                back = rotate(py.array(py.imread(self.path + f)),
                              angle)[y0:y1, x0:x1]
        OD = self.calcOD(wiAt, noAt, back)
        return OD

    def plotFirstOD(self,
                    area=(0, 1392, 0, 1040),
                    angle=0, noise=None, theTitle=None):
        """
        Argument:
            - area : 1D array. AOI in the form (x0, x1, y0, y1).
            - angle : integer. Angle of rotation of the OD image.
            - noise : 1D array. Delimiation of the noisy area of the OD image.
            - theTitle : string. Title of the resulting figure.
        Return:
            - None.

        Used to plot the resulting OD image from the "computeFirstOD"
        method."""

        OD = self.computeFirstOD(area, angle)
        self.setAngle(angle)
        self.setArea(area)
        py.figure()
        py.imshow(OD, cmap='jet')
        if noise is not None:
            y0, y1 = noise
            py.axhline(y=y0, color='red', linestyle='--')
            py.axhline(y=y1, color='red', linestyle='--')
            self.noise = noise
        py.xlabel("momentum")
        py.title(theTitle)
        py.colorbar()
        py.show()

    def findNoisyArea(self):
        """
        Argument :
            - None
        Return :
            - None

        Use to detect the noisy area (the area without atoms) by an edge
        detection technique.
        """
        OD = self.computeFirstOD()
        yProfile = py.sum(OD, axis=1)
        derivative = py.gradient(yProfile)
        N = 10
        # because the derivative is usually very noisy, a sliding average is
        # performed in order to smooth the signal. This is done by a
        # convolution with a gate function of size "N".
        res = py.convolve(derivative, py.ones((N,)) / N, mode='valid')
        mean = res.mean()
        # index of the maximum value of the signal.
        i = res.argmax()
        while res[i] >= mean:
            # once "i" is greater or equal to the mean of the derivative,
            # we have found the upper bound of the noisy area.
            i -= 1
        # index of the minimum value of the signal.
        upBound = i - 50
        i = res.argmin()
        while res[i] < mean:
            # once "i" is smaller or equal to the mean of the derivative,
            # we have found the lower bound of the noisy area.
            i += 1
        downBound = i + 50
        self.setNoisyArea((upBound, downBound))
        return OD

    def setNoisyArea(self, noise):
        self.plotFirstOD(noise=noise, theTitle="Noisy area")

    def findAreaOfInterest(self):
        """
        Argument :
            - None
        Return :
            - None

        Use to detect the AOI by an edge detection technique.
        """
        OD = self.findNoisyArea()
        max = 0
        bestAngle = None
        for angle in range(-10, 10, 1):
            # the best angle is the one for wich the maximum peak is reached
            # on the yProfile (integral on the horizontal direction) ...
            image = rotate(OD, angle)
            yProfile = py.sum(image, axis=1)
            newMax = yProfile.max()
            if newMax > max:
                max = newMax
                bestAngle = angle
        # ... once found, the resulting OD image is kept and used to find the
        # the top and bottom bounds by edge detection.
        bestOD = rotate(OD, bestAngle)
        YProfile = py.sum(bestOD, axis=1)
        derivative = py.gradient(YProfile)
        N = 10
        # because the derivative is usually very noisy, a sliding average is
        # performed in order to smooth the signal. This is done by a
        # convolution with a gate function of size "N".
        res = py.convolve(derivative, py.ones((N,)) / N, mode='valid')
        mean = res.mean()
        # index of the maximum value of the signal.
        i = res.argmax()
        while res[i] > mean:
            # once "i" is greater or equal to the mean of the derivative,
            # we have found the upper bound of the AOI.
            i -= 1
        # for security we take an extra 50 pixels.
        y0 = int(i - 50)
        # index of the minimum value of the signal.
        i = res.argmin()
        while res[i] < mean:
            # once "i" is smaller or equal to the mean of the derivative,
            # we have found the lower bound of the AOI.
            i += 1
        # Again, for security, we take an extra 50 pixels.
        y1 = int(i + 50)
        # The horizontal bound are taken to be maximal, but for security, we
        # take an extra 50 pixels.
        x0 = 50
        x1 = py.shape(OD)[0] - 50
        self.setAreaOfInterest(area=(x0, x1, y0, y1), angle=bestAngle)
        return bestOD[y0:y1, x0:x1]

    def findOrders(self, n=9, sep=70, width=0.4):
        """
        Argument :
            - n : integer. Number of orders seen during the experiment.
            - sep : integer. Separation between each orders.
            - width : float. Width of each order in percentage of the
            separation between the orders.
        Return :
            - None.

        Tries to find the location of each order. First, it find the AOI.
        Then it compute each OD in order to plot all the OD at once for the
        user to check if the location is correct.
        """
        OD = self.findAreaOfInterest()
        xProfile = py.sum(OD, axis=0)
        max = xProfile.argmax()
        w = width * sep
        orders = {}
        for i in range(-int(n / 2), int(n / 2) + 1):
            orders[i] = (int(max - i * sep - w / 2),
                         int(max - i * sep + w / 2))
        self.computeAllOD()
        self.computeAllProfile()
        self.plotAllODAtOnce([v for k, v in orders.items()])

    def setAreaOfInterest(self, area=(0, 1392, 0, 1040), angle=0):
        """
        Argument:
            - area : 1D array. AOI in the form (x0, x1, y0, y1).
            - angle : integer. Angle of rotation of the OD image.
        Return:
            - None.

        Wrap of the "plotFirstOD" method. Used to set the "area" attribut."""

        t = "angle : {}, area : {}".format(angle, area)
        self.plotFirstOD(area=area, angle=angle, theTitle=t)

    def plotOD(self, index=None, coords=None, croped=True):
        """Method used to show optical density image.
           If 'index' is given, it shows the corresponding density images.
           If 'index' is 'None', it shows all the  images.

           _Note_ : 'generateAllOD' has to be called first."""

        if croped is True:
            ODArray = self.ODArrayCroped
        else:
            ODArray = self.ODArray
        if index is None:
            # if index is None, we show all the images in a for loop (so we
            # have to call py.show() at every iteration)
            for OD, var in zip(ODArray, self.variableArray):
                py.figure()
                py.title(self.variable + " : " + str(var))
                py.imshow(OD, cmap='jet')
                py.xlabel("momentum")
                py.title(self.variable + " : " + str(var))
                if coords is not None:
                    # if coords is not None, we show the order delimitation
                    for x0x1 in coords:
                        x0, x1 = x0x1
                        py.axvline(x=x0, color='red', linestyle='--')
                        py.axvline(x=x1, color='red', linestyle='--')
                py.colorbar()
                py.show()
        else:
            # if index is given we show only the wanted image
            py.figure()
            py.title(self.variable + " : " + str(self.variableArray[index]))
            py.imshow(ODArray[index], cmap='jet')
            py.xlabel("momentum")
            if coords is not None:
                for x0x1 in coords:
                    x0, x1 = x0x1
                    py.axvline(x=x0, color='red', linestyle='--')
                    py.axvline(x=x1, color='red', linestyle='--')
            py.colorbar()
            py.show()

    def plotAllODAtOnce(self, coords=None, save=True):
        """
        Argument :
            - coords : 1D array. Coordinates of each orders.
            - save : boolean. Decision to save the resulting image.
        Return :
            - None.

        Show all optical density images in one image by
        concatenating them."""
        fig = py.figure()
        all = py.concatenate(self.ODArrayCroped, axis=0)
        py.imshow(all, cmap='jet', aspect=(2 / 10))
        if coords is not None:
            self.setOrdes(coords)
            for x0x1 in self.coords:
                x0, x1 = x0x1
                py.axvline(x=x0, color='red', linestyle='--')
                py.axvline(x=x1, color='red', linestyle='--')
        py.xlabel("momentum")
        py.ylabel(self.variable)
        py.colorbar()
        py.show()
        if save is True:
            fig.savefig('allODAtOnce.svg', format='svg')

    def computeAllProfile(self):
        """
        Argument :
            - None.
        Return :
            - None.

        Fill the "profileArray" attribute with the integral over the vertical
        direction of each OD image in the "ODArrayCroped" attribute.

        _Note_ : the "computeAllOD" method has to be called fist.
        """
        # the profile array is initialized at each call, otherwise it would
        # keep accumulate previous versions.
        self.profileArray = []
        for OD in self.ODArrayCroped:
            profile = py.sum(OD, axis=0)
            # the profile is normalize by its integral
            profile /= py.sum(profile)
            self.profileArray.append(profile)

    def plotProfile(self, index=None, coords=None):
        """
        Argument :
            - index : integer. Index of the profile to plot. If None is given
            all the profiles are ploted.
            - coords : 1D array. Coordinates of each orders.
        Return :
            - None.
        """
        if index is None:
            # if no index is given, all the profiles are plotted ...
            for profile, var in zip(self.profileArray, self.variableArray):
                py.figure()
                py.title(self.variable + " : " + str(var))
                py.xlabel("momentum")
                py.ylabel("Optical denisty")
                py.xticks(range(len(self.variableArray)), self.variableArray)
                py.plot(profile, 'b')
                if coords is not None:
                    for x0x1 in coords:
                        x0, x1 = x0x1
                        py.axvline(x=x0, color='red', linestyle='--')
                        py.axvline(x=x1, color='red', linestyle='--')
                py.show()
        else:
            # ... otherwise, the profile with the given index is plotted.
            py.figure()
            py.title(self.variable + " : " + str(self.variableArray[index]))
            py.xlabel("momentum")
            py.ylabel("Optical denisty")
            py.plot(self.profileArray[index], 'b')
            if coords is not None:
                for x0x1 in coords:
                    x0, x1 = x0x1
                    py.axvline(x=x0, color='red', linestyle='--')
                    py.axvline(x=x1, color='red', linestyle='--')
            py.show()

    def plotAllProfileAtOnce(self, coords=None):
        """
        Argument :
            - coords : 1D array. Coordinates of each orders.
        Return :
            - None.

        Shows all profiles in one image."""

        nPlots = len(self.profileArray)
        fig, axes = py.subplots(nrows=nPlots, ncols=1, sharex=True)
        fig.text(0.5, 0.04, 'common X', ha='center')
        fig.text(0.04, 0.5, 'common Y', va='center', rotation='vertical')
        for ax, i in zip(axes, range(nPlots)):
            ax.plot(self.profileArray[i])
            if coords is not None:
                for x0x1 in coords:
                    x0, x1 = x0x1
                    py.axvline(x=x0, color='red', linestyle='--')
                    py.axvline(x=x1, color='red', linestyle='--')
        py.xlabel("momentum")
        py.ylabel("Optical denisty")
        py.xticks(range(len(self.variableArray)), self.variableArray)
        py.show()

    def computeEvolutionOfOrder(self):
        """
        Argument :
            - None.
        Return :
            - None.

        Evaluate the population of each order for a given optical density
        image."""

        n = len(self.coords)
        indexesOfOrder = [i - int(n / 2) for i in range(n)]
        self.orders = {i: [] for i in indexesOfOrder}
        # Compute the integral of each order
        normalizeFactor = 0
        for profile in self.profileArray:
            for x0x1, i in zip(self.coords, indexesOfOrder):
                x0, x1 = x0x1[0], x0x1[1]
                integral = py.sum(profile[x0:x1])
                self.orders[i].append(integral)
                normalizeFactor += integral
        # Normalize all the evolution by the total number of atoms
        for i in indexesOfOrder:
            self.orders[i] /= normalizeFactor

    def plotOrder(self, index):
        """
        Argument :
            - index : integer. Index of the order to plot.
        Return :
            - None.

        Used to plot the evolution of a specific order.
        """
        py.figure()
        py.title("Order " + str(index))
        py.plot(self.orders[index], 'b+--')
        py.xlabel(self.variable)
        py.ylabel("Normalized denisty")
        py.xticks(range(len(self.variableArray)), self.variableArray)
        py.show()

    def plotAllOrderAtOnce(self, save=True):
        """
        Argument :
            - save : boolean. Decision to save the resulting image.
        Return :
            - None.

        Used to plot the evolution of all order.
        """
        fig = py.figure()
        # list of axis, index of the axis, number of orders
        ax, axis, n = [], 0, len(self.coords)
        # number of rows of the subplot. The subplot has allways
        nrows = int(n / 4) + 1 * (n % 4 != 0)
        for index, order in self.orders.items():
            if axis is 0:
                ax.append(py.subplot(nrows, 4, axis + 1))
            else:
                ax.append(py.subplot(nrows, 4, axis + 1,
                                     sharex=ax[0], sharey=ax[0]))
            ax[axis].set_py.title('Order ' + str(index))
            py.plot(order, 'b+--')
            axis += 1
        py.xticks(range(len(self.variableArray)), self.variableArray)
        py.tight_layout()
        py.show()
        if save is True:
            fig.savefig('allOrderAtOnce.svg', format='svg')
