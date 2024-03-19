# -*- coding: utf-8 -*-

import numpy as np
import numpy.random as random
from matplotlib.pyplot import plot, xlim, ylim, xlabel, ylabel, savefig


class Pnem:

    def __init__(self, tauxN, tauxM, Ninit):
        self.tauxN = np.vectorize(tauxN)
        self.tauxM = np.vectorize(tauxM)
        self.Ninit = Ninit

    def simulation(self, tmax):
        self.N = [self.Ninit]
        self.t = [0]
        ttot = 0

        if type(self.Ninit) == int:
            end = 0
            calculNbrInd = lambda x: x
        else:
            end = [0 for _ in range(len(self.Ninit))]
            calculNbrInd = lambda x: sum(x)
        while ttot <= tmax:
            nbrInd = calculNbrInd(self.N[-1])
            if nbrInd <= 0:
                self.t.append(tmax)
                self.N.append(end)
                break
            else:
                naissance = self.tauxN(nbrInd)
                mort = self.tauxM(nbrInd)
                tpsArret = random.exponential(1/(mort + naissance))
                ttot = ttot + tpsArret
                if ttot <= tmax:
                    self.t.append(ttot)
                    self.event(ttot)

    def event(self, t):
        nbrInd = self.N[-1]
        naissance, mort = self.getNM()
        amplitude = random.random()
        self.N.append(nbrInd -1 + 2*int(amplitude <= (naissance)/(naissance+mort)))

    def getNM(self):
        nbrInd = self.N[-1]
        naissance = self.tauxN(nbrInd)
        mort = self.tauxM(nbrInd)
        return naissance, mort

    def plotp(self, label = None):
        x = [k for t in zip(self.t, self.t) for k in t][1:]
        y = [k for t in zip(self.N, self.N) for k in t][:-1]
        plot(x, y, label=label)
        xlim(0)
        ylim(0)
        xlabel("Time")
        ylabel("Population size")

    def savep(self, label = None, path = "here.png"):
        x = [k for t in zip(self.t, self.t) for k in t][1:]
        y = [k for t in zip(self.N, self.N) for k in t][:-1]
        plot(x, y, label=label)
        xlim(0)
        ylim(0)
        xlabel("Time")
        ylabel("Population size")
        savefig(path)
