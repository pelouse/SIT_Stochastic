# -*- coding: utf-8 -*-

from Pnem import Pnem

import numpy as np
from matplotlib.pyplot import plot, legend, xlim, ylim


class Controle_lineaire(Pnem):

    def __init__(self, lbd, mu, R, Ninit):
        self.lbd = lbd
        self.mu = mu
        self.R = R
        super().__init__(tauxN=lambda x: x*lbd/(1+R),
                       tauxM=lambda x: x*mu,
                       Ninit=Ninit)

    def plotp(self, theorique=True, label=None):
        super().plotp(label)
        if theorique:
            lbdr = self.lbd/(1+self.R)
            precision = 0.1
            x = np.arange(0, self.t[-1] + precision, precision)
            def esperance(t):
                return self.Ninit*np.exp((lbdr - self.mu)*t)
            plot(x, esperance(x))

        legend()
        xlim(0,self.t[-1])
        ylim(0,1.2*max(self.N))


