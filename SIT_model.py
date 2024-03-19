# -*- coding: utf-8 -*-

from Pnem import Pnem

import numpy as np
import numpy.random as random
from matplotlib.pyplot import subplots

class Pnem4(Pnem):
    def __init__(self, m, betaE, nuE, deltaE, K, nu, deltaM, gammaS, deltaS, deltaF, sterile=None, R=False):
        self.m = m
        if type(sterile) == int:
            m[3] = sterile
        self.betaE = betaE
        self.deltaE = deltaE
        self.deltaM = deltaM
        self.deltaF = deltaF
        self.deltaS = deltaS
        self.nuE = nuE
        self.nu = nu
        self.K = K
        self.gammaS = gammaS
        self.tauxN = lambda x: nuE*x[0]
        self.tauxM = lambda x: self.muE*x[0] + self.muM*x[1] + self.muF*x[2] + self.muS*x[3]
        self.sterile = sterile
        self.R = R

        self.N = [self.m.copy()]
        self.t = [0]
        self.actualT = 0
        
        self.totalSterile = self.m[3]

    def getN(self, t):
        if t > self.t[-1]:
            print("t n'est pas atteint")
            return
        elif t < 0:
            print("t est inferieur a 0")
        else:
            for k in range(len(self.t)):
                if t <= self.t[k]:
                    return self.N[k]

    def getNinterval(self, tmin, tmax):
        if tmax > self.t[-1]:
            print("tmax trop grand")
        elif tmin < 0:
            print("tmin trop petit")
        else:
            index1 = 0
            for k in range(len(self.t)):
                if tmin <= self.t[k]:
                    index1 = k
                    break
            for k in range(index1, len(self.t)):
                if tmax <= self.t[k]:
                    return self.N[index1:k]

    def getT(self, condition):
        for k in range(len(self.t)):
            if condition(self.N[k], self.t[k]):
                return self.t[k]
        print("condition jamais verifiee")

    def getAllT(self, condition):
        tlist = []
        for k in range(len(self.t)):
            if condition(self.N[k], self.t[k]):
                tlist.append(self.t[k])
        return tlist

    def simulation(self, tmax, untilExtinction = False):
        if type(self.sterile) == int:
            ttot = self.t[-1]
            while ttot <= tmax:
                nbrInd = self.N[-1]
                if sum(nbrInd) <= 0:
                    break
                else:
                    param = sum(self.getPs())
                    tpsArret = random.exponential(1/param)
                    ttot = ttot + tpsArret
                    self.t.append(ttot)
                    self.event_constant(ttot)
            self.actualT = tmax
        else:
            sterile = self.sterile.copy()
            ttot = self.t[-1]
            while ttot <= tmax:
                nbrInd = self.N[-1]
                if sum(nbrInd) == 0 or (untilExtinction and nbrInd[0] <= 0):
                    break
                else:
                    param = sum(self.getPs())
                    if param <= 0:
                        print(self.N[-1])
                        print(self.t[-1])
                    tpsArret = random.exponential(1/param)
                    ttot = ttot + tpsArret
                    self.t.append(ttot)
                    if len(sterile) > 0 and sterile[0][0] <= ttot:
                        if self.R:
                            self.event_non_constant(ttot, sterile[0][1]*nbrInd[1])
                            self.totalSterile += self.sterile[0][1]*nbrInd[1]
                        else:
                            self.event_non_constant(ttot, sterile[0][1])
                            self.totalSterile += self.sterile[0][1]
                        sterile.pop(0)
                    else:
                        self.event_non_constant(ttot)
            self.actualT = tmax

    def getPs(self):
        E, M, F, Ms = self.N[-1]
        betaE = self.betaE
        muE = self.deltaE
        muM = self.deltaM
        muF = self.deltaF
        muS = self.deltaS
        lbdE = self.nuE
        lbd = self.nu
        K = self.K
        gammaS = self.gammaS

        p1 = betaE*(1-E/K)*F
        p2 = (1-lbd)*lbdE*E
        p3 = lbd*lbdE*M/(M + gammaS*Ms)*E if (M + gammaS*Ms) != 0 else 0
        p4 = (muE + lbd*lbdE*gammaS*Ms/(M + gammaS*Ms))*E if (M + gammaS*Ms) != 0 else 0
        p5 = 0
        p6 = muM*M
        p7 = muF*F
        p8 = muS*Ms
        return p1, p2, p3, p4, p5, p6, p7, p8

    def event_non_constant(self, t, addSterile=0):
        E, M, F, Ms = self.N[-1]
        p1, p2, p3, p4, p5, p6, p7, p8 = self.getPs()
        psum = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
        # if psum == 0:
        #     self.N.append([0, 0, 0, 0])
        #     return
        Ms += addSterile

        tirage = random.random()
        if tirage <= p1/psum:
            self.N.append([E+1, M, F, Ms])
        elif tirage <= (p1 + p2)/psum:
            self.N.append([E-1, M+1, F, Ms])
        elif tirage <= (p1 + p2 + p3)/psum:
            self.N.append([E-1, M, F+1, Ms])
        elif tirage <= (p1 + p2 + p3 + p4)/psum:
            self.N.append([E-1, M, F, Ms])
        elif tirage <= (p1 + p2 + p3 + p4 + p5)/psum:
            self.N.append([E-1, M, F, Ms])
        elif tirage <= (p1 + p2 + p3 + p4 + p5 + p6)/psum:
            self.N.append([E, M-1, F, Ms])
        elif tirage <= (p1 + p2 + p3 + p4 + p5 + p6 + p7)/psum:
            self.N.append([E, M, F-1, Ms])
        else:
            self.N.append([E, M, F, Ms-1])

    def event_constant(self, t):
        E, M, F, Ms = self.N[-1]
        p1, p2, p3, p4, p5, p6, p7, p8 = self.getPs()
        psum = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8

        tirage = random.random()
        if tirage <= p1/psum:
            self.N.append([E+1, M, F, Ms])
        elif tirage <= (p1 + p2)/psum:
            self.N.append([E-1, M+1, F, Ms])
        elif tirage <= (p1 + p2 + p3)/psum:
            self.N.append([E-1, M, F+1, Ms])
        elif tirage <= (p1 + p2 + p3 + p4)/psum:
            self.N.append([E-1, M, F, Ms])
        elif tirage <= (p1 + p2 + p3 + p4 + p5)/psum:
            self.N.append([E-1, M, F, Ms])
        elif tirage <= (p1 + p2 + p3 + p4 + p5 + p6)/psum:
            self.N.append([E, M-1, F, Ms])
        elif tirage <= (p1 + p2 + p3 + p4 + p5 + p6 + p7)/psum:
            self.N.append([E, M, F-1, Ms])
        else:
            self.N.append([E, M, F, Ms])

    def plotp(self, number = [0, 1, 2, 3]):
        
        label=["E", "M", "F", "Ms"]
        color=["dodgerblue", "darkorange", "green", "tomato"]
        
        index = np.searchsorted(self.t, self.actualT)
        x = self.t[:index]
        y = np.array(self.N[:index])
        
        # x = np.append(x, self.actualT)
        # y = np.append(y, y[-1])
        
        fig, host = subplots()
        msplot = host.twinx()
        host.set_xlim(0, self.actualT)
        host.set_xlabel("time")
        host.set_ylabel("quantity for E, M and F")
        msplot.set_ylabel("quantity for Ms")
        
        
        for num in number:
            if num == 3:
                msplot.plot(x, y[:, 3], label=label[num], color=color[num])
            else:
                host.plot(x, y[:, num], label=label[num], color=color[num])
        # ylim(0, 1.3*max(Ebarre, Mbarre, Fbarre))
        # xlabel("Time")
        # ylabel("Population size")
        
        host.set_ylim(bottom=0)
        msplot.set_ylim(bottom=0)
        
        lines_host, labels_host = host.get_legend_handles_labels()
        lines_msplot, labels_msplot = msplot.get_legend_handles_labels()
        fig.legend(lines_host + lines_msplot, labels_host + labels_msplot, loc='upper center')
        
        # fig.tight_layout()
        # fig.savefig(r'C:\Users\talie\OneDrive\Documents\cours\PROJET M2\images\blablabla.png')
        
        
    def esp_temp_retour(self, xarret):
        precision = 200
        esperance = 0
        if type(self.sterile) == int:
            print("erreur, la quantite de moustiques sterile ne peut pas etre constante")
        else:
            self.t = [0]
            esperance = 0
            for _ in range(precision):
                self.t = [0]
                self.N = [self.m.copy()]
                sterile = self.sterile.copy()
                ttot = self.t[-1]
                while self.N[-1][2] >= xarret:
                    param = sum(self.getPs())
                    tpsArret = random.exponential(1/param)
                    ttot = ttot + tpsArret
                    self.t.append(ttot)
                    if len(sterile) > 0 and sterile[0][0] <= ttot:
                        self.event_non_constant(ttot, sterile[0][1])
                        sterile.pop(0)
                    else:
                        self.event_non_constant(ttot)
                    # if(self.N[-1][3] <= 0):
                    #     print("break")
                    #     print(t[-1])
                    #     break
                    # if self.N[-1][2] <= xarret or t[-1] > 40:
                    #     pass
                        # print(self.N[-1][2], t[-1])
                while self.N[-1][2] < xarret:
                    param = sum(self.getPs())
                    tpsArret = random.exponential(1/param)
                    ttot = ttot + tpsArret
                    self.t.append(ttot)
                    if len(sterile) > 0 and sterile[0][0] <= ttot:
                        self.event_non_constant(ttot, sterile[0][1])
                        sterile.pop(0)
                    else:
                        self.event_non_constant(ttot)
                esperance = esperance + self.t[-1]
            esperance = esperance/precision
        return esperance

    def minFemale(self, dropInterval, dropQuantity, dropNumber):
        if type(self.sterile) == int:
            print("erreur, la quantite de moustiques sterile ne peut pas etre constante")
        else:
            self.t = [0]
            self.N = [self.m.copy()]
            ttot = self.t[-1]
            sterile = [(k*dropInterval, dropQuantity) for k in range(dropNumber)]
            minFtab = []
            minFt = []
            S = dropQuantity + self.N[0][3]
            while ttot < dropInterval*dropNumber and len(sterile) > 0:
                minF = self.N[0][2]
                mint = ttot
                F = minF
                minSterile = S
                k = 0
                while S <= minSterile and (len(sterile) > 0 or S > 0):
                    k += 1
                    param = sum(self.getPs())
                    tpsArret = random.exponential(1/param)
                    ttot = ttot + tpsArret
                    self.t.append(ttot)
                    if len(sterile) > 0 and sterile[0][0] <= ttot:
                        if self.R:
                            self.event_non_constant(ttot, sterile[0][1]*self.N[-1][1])
                        else:
                            self.event_non_constant(ttot, sterile[0][1])
                        sterile.pop(0)
                    else:
                        self.event_non_constant(ttot)
                    F = self.N[-1][2]
                    S = self.N[-1][3]
                    if F < minF:
                        minF = F
                        mint = ttot
                    if S < minSterile:
                        minSterile = S
                minFtab.append(minF)
                minFt.append(mint)
            self.actualT = ttot
            return minFt, minFtab

    def minFemale2(self, dropInterval, dropQuantity, dropNumber, T = 0):
        self.sterile = [(k*dropInterval, dropQuantity) for k in range(dropNumber+1)]
        self.N = [self.m.copy()]
        self.t = [0]
        self.simulation(max(T, dropInterval*dropNumber))

        minFemale = np.zeros(len(self.sterile) - 1)
        minT = np.zeros(len(minFemale))
        for k in range(len(self.sterile) - 1):
            N = np.array(self.getNinterval(self.sterile[k][0], self.sterile[k+1][0]))
            minFemale[k] = min(N[:, 2])
            minT[k] = self.getT(lambda N, t: N[2] == minFemale[k] and t >= self.sterile[k][0])
        
        return minT, minFemale

    def plotMinFemale(self, dropInterval, dropQuantity, dropNumber, T = 0, number = [2, 3]):
        dropNumber += 1
        px, py = self.minFemale(dropInterval, dropQuantity, dropNumber)
        if T > dropInterval*dropNumber:
            self.simulation(200)

        # fonction plot

        label=["E", "M", "F", "Ms"]
        color=["dodgerblue", "darkorange", "green", "tomato"]
        
        index = np.searchsorted(self.t, self.actualT)
        x = self.t[:index]
        y = np.array(self.N[:index])
        
        fig, host = subplots()
        msplot = host.twinx()
        host.set_xlim(0, self.actualT)
        host.set_xlabel("time")
        host.set_ylabel("quantity for E, M and F")
        msplot.set_ylabel("quantity for Ms")
        
        
        for num in number:
            if num == 3:
                msplot.plot(x, y[:, 3], label=label[num], color=color[num])
            else:
                host.plot(x, y[:, num], label=label[num], color=color[num])
       
        host.set_ylim(bottom=0)
        msplot.set_ylim(bottom=0)
        
        lines_host, labels_host = host.get_legend_handles_labels()
        lines_msplot, labels_msplot = msplot.get_legend_handles_labels()
        fig.legend(lines_host + lines_msplot, labels_host + labels_msplot, loc='upper center')
        
        # scatter
        
        for k in range(len(px)):
            host.scatter(px[k], py[k], s = 100, color="red", marker="+")
