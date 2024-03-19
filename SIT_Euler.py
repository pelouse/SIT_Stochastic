# -*- coding: utf-8 -*-

from SIT_model import Pnem4

import numpy as np
import numpy.random as random
from matplotlib.pyplot import subplots, show, hist, title, legend


class PnemEuler(Pnem4):

    def __init__(self, m, betaE, nuE, deltaE, K, nu, deltaM, gammaS, deltaS, deltaF, sterile=None, R=False, sizeMin=500):
        self.sizeMin = sizeMin
        super().__init__(m, betaE, nuE, deltaE, K, nu, deltaM, gammaS, deltaS, deltaF, sterile, R)
        self.totalSterile = self.m[3]
        self.tEuler, self.yEuler = self.Euler(0.01)
        self.t = [self.tEuler[-1]]
        self.N = [[int(self.yEuler[-1][k]) for k in range(4)]]
        print(self.t, self.N)
        

    def equadiff_function(self, E, M, F, Ms, t):
        betaE = self.betaE
        K = self.K
        nuE = self.nuE
        deltaE = self.deltaE
        nu = self.nu
        deltaM = self.deltaM
        gammaS = self.gammaS
        deltaF = self.deltaF
        deltaF =self.deltaF
        deltaS = self.deltaS
        return np.array((betaE*F*(1-E/K) - (nuE + deltaE)*E,
                (1-nu)*nuE*E - deltaM*M,
                nu*nuE*E*M/(M+gammaS*Ms) - deltaF*F,
                -deltaS*Ms))

    def Euler(self, eps):
        t = np.arange(0, self.sterile[-1][0] + eps, eps)
        d = len(t)
        y = np.zeros((d, len(self.m)))
        y[0] = self.m.copy()
        if self.R:
            for k in range(1, d):
                y[k] = y[k-1] + eps*self.equadiff_function(*y[k-1], t[k-1])
                if len(self.sterile) > 0 and t[k] >= self.sterile[0][0]:
                    y[k, 3] += y[k-1][1]*self.sterile[0][1]
                    self.totalSterile += y[k-1][1]*self.sterile[0][1]
                    self.sterile.pop(0)
                if y[k, 0] < self.sizeMin:
                    return t[:k], y[:k]
        else:
            for k in range(1, d):
                y[k] = y[k-1] + eps*self.equadiff_function(*y[k-1], t[k-1])
                if len(self.sterile) > 0 and t[k] >= self.sterile[0][0]:
                    y[k, 3] += self.sterile[0][1]
                    self.sterile[0][1]
                    self.totalSterile += self.sterile[0][1]
                    self.sterile.pop(0)
                if y[k, 0] < self.sizeMin:
                    return t[:k], y[:k]
        return t, y

    def plotp(self, number=[0, 1, 2, 3]):
        label = ["E_es", "M_es", "F_es", "Ms_es"]
        color=["dodgerblue", "darkorange", "green", "tomato"]
        
        fig, host = subplots()
        msplot = host.twinx()
        host.set_xlim(0, max(self.actualT, self.tEuler[-1]))
        host.set_xlabel("time")
        host.set_ylabel("quantity for E, M and F")
        msplot.set_ylabel("quantity for Ms")
        
        if len(self.N) > 1:
        
            index = np.searchsorted(self.t, self.actualT)
            x = self.t[:index]
            y = np.array(self.N[:index])[:, number]
            
            for num in number:
                if num == 3:
                    msplot.plot(x, y[:, 3], label=label[num], color=color[num])
                    msplot.plot(self.tEuler, self.yEuler[:, 3], color=color[num])
                else:
                    host.plot(x, y[:, num], label=label[num], color=color[num])
                    host.plot(self.tEuler, self.yEuler[:, num], color=color[num])
        
        else:
            for num in number:
                if num == 3:
                    msplot.plot(self.tEuler, self.yEuler[:, 3], color=color[num])
                else:
                    host.plot(self.tEuler, self.yEuler[:, num], color=color[num])
            
        
        # self.t = [0]
        # self.N = [self.m.copy()]
        # self.sterile = [(k*5, 100000) for k in range(100)]
        # self.simulation(400)
        # label = ["E_ts", "M_ts", "F_ts", "Ms_ts"]
        # color=["blue", "orange", "violet", "red"]
        # index = np.searchsorted(self.t, self.actualT)
        # x = self.t[:index]
        # y = np.array(self.N[:index])[:, number]
        # for num in number:
        #     if num == 3:
        #         msplot.plot(x, y[:, 3], label=label[num], color=color[num])
        #     else:
        #         host.plot(x, y[:, num], label=label[num], color=color[num])
        
        
        host.set_ylim(bottom=0)
        msplot.set_ylim(bottom=0)
        
        
        lines_host, labels_host = host.get_legend_handles_labels()
        lines_msplot, labels_msplot = msplot.get_legend_handles_labels()
        fig.legend(lines_host + lines_msplot, labels_host + labels_msplot, loc='upper right')
        
        # fig.tight_layout()
        # fig.savefig(r'C:\Users\talie\OneDrive\Documents\cours\PROJET M2\images\blablabla.png')

        show()
        
    def simulation(self, tmax, untilExtinction = False):
        if tmax > self.t[-1]:
            super().simulation(tmax, untilExtinction)


    def extinction(self, dropInterval, dropQuantity, tmax = 1000):
        sterile = [(dropInterval*k, dropQuantity) for k in range(int(tmax/dropInterval)+2)]
        self.sterile = sterile
        self.m = [self.m[0], self.m[1], self.m[2], 0]
        self.totalSterile = 0
        self.tEuler, self.yEuler = self.Euler(0.01)
        self.t = [self.tEuler[-1]]
        self.N = [[int(self.yEuler[-1][k]) for k in range(4)]]
        self.simulation(tmax, untilExtinction = True)
        if self.t[-1] >= tmax:
            return (self.totalSterile, 0)
        else:
            return (self.totalSterile, self.t[-1])

    def tabExtinction(self, tabInterval, tabQuantity, tmax = 1000):
        sizeInterval = len(tabInterval)
        sizeQuantity = len(tabQuantity)
        datas = np.zeros((sizeInterval, sizeQuantity, 2))
        for interval in range(sizeInterval):
            for quantity in range(sizeQuantity):
                datas[interval, quantity] = self.extinction(tabInterval[interval],
                                                            tabQuantity[quantity],
                                                            tmax)
        return datas

    def hist_ext(self, dropInterval, dropQuantity, tmax = 1000, precision = 100):
        steril_list = np.zeros(precision)
        time_list = np.zeros(precision)
        for k in range(precision):
            temp = self.extinction(dropInterval, dropQuantity, tmax)
            steril_list[k] = temp[0]
            time_list[k] = temp[1]
        
        indexes = np.where(time_list > 0)
        
        steril_hist = steril_list[indexes]
        time_hist = time_list[indexes]
        
        proba = np.mean([1 if ext > 0 else 0 for ext in time_list])
        
        hist(steril_hist, bins=50, density=True, label="extinction proba = " + str(proba))
        title("sterile histogram")
        legend()
        show()
        hist(time_hist, bins=50, density=True, label="extinction proba = " + str(proba))
        title("extinction time histogram")
        legend()
        show()
        
    def dichotomy_min_steril(self, dropInterval, dropQuantity, nbrIter = 2, tmax = 1000):
        nbr_interval = 3
        minDropInterval, maxDropInterval = dropInterval
        minDropQuantity, maxDropQuantity = dropQuantity
        r = np.zeros((nbrIter, 3))
        for k in range(nbrIter):
            interval = np.linspace(minDropInterval, maxDropInterval, nbr_interval)
            quantity = np.linspace(minDropQuantity, maxDropQuantity, nbr_interval)
            df = self.tabExtinction(interval, quantity, tmax=tmax)
            minS = 10_000_000
            indexes = (0, 0)
            for i in range(nbr_interval):
                for j in range(nbr_interval):
                    # print(df[i, j, 1], df[i, j, 0])
                    if df[i, j, 1] != 0 and df[i, j, 0] < minS:
                        minS = df[i, j, 0]
                        indexes = (i, j)
            i, j = indexes
            print(indexes)
            if i == 0:
                minDropInterval = interval[i]
                maxDropInterval = interval[i+1]
                if j == 0:
                    minDropQuantity = quantity[j]
                    maxDropQuantity = quantity[j+1]
                elif j == nbr_interval-1:
                    minDropQuantity = quantity[j-1]
                    maxDropQuantity = quantity[j]
                else:
                    if df[i+1, j+1, 1] != 0 and df[i+1, j+1, 0] < df[i+1, j-1, 0]:
                        minDropQuantity = quantity[j]
                        maxDropQuantity = quantity[j+1]
                    else:
                        minDropQuantity = quantity[j-1]
                        maxDropQuantity = quantity[j]
            elif indexes[0] == nbr_interval-1:
                maxDropInterval = interval[i]
                minDropInterval = interval[i-1]
                if j == 0:
                    minDropQuantity = quantity[j]
                    maxDropQuantity = quantity[j+1]
                elif j == nbr_interval-1:
                    minDropQuantity = quantity[j-1]
                    maxDropQuantity = quantity[j]
                else:
                    if df[i-1, j-1, 1] != 0 and df[i-1, j-1, 0] < df[i-1, j+1, 0]:
                        minDropQuantity = quantity[j-1]
                        maxDropQuantity = quantity[j]
                    else:
                        minDropQuantity = quantity[j]
                        maxDropQuantity = quantity[j+1]
            else:
                if j == 0:
                    minDropQuantity = quantity[j]
                    maxDropQuantity = quantity[j+1]
                elif j == nbr_interval-1:
                    minDropQuantity = quantity[j-1]
                    maxDropQuantity = quantity[j]
                else:
                    minCorner = min(df[i-1, j-1, 0]*df[i-1, i-1, 1]!=0,
                                    df[i+1, j-1, 0]*df[i+1, j-1, 1]!=0,
                                    df[i-1, j+1, 0]*df[i-1, j+1, 1]!=0,
                                    df[i+1, j+1, 0]*df[i+1, j+1, 1]!=0)
                    if minCorner == df[i-1, j-1, 0]:
                        minDropInterval = interval[i-1]
                        maxDropInterval = interval[i]
                        minDropQuantity = quantity[j-1]
                        maxDropQuantity = quantity[j]
                    elif minCorner == df[i-1, j+1, 0]:
                        minDropInterval = interval[i-1]
                        maxDropInterval = interval[i]
                        minDropQuantity = quantity[j]
                        maxDropQuantity = quantity[j+1]
                    elif minCorner == df[i+1, j-1, 0]:
                        minDropInterval = interval[i]
                        maxDropInterval = interval[i+1]
                        minDropQuantity = quantity[j-1]
                        maxDropQuantity = quantity[j]
                    elif minCorner == df[i+1, j+1, 0]:
                        minDropInterval = interval[i]
                        maxDropInterval = interval[i+1]
                        minDropQuantity = quantity[j]
                        maxDropQuantity = quantity[j+1]
            r[k] = [minS, interval[i], quantity[j]]
        return r

    def opt_sto(self, dropInterval, dropQuantity, varianceInterval, varianceQuantity, nbrIter, nbrMC, minExtPr, tmax=1000):
        
        def mean(tab):
            if len(tab) == 0:
                return 0
            else:
                return np.mean(tab)

        def calculJ(interval, quantity):
            extinctions = np.zeros(nbrMC)
            values = np.zeros(nbrMC)
            for k in range(nbrMC):
                temp = self.extinction(interval, quantity, tmax)
                values[k] = temp[0]
                extinctions[k] = temp[1]
            proba = mean([1 if ext > 0 else 0 for ext in extinctions])
            value = mean(values[np.where(extinctions > 0)])
            return value, proba
            
        minJ, proba = calculJ(dropInterval, dropQuantity)
        self.plotp(number=[0,1,2])
        minParams = (dropInterval, dropQuantity)
        params = np.zeros((nbrIter, 2))
        params[0] = minParams
        for k in range(nbrIter-1):
            X = random.normal(0, varianceInterval)
            Y = random.normal(0, varianceQuantity)
            newParams = (np.abs(minParams[0] + X), int(np.abs(minParams[1] + Y)))
            print("try params : ", newParams)
            newJ, proba = calculJ(*newParams)
            # print("params : ", minParams, " ---> ", minJ)
            # print("new params : ", newParams, " ---> ", newJ)
            print("proba : ", proba, " value : ", newJ)
            if proba >= minExtPr and newJ < minJ:
                print("accept params : ", newParams)
                varianceInterval = varianceInterval
                varianceQuantity = varianceQuantity
                params[k+1] = (newParams[0], newParams[1])
                minParams = newParams
                minJ = newJ
            else:
                params[k+1] = (minParams[0], minParams[1])
        
        return params

    def dicho(self, intervalArray, quantityArray, nbrMC, minExtPr, tmax=1000):
        def mean(tab):
            if len(tab) == 0:
                return 0
            else:
                return np.mean(tab)

        def calculJ(interval, quantity):
            extinctions = np.zeros(nbrMC)
            values = np.zeros(nbrMC)
            for k in range(nbrMC):
                temp = self.extinction(interval, quantity, tmax)
                values[k] = temp[0]
                extinctions[k] = temp[1]
            proba = mean([1 if ext > 0 else 0 for ext in extinctions])
            value = mean(values[np.where(extinctions > 0)])
            return value, proba
        
        leni, lenq = len(intervalArray), len(quantityArray)
        
        
        minValue = None
        minInterval = None
        minQuantity = None
        for i in range(leni):
            for q in range(lenq):
                v, p = calculJ(intervalArray[i], quantityArray[q])
                print("try for params : ", intervalArray[i], ", ", quantityArray[q])
                print("value : ", v, " proba : ", p)
                if p >= minExtPr and (minValue is None or v < minValue):
                    print("accepted params : ", intervalArray[i], ", ", quantityArray[q])
                    minInterval = intervalArray[i]
                    minQuantity = quantityArray[q]
                    minValue = v
        
        return minInterval, minQuantity
