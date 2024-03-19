# -*- coding: utf-8 -*-
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from SIT_model import Pnem4
from SIT_Euler import PnemEuler
from deterministic import plot_Euler, plot_Euler_constant
from Linear_model import Controle_lineaire

import time as time
import numpy as np



betaE = 10
nuE = 0.05
deltaE = 0.03
K = 1000
nu = 0.49
deltaM = 0.1
gammaS = 1
deltaS = 0.12
deltaF = 0.04

model_params = (betaE, nuE, deltaE, K, nu, deltaM, gammaS, deltaS, deltaF)

Ebarre = K*(1-(deltaF*(nuE+deltaE))/(nu*betaE*nuE))
Mbarre = (1-nu)*nuE*Ebarre/deltaM
Fbarre = nu*nuE*Ebarre/deltaF







if __name__ == "__main__":

    timeBegin = time.time()





    
    # sterile = [(k*5,  40) for k in range(200)]
    
    # # sterile = 20_000
    # T = 350
    # # t, y = Euler(equadiff_function, (e0, m0, f0, ms0), T, 0.01, sterile.copy())

    # # plotEuler(t, y)

    # eee4 = PnemEuler([int(Ebarre), int(Mbarre), int(Fbarre), 10],
    #               *model_params,
    #               sterile = sterile.copy(),
    #               R = True)

    # eee4.simulation(1000)
    # eee4.plotp([0, 1, 2])

    # # eee4.hist_ext(5, 150, precision=3000)
    
    tab = [ 1.5, 2, 5, 10 ]

    for R in tab:
        a = Controle_lineaire(4, 2, R, 1000)
        a.simulation(10)
        a.plotp(False, "R = " + str(R))





    print("temps ecoule : ", time.time() - timeBegin)
