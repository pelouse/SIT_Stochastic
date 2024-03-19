# -*- coding: utf-8 -*-

import numpy as np
from matplotlib.pyplot import subplots


def equadiff_function(betaE, nuE, deltaE, K, nu, deltaM, gammaS, deltaS, deltaF, E, M, F, Ms, t):
    return np.array((betaE*F*(1-E/K) - (nuE + deltaE)*E,
            (1-nu)*nuE*E - deltaM*M,
            nu*nuE*E*M/(M+gammaS*Ms) - deltaF*F,
            -deltaS*Ms))


def equadiff_function_constant_sterile(betaE, nuE, deltaE, K, nu, deltaM, gammaS, deltaS, deltaF, E, M, F, Ms, t):
    return np.array((betaE*F*(1-E/K) - (nuE + deltaE)*E,
            (1-nu)*nuE*E - deltaM*M,
            nu*nuE*E*M/(M+gammaS*Ms) - deltaF*F,
            0))


def Euler(f, x0, tmax, eps, sterile):
    t = np.arange(0, tmax + eps, eps)
    d = len(t)
    y = np.zeros((d, len(x0)))
    y[0] = x0
    for k in range(1, d):
        y[k] = y[k-1] + eps*f(*y[k-1], t[k-1])
        if len(sterile) > 0 and t[k] >= sterile[0][0]:
            y[k, 3] += sterile[0][1]
            sterile.pop(0)
    return (t, y)


def Euler_constant(f, x0, tmax, eps):
    t = np.arange(0, tmax + eps, eps)
    d = len(t)
    if type(x0) == int:
        y = np.zeros(d)
        y[0] = x0
        for k in range(1, d):
            y[k] = y[k-1] + eps*f(y[k-1], t[k-1])
    else:
        y = np.zeros((d, len(x0)))
        y[0] = x0
        for k in range(1, d):
            y[k] = y[k-1] + eps*f(*y[k-1], t[k-1])
    return (t, y)



def plot_Euler(params, x0, tmax, sterile, number=[0,1,2,3]):
    t, y = Euler(lambda E, M, F, Ms, t: equadiff_function(*params, E, M, F, Ms, t),
                 x0, tmax, 1, sterile)

    label=["E", "M", "F", "Ms"]
    color=["dodgerblue", "darkorange", "green", "tomato"]

    fig, host = subplots()
    msplot = host.twinx()
    host.set_xlim(0, tmax)
    host.set_xlabel("time")
    host.set_ylabel("quantity for E, M and F")
    msplot.set_ylabel("quantity for Ms")
    
    
    for num in number:
        if num == 3:
            msplot.plot(t, y[:, 3], label=label[num], color=color[num])
        else:
            host.plot(t, y[:, num], label=label[num], color=color[num])
    
    host.set_ylim(bottom=0)
    msplot.set_ylim(bottom=0)
    
    lines_host, labels_host = host.get_legend_handles_labels()
    lines_msplot, labels_msplot = msplot.get_legend_handles_labels()
    fig.legend(lines_host + lines_msplot, labels_host + labels_msplot, loc='upper center')


def plot_Euler_constant(params, x0, tmax, number=[0,1,2,3]):
    t, y = Euler_constant(lambda E, M, F, Ms, t: equadiff_function_constant_sterile(*params, E, M, F, Ms, t),
                 x0, tmax, 1)

    label=["E", "M", "F", "Ms"]
    color=["dodgerblue", "darkorange", "green", "tomato"]

    fig, host = subplots()
    msplot = host.twinx()
    host.set_xlim(0, tmax)
    host.set_xlabel("time")
    host.set_ylabel("quantity for E, M and F")
    msplot.set_ylabel("quantity for Ms")
    
    
    for num in number:
        if num == 3:
            msplot.plot(t, y[:, 3], label=label[num], color=color[num])
        else:
            host.plot(t, y[:, num], label=label[num], color=color[num])
    
    host.set_ylim(bottom=0)
    msplot.set_ylim(bottom=0)
    
    lines_host, labels_host = host.get_legend_handles_labels()
    lines_msplot, labels_msplot = msplot.get_legend_handles_labels()
    fig.legend(lines_host + lines_msplot, labels_host + labels_msplot, loc='upper center')



