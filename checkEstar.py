import numpy as np

cal2joule = 4.184
kB = 0.00198720425
ngrid = 100
T = 87
beta = 1/(kB*T)

if __name__ == '__main__':

    Energy = np.loadtxt("Energy.dat") / cal2joule
    Edim = len(Energy)
    Emin = np.min(Energy)
    Emax = np.max(Energy)
    Func  = np.exp(beta*(Energy-Emax))
    Func2 = np.exp(2*beta*(Energy-Emax))

    def calcErr2(Estar):
        Heaviside = (Energy <= Estar)
        Avg = np.mean(Heaviside*Func)
        Avg2 = np.mean(Heaviside*Func2)
        Err2 = (Avg2/Avg**2 - 1.0)/Edim
        return Err2

    def shiftEstar(Estar):
        Heaviside = (Energy <= Estar)
        Avg = np.mean(Heaviside*Func)
        Avg2 = np.mean(Heaviside*Func2)
        newEstar = (np.log(2) + np.log(Avg2) - np.log(Avg))/beta + Emax
        return newEstar

    def findEstar(Estar0):
        Estar = Estar0
        while True:
            Estar_prev = Estar
            Estar = shiftEstar(Estar)
            dEstar = abs(Estar - Estar_prev)
            print("               Estar  = %20.12f  Err2 = %16.8f" % (Estar,calcErr2(Estar)))
            if dEstar <= 1e-6:  return Estar

    def cutoff(Estar):
        return np.mean(Energy > Estar)*100

    print("GRID SEARCH")
    dE = (Emax-Emin)/(ngrid-1)
    bestEstar = Emax
    bestErr2 = 1e99
    for i in range(ngrid):
        Estar = Emin + i*dE
        if i==ngrid-1: Estar = Emax
        Err2 = calcErr2(Estar)
        print("      Estar = %20.12f  Err2 = %16.8f  cut = %7.3f" % (Estar,Err2,cutoff(Estar)))
        if Err2 < bestErr2:
            bestEstar = Estar
            bestErr2 = Err2
    print("BEST: Estar = %20.12f  Err2 = %16.8f  cut = %7.3f" % (bestEstar,bestErr2,cutoff(bestEstar)))
    print()

    print("MINIMIZATION")
    for ratio in np.linspace(0,1,11):
        Estar0 = (1.0-ratio)*Emin + ratio*Emax
        print("  ratio = %3.1f  Estar0 = %20.12f" % (ratio, Estar0))
        Estar = findEstar(Estar0)
        print("               Estar  = %20.12f  Err2 = %16.8f  cut = %7.3f" % (Estar,calcErr2(Estar),cutoff(Estar)))
    print()
