from math import exp, fabs
import scipy.constants as const
import scipy.integrate as _int

import os

class CVC():
    def __init__(self, plasma, gun, probe, lay, order = 4):
        self.plasma = plasma #Plasma parameters: Te, Ti, n, mi
        self.gun    = gun    #Electron gun: I, U
        self.probe  = probe  #Probe characteristics: S
        self.lay    = lay    #Lay characteristics: d and tunneling probability - D
        self.order  = order
        self.data   = self._CVC()

    def I(self, U):
        return {'I1': self._I(U, self.data['CVC1']), 'I2': self._I(U, self.data['CVC2'])}

    def _I(self, U, CVC):
        for i in range(len(CVC['U'])):
            if (CVC['U'][i] > U):
                x1 = CVC['U'][i-1]
                y1 = CVC['I'][i-1]
                x2 = CVC['U'][i]
                y2 = CVC['I'][i]
                k = (y2 - y1)/(x2 - x1)
                b = (y1*x2 - y2*x1)/(x2 - x1)
                return (k*U + b)

    def _CVC(self):
        CVC1  = {'U': [], 'I': []}
        CVC2  = {'U': [], 'I': []}
        Error = []
        _U = 0
        dU = 0.01
        while (_U < self.gun['U']/5):
            print ('%.2f' % (100*_U/(2*self.gun['U'])) + ' %')
            CVC1['U'].append(_U)
            CVC1['I'].append(self.plasmaCurrent(_U))
            cvc = self.totalCurrent(_U, CVC1['I'][-1])
            CVC2['U'].append(cvc['U'])
            CVC2['I'].append(cvc['I'])
            Error.append(cvc['error'])
            _U += dU
            os.system('cls' if os.name == 'nt' else 'clear')
        return {'CVC1': CVC1, 'CVC2': CVC2, 'error': sum(Error)/len(Error)}

    #Total current
    def totalCurrent(self, U, I):
        dU = self.intersection(self.plasmaCurrent(U), self.order)
        I_FE = self.FECurrent(dU)
        return {'U': U - dU, 'I': self.lay['D']*I_FE, 'error': I - (1 - self.lay['D'])*I_FE}

    #Current from plasma
    def plasmaCurrent(self, U):
        return self.plasmaElectronCurrent(U) - self.ionCurrent(U) + self.hotElectronCurrent(U)

    #Ion current and secondary ion electron emission current

    #Ion saturation current (Ii = 0.52*S*e*n*sqrt(kTe/M))
    def ionCurrent(self, U):
        Vi = ((8*const.k*self.plasma['Ti'])/(const.pi*self.plasma['mi']))**0.5
        Ii = 0.25*const.e*self.plasma['n']*Vi*self.probe['S']
        return Ii*(1+self.gamma(U))
    def gamma(self, E):
        return 0.3*(0.001*E)**0.53

    #Electron current (Maxwellian part)
    def plasmaElectronCurrent(self, U):
        Ve = ((8*const.k*self.plasma['Te'])/(const.pi*const.m_e))**0.5
        return 0.25*const.e*self.plasma['n']*Ve*self.probe['S']*exp(-(const.e*U)/(const.k*self.plasma['Te']))

    #Electron current from the gun and secondary electron-electron emission current
    def hotElectronCurrent(self, U):
        if (U < self.gun['U']):
            self._U = U
            return _int.quad(self._hotElectronCurrent, U, self.gun['U'])[0]
        else:
            return 0
    def _hotElectronCurrent(self, E):
        return self.gun['I']*self.electronGunDistribution(E)*(1 - self.sigma(E - self._U))
    def electronGunDistribution(self, E):
        return 1/self.gun['U']
    def sigma(self, E):
        _sigmaMax  = 6.4
        _eqPower   = 1.2
        _maxEnergy = 600
        return _sigmaMax*(E/_maxEnergy)**0.5*exp(-(E**_eqPower - _maxEnergy**_eqPower)/(2*_eqPower*_maxEnergy**_eqPower))

    #Field emission current in the thin lay of Al2O3 (work function - Phi = 4.28)
    def FECurrent(self, dU):
        if (dU == 0):
            return 0
        else:
            Phi = 4.28
            E = fabs(dU)/(self.lay['d']*1e2)
            u = (4.39/Phi**0.5) - 2.82e7*(Phi**1.5/E)
            j = 1.4e-6*(E**2/Phi)*10**u #[A/cm^-3]
            return 1e4*self.probe['S']*j*(dU/fabs(dU))
    def intersection(self, I, order = 4):
        dU = 0
        step = 1
        for i in range(order):
            while (fabs(I) > (1 - self.lay['D'])*self.FECurrent(dU)):
                dU += step
            dU -= step
            step /= 10**i
        if (I >= 0):
            return (dU-step)
        else:
            return -(dU-step)
