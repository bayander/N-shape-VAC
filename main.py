import CVC

plasma = {'Te': 2*11600, 'n': 1e17, 'Ti': 300, 'mi': 1.67e-27}
gun    = {'I': 0.1, 'U': 100}
probe  = {'S': 3.14*1.5e-2**2}
lay    = {'d': 1e-9, 'D': 0.8}

CVC = CVC.CVC(plasma, gun, probe, lay, 4)
data = CVC.data

print('main error = ', data['error'])

import pylab
pylab.plot(data['CVC1']['U'], data['CVC1']['I'])
pylab.plot(data['CVC2']['U'], data['CVC2']['I'])

U = 30
I = CVC.I(U)
pylab.plot(U, I['I1'], 'o')
pylab.plot(U, I['I2'], 'o')

minY = min(data['CVC2']['I'])
pylab.axis([0, 2*gun['U'], minY, -minY])
pylab.show()
