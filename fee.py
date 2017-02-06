import pylab

def fee(dU):
    E = dU/h
    u = (4.39/(Phy**0.5)) - (2.82e7*(Phy**1.5)/E)
    return 1.4e-6*(E**2)*(10**u)/Phy

#Film properties
h   = 20e-7   #[nm]
Phy = 4.28 #[V]

import numpy
x = numpy.linspace(1e-9, 52.0, 1e4)
j = [fee(_x) for _x in x]

with open('20nm.txt', 'w') as file:
    for i in range(len(x)):
        file.write(str(x[i]) + ' ' + str(j[i]) + '\n')
    
pylab.plot(x, j)
pylab.show()
