import pylab, os

def f(x, X, Y):
    for i in range(len(X)):
        if (X[i] > x):
            x1 = X[i-1]
            y1 = Y[i-1]
            x2 = X[i]
            y2 = Y[i]
            k = (y2 - y1)/(x2 - x1)
            b = (y1*x2 - y2*x1)/(x2 - x1)
            return (k*x + b)

print "Plotting..."

for file in os.listdir(os.getcwd()):
    if file.endswith(".txt"):
        _file = open(file)
        U = []
        I = []
        for line in _file:
            U.append(float(line.split('\t')[0]))
            I.append(float(line.split('\t')[1]))
        _file.close()

        x_step = 1.0
        x_min  = min(U)
        x_max  = max(U)
        x = [(x_min + i*x_step) for i in range(int((x_max - x_min)/x_step))]
        y = [f(_x, U, I) for _x in x]
        pylab.plot(x, y, '.-')
pylab.show()
