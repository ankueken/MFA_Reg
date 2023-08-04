from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

def Bspline(y,n):
    #y = [0,1,3,4,3,5,7,5,2,3,4,8,9,8,7] # MID data
    #x = range(0, int(n))
    x = [0,5,10,20,40]

    tck = interpolate.splrep(x, y, k=2)
    x_new = np.linspace(min(x), max(x), int(n))
    spl = interpolate.BSpline(*tck)
    y_fit = spl(x_new)
    return(y_fit)

y_fit = Bspline(y,n)


#fig, ax=plt.subplots(dpi=300)
#plt.title("BSpline curve fitting")
#plt.plot(x, y, 'ro', label="original")
#plt.plot(x_new, y_fit, '-c', label="B-spline")
#plt.legend(loc='best', fancybox=True, shadow=True)
#plt.grid()
#plt.show() 