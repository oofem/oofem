from __future__ import print_function
import liboofem
import numpy as np
# instantiate uniform field
f=liboofem.UniformGridField()

# set 2d geometry
f.setGeometry(lo=(0,0),hi=(1,1),div=(2,2))
# set data: div is 2x2, hence 3x3=9 node grid with 9 values
f.setValues([0,.5,0, .5,1,.5, 0,.5,.5])
# grid which is larger and denser than field domain
# closest-point values are returned for out-of-domain points
X,Y=np.meshgrid(np.arange(-.5,1.5,.1),np.arange(-.5,1.5,.1))
@np.vectorize
def fEval(x,y):
    return f.evaluateAtPos((x,y))[0]
# evaluate field in all meshgrid points
Z=fEval(X,Y)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, linewidth=1, rstride=1, cstride=1,cmap=cm.coolwarm, shade=True)
plt.show()
