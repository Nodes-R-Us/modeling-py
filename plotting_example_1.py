import matplotlib.pyplot as plt
from numpy import arange
from math import *


# plt.title(r'$5^{th}$ degree Maclaurin polynomial vs $\sin(x)$')
def taylor_sine(x):
    return(x-(x**3)/(2*3) + (x**5)/(2*3*4*5))

# taylor_label = r'$x- \frac{x^3}{3!}+\frac{x^5}{5!}$'

n = 32 # number of intervals per unit
delta_theta = pi/n

thetas = [-pi+k*delta_theta for k in range(3*n+1)]

taylor_ys = [taylor_sine(t) for t in thetas]

sin_ys = [sin(t) for t in thetas]

# Setting up axes and plotting the graphs

plt.title(r'$5^{th}$ degree Maclaurin polynomial vs $\sin(x)$')
taylor_label = r'$x- \frac{x^3}{3!}+\frac{x^5}{5!}$'
plt.axis([-pi, 2*pi , -2, 2], 'equal')

# draw an x-axis default hline at y=0 that spans the xrange
plt.axhline(y=0, color='black')
# draw a y-axis default vline at x=0 that spans the yrange
plt.axvline(x=0, color='black')

# plotting the taylor poly graph 
plt.plot(thetas, taylor_ys, color='red', label=taylor_label)

# plotting the sine function graph 
plt.plot(thetas, sin_ys, color='blue', label= r'$\sin(x)$')

plt.legend(loc='lower right')

# Show the plot
plt.show()

# This doesn't seem to work?
plt.savefig("taylor_sin.png")


