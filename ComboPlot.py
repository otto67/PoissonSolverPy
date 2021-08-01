# Module to plot a vector 

import matplotlib.pyplot as plt
import numpy as np

def plot(solu, delta):

        x = y = np.arange(0.0, 1.0, delta)
        X, Y = np.meshgrid(x, y)

        fig = plt.figure(figsize=(20, 20))
        fig.suptitle("Simulation results \n Note that FEM solver only works for rhs=0")

        ax = fig.add_subplot(2, 2, 1, projection='3d')

        ax.plot_surface(X, Y, solu, cmap='viridis', edgecolor='none')
        ax.set_title('Surface plot u(x,y)')

        ax = fig.add_subplot(2, 2, 2)
        ax.set_title('2D projection of u(x,y)')
        plt.imshow(solu, cmap='hot')
        plt.colorbar()

        ax = fig.add_subplot(2, 2, 3)
        ax.set_title('Fake test vector plot of u(x,y)')
        ax.set_ylabel('Damped oscillation')
        plt.quiver(X, Y, solu, solu)

        ax = fig.add_subplot(2, 2, 4)
        ax.set_title('Contour lines for u(x,y)')
        plt.contour(X, Y, solu)

        plt.show()