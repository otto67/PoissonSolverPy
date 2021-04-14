import matplotlib.pyplot as plt
import numpy as np

# import mathlib as ml
# from mpl_toolkits import mplot3d

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

class Poisson:

    nnoy, nnox = 41, 41

    def __init__(self, nnoy=41, nnox=41):
        self.nnoy = nnoy
        self.nnox = nnox
        self.dx, self.dy = 1 / (nnox - 1), 1 / (nnoy - 1)
        self.A = np.zeros((nnox * nnoy, nnoy * nnox))
        self.b = np.zeros(nnox * nnoy)
        self.phi = np.zeros((nnoy * nnox))
        self.u = np.zeros(nnox)
        self.analy = np.zeros(nnox)
        self.solu = np.zeros((nnoy, nnox))

    def f(self, x, y):
        return x*x

    def analytic(self, x, y):
        return (x**4/12) + (x/12)


    def assembleLinSys(self):
        h1, h2 = self.dx, self.dy
# Fill matrix and vector. First row and last rows are assumed to be essential boundary conditions
        for i in range(self.nnox, self.nnox * (self.nnoy - 1)):
            self.A[i, i] = -2 * ((1 / (h1 ** 2)) + (1 / (h2 ** 2)))
            self.A[i, i - self.nnox] = 1 / (h2 ** 2)
            self.A[i, i + self.nnox] = 1 / (h2 ** 2)
            self.A[i, i - 1] = 1 / (h1 ** 2)
            self.A[i, i + 1] = 1 / (h1 ** 2)
        for i in range(0, self.nnoy):
            for j in range(0, self.nnox):
                tmp1, tmp2 = j*h1, i * h2
                self.b[i * self.nnox + j] = self.f(tmp1, tmp2)


    def essBC(self, x, y):
        return self.analytic(x, y)


    def fillEssBC(self):
        # Boundary conditions for upper and lower boundaries of a unit square
        for i in range(0, self.nnoy):
            for j in range(0, self.nnox * self.nnoy):
                self.A[i, j] = 0.0
                self.A[self.nnoy * (self.nnox - 1) + i, j] = 0.0
            self.A[i, i] = 1.0
            self.A[self.nnoy * (self.nnox - 1) + i, self.nnoy * (self.nnox - 1) + i] = 1.0
            self.b[i] = self.essBC(i * self.dx, 0)
            self.b[self.nnoy * (self.nnox - 1) + i] = self.essBC(i * self.dx, 1)

        # Boundary conditions for left and right boundaries of A unit square
        for i in range(0, self.nnoy):
            for j in range(0, self.nnox * self.nnoy):
                self.A[i * self.nnox, j] = 0.0
                self.A[(self.nnoy - 1) + (i * self.nnoy), j] = 0.0
            self.A[i * self.nnoy, i * self.nnoy] = 1.0
            self.A[(self.nnoy - 1) + i * self.nnoy, (self.nnoy - 1) + i * self.nnoy] = 1.0
            self.b[i * self.nnoy] = self.essBC(0, i * self.dy)
            self.b[(self.nnoy - 1) + (i * self.nnoy)] = self.essBC(1, i * self.dy)

    def solve(self):

        self.assembleLinSys()
        self.fillEssBC()
        self.phi = np.linalg.solve(self.A, self.b)

        for i in range(0, self.nnox):
            self.u[i] = self.phi[int(self.nnoy / 2) * self.nnox + i]
            self.analy[i] = self.analytic(i * self.dx, 0.5)

        for i in range(0, self.nnoy):
            for j in range(0, self.nnox):
                self.solu[i, j] = self.phi[i * self.nnox + j]
        self.plot()

    def plot(self):

        x = y = np.arange(0.0, 1.0, 1 / self.nnoy)
        X, Y = np.meshgrid(x, y)


        fig = plt.figure(figsize=(20, 20))
        fig.suptitle('A tale of 4 subplots')

        ax = fig.add_subplot(2, 2, 1, projection='3d')

        ax.plot_surface(X, Y, self.solu, cmap='viridis', edgecolor='none')
        ax.set_title('Surface plot u(x,y)')

        ax = fig.add_subplot(2, 2, 2)
        ax.set_title('2D projection of u(x,y)')
        plt.imshow(self.solu, cmap='hot')
        plt.colorbar()

        ax = fig.add_subplot(2, 2, 3)
        ax.set_title('Fake test vector plot of u(x,y)')
        ax.set_ylabel('Damped oscillation')
        plt.quiver(X, Y, self.solu, self.solu)

        ax = fig.add_subplot(2, 2, 4)
        ax.set_title('Contour lines for u(x,y)')
        plt.contour(X, Y, self.solu)

        plt.show()


# Solve to verify that analytical solution is reproduced
if __name__ == '__main__':

    solver = Poisson(41, 41)
    solver.solve()
