import matplotlib.pyplot as plt
import numpy as np

class Poisson:

    def __init__(self, nnox=41, nnoy=41, xmax=1, ymax=1, xmin=0, ymin=0):
        self.nnoy = nnoy
        self.nnox = nnox
        self.xmax = xmax
        self.ymax = ymax
        self.xmin = xmin
        self.ymin = ymin

        self.dx, self.dy = (xmax - xmin) / (nnox - 1), (ymax - ymin) / (nnoy - 1)
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

    # Fill matrix and right hand side vector of linear equation system 
    def assembleLinSys(self):

        h1, h2 = self.dx, self.dy
        for i in range(self.nnox, self.nnox * (self.nnoy - 1)):
            self.A[i, i] = -2 * ((1 / (h1 ** 2)) + (1 / (h2 ** 2)))
            self.A[i, i - self.nnox] = 1 / (h2 ** 2)
            self.A[i, i + self.nnox] = 1 / (h2 ** 2)
            self.A[i, i - 1] = 1 / (h1 ** 2)
            self.A[i, i + 1] = 1 / (h1 ** 2)

        for i in range(0, self.nnoy):
            for j in range(0, self.nnox):
                x, y = j*h1, i * h2
                self.b[i * self.nnox + j] = self.f(x, y)

    # Actual (essential) boundary conditions
    def essBC(self, x, y):
        return self.analytic(x, y)

    # Modify matrix and right hand side vector to include essential BC's
    def fillEssBC(self):
        # Boundary conditions for upper and lower boundaries of a rectangle
        for i in range(0, self.nnox):
            for j in range(0, self.nnox * self.nnoy):
                self.A[i, j] = 0.0
                self.A[self.nnox * (self.nnoy - 1) + i, j] = 0.0 
            self.A[i, i] = 1.0
            self.A[((self.nnoy - 1) * self.nnox) + i, (self.nnox * (self.nnoy - 1)) + i] = 1.0
            self.b[i] = self.essBC(i * self.dx, self.ymin)
            self.b[self.nnox * (self.nnoy - 1) + i] = self.essBC(i * self.dx, self.ymax)

        # Boundary conditions for left and right boundaries of a rectangle
        for i in range(0, self.nnoy):
            for j in range(0, self.nnox * self.nnoy):
                self.A[i * self.nnox, j] = 0.0 
                self.A[(self.nnox - 1) + (i * self.nnox), j] = 0.0
            self.A[i * self.nnox, i * self.nnox] = 1.0
            self.A[(self.nnox - 1) + i * self.nnox, (self.nnox - 1) + i * self.nnox] = 1.0
            self.b[i * self.nnox] = self.essBC(self.xmin, i * self.dy)
            self.b[(self.nnox - 1) + (i * self.nnox)] = self.essBC(self.xmax, i * self.dy)
        
    def solve(self):

        self.assembleLinSys()
        self.fillEssBC()
        self.phi = np.linalg.solve(self.A, self.b)

        # Solution along an y center line. For comparison with analytic solution
        for i in range(0, self.nnox):
            self.u[i] = self.phi[int(self.nnoy / 2) * self.nnox + i]
            self.analy[i] = self.analytic(i * self.dx, 0.5)

        # Nodal values of solution
        for i in range(0, self.nnoy):
            for j in range(0, self.nnox):
                self.solu[i, j] = self.phi[i * self.nnox + j]
        self.plot()

    def plot(self):

        x = np.arange(self.xmin, self.xmax, (self.xmax - self.xmin) / self.nnox)
        y = np.arange(self.ymin, self.ymax, (self.ymax - self.ymin) / self.nnoy)
        X, Y = np.meshgrid(x, y)


        fig = plt.figure(figsize=(20, 20))
        fig.suptitle('Plots of solution')

        ax = fig.add_subplot(2, 2, 1, projection='3d')

        ax.plot_surface(X, Y, self.solu, cmap='viridis', edgecolor='none')
        ax.set_title('Surface plot u(x,y)')

        ax = fig.add_subplot(2, 2, 2)
        ax.set_title('2D projection of u(x,y)')
        plt.imshow(self.solu, cmap='hot')
        plt.colorbar()

        ax = fig.add_subplot(2, 2, 3)
        ax.set_title('Fake test vector plot of u(x,y)')
        # ax.set_ylabel('Damped oscillation')
        plt.quiver(X, Y, self.solu, self.solu)

        ax = fig.add_subplot(2, 2, 4)
        ax.set_title('Contour lines for u(x,y)')
        plt.contour(X, Y, self.solu)

        plt.show()

    # Compare nodal values of solution
    def compare2analytic(self):

        l2error = 0.0        
        for i in range(self.nnoy):
            for j in range (self.nnox):
                l2error += (self.solu[i,j] - self.analytic(j * self.dx, i*self.dy))**2 
        l2error = l2error**(0.5)

        print("L2 error is " + str(l2error)) 

# Solve to verify that analytical solution is reproduced
if __name__ == '__main__':

    solver = Poisson(68, 38, 2, 1, -1, -1)
    solver.solve()
    solver.compare2analytic()
