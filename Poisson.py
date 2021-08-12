import matplotlib.pyplot as plt
import numpy as np
import comboplot as plotter
from RHS import RHS
from BC import BC

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
        self.analy = np.zeros(nnox)
        self.solu = np.zeros((nnoy, nnox))

        self.rhs = RHS()
        self.bc = BC()

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
                self.b[i * self.nnox + j] = self.rhs.f(x, y)

    # Modify matrix and right hand side vector to include essential BC's
    def fillEssBC(self):
        # Boundary conditions for upper and lower boundaries of a rectangle
        for i in range(0, self.nnox):
            for j in range(0, self.nnox * self.nnoy):
                self.A[i, j] = 0.0
                self.A[self.nnox * (self.nnoy - 1) + i, j] = 0.0 
            self.A[i, i] = 1.0
            self.A[((self.nnoy - 1) * self.nnox) + i, (self.nnox * (self.nnoy - 1)) + i] = 1.0
            self.b[i] = self.bc.essBC(i * self.dx, self.ymin, 4)
            self.b[self.nnox * (self.nnoy - 1) + i] = self.bc.essBC(i * self.dx, self.ymax, 2)

        # Boundary conditions for left and right boundaries of a rectangle
        for i in range(0, self.nnoy):
            for j in range(0, self.nnox * self.nnoy):
                self.A[i * self.nnox, j] = 0.0 
                self.A[(self.nnox - 1) + (i * self.nnox), j] = 0.0
            self.A[i * self.nnox, i * self.nnox] = 1.0
            self.A[(self.nnox - 1) + i * self.nnox, (self.nnox - 1) + i * self.nnox] = 1.0
            self.b[i * self.nnox] = self.bc.essBC(self.xmin, i * self.dy, 3)
            self.b[(self.nnox - 1) + (i * self.nnox)] = self.bc.essBC(self.xmax, i * self.dy, 1)
        
    def solve(self):

        self.assembleLinSys()
        self.fillEssBC()
        self.phi = np.linalg.solve(self.A, self.b)

        
        # Nodal values of solution
        for i in range(0, self.nnoy):
            for j in range(0, self.nnox):
                self.solu[i, j] = self.phi[i * self.nnox + j]
        
        plotter.plot(self.solu, 1/self.nnox)

    
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

    solver = Poisson(51,51, 1, 1, 0, 0)
    solver.solve()
    solver.compare2analytic()
