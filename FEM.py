import numpy as np
from DiscPrms import *
import matplotlib.pyplot as plt
import datetime

import PySimpleGUI as sg
import os.path

class FEM:

    def __init__(self, params, grid_):
        self.prms = params
        self.grid = grid_

    def fillEssBC(self):
        for i in range(self.grid.nno):
            for j in range(self.grid.nno):
                if self.grid.isBoNode(i, j):
                    info = self.grid.boNodeMap[i][j]
                    self.grid.addBC(i, j, self.essBC(info[0], info[1]))

    def essBC(self, x_val, y_val, bound_o):
        raise NotImplementedError()

    def integrands(self, elem, elmat, elvec):
        raise NotImplementedError()

    def integrandsBound(self, elem, elmat, elvec):
        raise NotImplementedError()

    def numItgOverBoundary(self, i, time, elmat, elmvec):
        elem = self.grid.elems[i]
        for pt in range(1, 3):
            elem.initAtItgPtBound(pt)
            self.integrandsBound(elem, elmat, elmvec)


    def numItgOverElm(self, i, time, elmat, elmvec):
        elem = self.grid.elems[i]
        # Only one integration rule implemented so far. For Gauss rule with M=2, weights are 1
        for pt in range(1, 5):
            elem.initAtItgPt(pt)
            self.integrands(elem, elmat, elmvec)


    def assembleLinEqSys(self, time):

        A_e = np.zeros((self.grid.elems[0].ngf, self.grid.elems[0].nbf))  # Assume equal elements throughout
        b_e = np.zeros(self.grid.elems[0].ngf)

        for elm in range(self.grid.n_elms):
            self.numItgOverElm(elm, time, A_e, b_e)
            self.numItgOverBoundary(elm, time, A_e, b_e)

            # Boundary conditions
            for r in range(self.grid.elems[elm].nbf):
                globdof = self.grid.dofmap[elm][r]
                if globdof in self.grid.essbcnodes.keys():
                    val = self.grid.essbcnodes.get(globdof)
                    b_e -= val*A_e[:,r]
                    A_e[r,:] = 0
                    A_e[:,r] = 0
                    A_e[r,r] = 1
                    b_e[r] = val

            # Assemble into global matrix and vector
            for i in range(self.grid.elems[elm].ngf):
                l2gi = self.grid.loc2glob(i, elm) # // self.grid.nno
                for j in range(self.grid.elems[elm].nbf):
                    l2gj = self.grid.loc2glob(j, elm) # % self.grid.nno
                    self.grid.A[l2gi,l2gj] += A_e[i, j]
                self.grid.b[l2gi] += b_e[i]

            A_e[:, :] = 0
            b_e[:] = 0


    def solve(self, time):
        if (time < 0):
            self.fillEssBC() # Fill ess bc into grid/dof data structures
            self.assembleLinEqSys(time)
            self.grid.phi = np.linalg.solve(self.grid.A, self.grid.b)
            sol = self.grid.interpolSolution()
            self.plot(sol)

    def plot(self, solu):

        x = y = np.arange(0.0, 1.0, 1 / self.prms.nno_y)
        X, Y = np.meshgrid(x, y)


        fig = plt.figure(figsize=(20, 20))
        fig.suptitle('NOTE: There is a bug in the FEM solver')

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


class Simulator(FEM):

    def __init__(self, params, grid_):
        super().__init__(params, grid_)
        self.rhs_ = []
        self.bcs_ = [0,0,0,0]

    def integrands(self, elem, elmat, elvec):
        x_val_loc = elem.coor_at_itg_pt[0]
        y_val_loc = elem.coor_at_itg_pt[1]
        nsd = 2

        for i in range(1, elem.ngf + 1):
            for j in range(1, elem.nbf + 1):
                for k in range(1, nsd+1):
                    elmat[i-1,j-1] += elem.dN(i, k)*elem.dN(j, k)*elem.detJxW
            elvec[i-1] += -self.f(x_val_loc, y_val_loc)*elem.N(i)*elem.detJxW


# Note: This function os not yet implemented,
# implying that du/dn = 0 on boundaries on which 
# essential (Dirichlet) boundary conditions are not prescribed 
    def integrandsBound(self, elem, elmat, elvec):
        x_val_loc = elem.coor_at_itg_pt[0]
        nsd = 2

        for i in range(1, elem.ngf + 1):
            for j in range(1, elem.nbf + 1):
                for k in range(1, nsd+1):
                    pass
                    # elmat[i-1,j-1] += elem.N(i)*elem.dN(j, k)*elem.detJxW

    def essBC(self, x_val, y_val, bo_ind):
        if (bo_ind == 1):
            return 0.0 # x_val**2
        if (bo_ind == 3):
            return 0.0 # x_val**2
        if (bo_ind == 2):
            return 0.0 # x_val**2  # x_val*x_val*x_val
        if (bo_ind == 4):
            return 0.0 # x_val**2  # x_val*x_val*x_val

        print("essBC : Something's wrong")
        return -999999999999

    # Right hand side
    def f(self, x, y):
        return 1



if __name__ == '__main__':

    parameters = DiscPrms(nnx=51, nny=51, dt=-1, t_max=-1)
    grid = Grid2d(parameters)
    boind_list = [1, 3]
    grid.setBoindWithEssBC(boind_list)
    sim = Simulator(parameters, grid)
    sim.solve(-1)