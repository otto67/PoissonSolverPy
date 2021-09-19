# Base class for solving a PDE with FEM
import numpy as np
from DiscPrms import *
import PySimpleGUI as sg
import comboplot as plotter
from RHS import RHS
from BC import BC

# Equation system assembly and solution for a general differential equation using FEM
# This class is not to be instantiated, as it does not implement any integrand for the
# integral
class FEM:

    def __init__(self, params, grid_):
        self.prms = params
        self.grid = grid_
        self.bc = BC()

    # Adds bc's to the dictionary of the grid object
    def fillEssBC(self):
        for i in range(self.prms.nno_x):
            for j in range(self.prms.nno_y):
                if self.grid.isBoNode(i, j):
                    info = self.grid.boNodeMap[i][j]
                    self.grid.addBC(i, j, self.bc.essBC(info[0], info[1], info[2]))

    # Integrand for the integral to be calculated
    def integrands(self, elem, elmat, elvec):
        raise NotImplementedError()

    # Integrand for boundary in case of neuman bc's
    def integrandsBound(self, elem, elmat, elvec):
        raise NotImplementedError()
    
    # Numerical integration over boundary
    # Assembles into local element matrix and vector
    def numItgOverBoundary(self, i, time, elm_mat, elm_vec):
        elem = self.grid.elems[i]
        for pt in range(1, 3):
            elem.initAtItgPtBound(pt)
            self.integrandsBound(elem, elm_mat, elm_vec)

    # Numerical integration over element
    # Assembles into local element matrix and vector
    def numItgOverElm(self, i, time, elmat, elmvec):
        elem = self.grid.elems[i]
        # Only one integration rule implemented so far. 
        # For Gauss rule with M=2, weights are 1
        for pt in range(1, 5):
            elem.initAtItgPt(pt)
            self.integrands(elem, elmat, elmvec)

    # Loops over all element and adds local element contributions to
    # global equation system
    def assembleLinEqSys(self, time):

        # Assume equal elements throughout
        A_e = np.zeros((self.grid.elems[0].ngf, self.grid.elems[0].nbf))  
        b_e = np.zeros(self.grid.elems[0].ngf)

        # integration loop
        for elm in range(self.grid.n_elms):
            self.numItgOverElm(elm, time, A_e, b_e)
        #    self.numItgOverBoundary(elm, time, A_e, b_e)

            # Boundary conditions
            for r in range(self.grid.elems[elm].nbf):
                globdof = self.grid.dofmap[elm][r] 
                if globdof in self.grid.ess_bc_nodes.keys():
                    val = self.grid.ess_bc_nodes.get(globdof)
                    b_e -= val*A_e[:,r]
                    A_e[r,:] = 0.0
                    A_e[:,r] = 0.0
                    A_e[r,r] = 1.0
                    b_e[r] = val
            
            # Assemble into global matrix and vector
            for i in range(self.grid.elems[elm].ngf):
                l2gi = self.grid.loc2glob(i, elm) 
                for j in range(self.grid.elems[elm].nbf):
                    l2gj = self.grid.loc2glob(j, elm) 
                    self.grid.A[l2gi,l2gj] += A_e[i, j]
                self.grid.b[l2gi] += b_e[i]

            A_e[:,:] = 0.0
            b_e[:] = 0.0


    def solve(self, time):
        if (time < 0):
            self.fillEssBC() # Fill ess bc into grid/dof data structures
            self.assembleLinEqSys(time)
            self.grid.phi = np.linalg.solve(self.grid.A, self.grid.b)
            sol = self.grid.interpolSolution()
            plotter.plot(sol, 1/self.prms.nno_x)

# Implements the integrand for a Poisson equation
# This particular subclass is used for veryfying the 
# implementation for an analytical solution

class PoissonFEM(FEM):

    def __init__(self, params, grid_):
        super().__init__(params, grid_)
        self.rhs = RHS()

    def integrands(self, elem, elmat, elvec):
        x_val_loc = elem.coor_at_itg_pt[0]
        y_val_loc = elem.coor_at_itg_pt[1]
        nsd = 2

        for i in range(1, elem.ngf + 1):
            for j in range(1, elem.nbf + 1):
                for k in range(1, nsd+1):
                    elmat[i-1,j-1] += elem.dN(i, k)*elem.dN(j, k)*elem.detJxW
            elvec[i-1] += -self.rhs.f(x_val_loc, y_val_loc)*elem.N(i)*elem.detJxW

# Note: This function is not yet implemented,
# implying that du/dn = 0 on boundaries on which 
# essential (Dirichlet) boundary conditions are not prescribed 
    def integrandsBound(self, elem, elmat, elvec):
        x_val_loc = elem.coor_at_itg_pt[0]
        nsd = 2

        for i in range(1, elem.ngf + 1):
            for j in range(1, elem.nbf + 1):
                for k in range(1, nsd+1):
                    pass # elmat[i-1,j-1] += elem.N(i)*elem.dN(j, k)*elem.detJxW

    # Analytic solution for test
    def analytic(self, x, y):
        return (x**4/12) + (x/12)
    
    # Compares nodal solution to analytic solution for test case 
    # and reports the L2 error
    def compare2analytic(self):

        sol = self.grid.interpolSolution()
        l2error = 0.0        
        for i in range(self.prms.nno_y):
            for j in range (self.prms.nno_x):
                l2error += (sol[i,j] - self.analytic(j * self.prms.dx, i*self.prms.dy))**2 
        l2error = l2error**(0.5)

        print("L2 error is " + str(l2error)) 


if __name__ == '__main__':

    parameters = DiscPrms(nnx=45, nny=45, dt=-1, t_max=-1)
    grid = Grid2d(parameters)
    boind_list = [1, 3]
    grid.setBoindWithEssBC(boind_list)
    sim = PoissonFEM(parameters, grid)
    sim.solve(-1)
    sim.compare2analytic()
