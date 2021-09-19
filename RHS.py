
# Class providing the right hand side of the equation
# for both FDM and FEM
class RHS:
    def __init__(self) -> None:
        pass

    def f(self, x, y):
        return x*x

    def attachRHS(self, rhs_):
        pass

# Polynomial RHS
class polynomialRHS(RHS):
        
    def __init__(self) -> None:
        self.rhs = []
    
    def attachRHS(self, rhs_):
        self.rhs = rhs_

    def f(self, x, y):
        if not self.rhs:
            return 0

        retval = 0.0
        if len(self.rhs) == 1:  # Constant right hand side
            retval = self.rhs[0]
        elif len(self.rhs) == 3: # Linear rhs
            retval = self.rhs[0]*x + self.rhs[1]*y + self.rhs[2]
        elif len(self.rhs) == 6: # Quadratic rhs
            retval = self.rhs[0]*x**2 + self.rhs[1]*y**2 + self.rhs[2]*x*y + \
                     self.rhs[3]*x + self.rhs[4]*y + self.rhs[5]

        return retval
