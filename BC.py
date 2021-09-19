# Class for handling essential boundary conditions
# for both FDM and FEM

class BC:
    def __init__(self) -> None:
        pass

# Analytic solution for f = x*x
# x,y - positions in domain, floats
# boind - boundary index, integer
    def essBC(self, x, y, boind):
        return (x**4/12) + (x/12)
    
# Boundary conditions prescribed through 
# PySimpleGUI
# bcs_ - List of boundary conditions for boundary indexes
    def attachBC(bcs_):
        pass

# Constant bc's on all boundaries
class constBC(BC):
    def __init__(self) -> None:
        self.bcs = []

    def attachBC(self, bcs_):
        self.bcs = bcs_

    def essBC(self, x, y, boind):

        if not self.bcs:
            return 0.0 

        if boind == 1:
            return self.bcs[0]
        if boind == 2:
            return self.bcs[1]
        if boind == 3:
            return self.bcs[2]
        if boind == 4:
            return self.bcs[3]  
        
        return 0.0
