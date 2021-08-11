from Poisson import Poisson
import DiscPrms as dsc 
import PySimpleGUI as sg
import FEM
from RHS import polynomialRHS
from BC import constBC

# Derived class to implement the right hand side  
# and boundary conditions for FDM
class PoissonSub(Poisson):

    def __init__(self, nno_x, nno_y, xmax, ymax, xmin, ymin):
        super().__init__(nno_x, nno_y, xmax, ymax, xmin, ymin)
        self.bc = constBC()
        self.rhs = polynomialRHS()

    def essBC(self, x, y, boind):
        retval = self.bc.essBC(x,y,boind)
        return retval


# Derived class to implement the right hand side  
# and boundary conditions for FEM solution
class FEMPoissonSub(FEM.PoissonFEM):

    def __init__(self, params_, grid_):
        super().__init__(params_, grid_)
        self.bc = constBC()
        self.rhs = polynomialRHS()
        
    # Only constant essbc's for now 
    def essBC(self, x, y, boind):
        retval = self.bc.essBC(x,y,boind)
        return retval

# Parses a float from a string 
def string2float(str):

    sign = 1
    lst = str.split('.')
    if len(lst) > 2:
        print("string2float: Error, ", str, " is not a valid floating point value")
        return 0.0

    if str[0] == '-':
        sign = -1
        str = str[1:]

    if len(lst) == 1:
        if str.isnumeric():
            return sign*float(str)

    if len(lst) == 2:
        if str.isnumeric() and lst[1].isnumeric():
            return sign*float(str)

    print(" S2F 2: Error, ", str, " is not a valid number ")
    return 0.0

# Functions to read from input form
def parse_rhs(arg, coeff):
    a = []
    lst = coeff.split(',')

    if arg == 'Constant':
        if len(lst) != 1:
            return []
        a.append(string2float(lst[0]))
        return a
    elif arg == 'Linear':
        if len(lst) != 3:
            print("Coefficient string is incorrect for linear rhs ", coeff)
            return []

        for i in range(3):
            a.append(string2float(lst[i].strip()))
        return a
    elif arg == 'Quadratic':
        if len(lst) != 6:
            print("Coefficient string is incorrect for quadratic rhs ", coeff)
            return []

        for i in range(6):
            a.append(string2float(lst[i].strip()))
        return a
    else:
        return []


def parse_bc(arg):
    a = []
    a.append(string2float(arg['bcright'].strip()))
    a.append(string2float(arg['bcupp'].strip()))
    a.append(string2float(arg['bcleft'].strip()))
    a.append(string2float(arg['bclow'].strip()))
  
    return a


def parse_domain(arg):
    temp = arg.strip()
    if temp[0] != '[':
        print("Invalid domain :", temp)
        return ()

    temp = temp[1:]
    pos = temp.find(',')
    x_min = string2float(temp[:pos].strip())
    temp = temp[(pos+1):]
    
    pos = temp.find(']')
    x_max = string2float(temp[:pos].strip())
       
    pos = temp.find('x')
    temp = temp[pos+1:].strip()  

    if temp[0] != '[':
        print("Invalid domain :", temp)
        return ()

    temp = temp[1:]
    pos = temp.find(',')
    y_min = string2float(temp[:pos])
    temp = temp[(pos+1):]
 
    pos = temp.find(']')
    y_max = string2float(temp[:pos].strip())

    return (x_min, x_max, y_min, y_max)

# Read input values from file
def read_from_file(values):

    filename = sg.PopupGetFile('Load input file', no_window=True)

    if not filename:
        print("No valid file provided")
        return
    try:
        file = open(filename, 'r')
        lines = file.readlines()
        validinputs = ('right hand side','number of nodes', 'domain','method of solution', 'rhs_coeff', 'bclow', 'bcupp', 'bcleft', 'bcright')
        for row in lines:
            if row.strip()[0] == '#':  # This is a comment
                print("Found a comment")
                continue
            else:
                ans = row.find(':')
                if ans < 0:
                    continue
                else:
                    entry = row[:ans].strip()
                    if entry.lower() in validinputs:
                        row = row[ans+1:].strip()
                        window.Element(entry).Update(row)

    except IOError:
        print("Error: File not found.")

if __name__ == '__main__':
    
    sg.theme('DarkAmber')  # Add a touch of color
    # All the stuff inside your window.
    layout = [[sg.Text('File', size=(12, 1)), sg.Button('Read input file'), sg.Button('Save input')],
              [sg.Text('Right hand side', size=(12, 1), auto_size_text=False, justification='left'),
               sg.InputOptionMenu(('Constant', 'Linear', 'Quadratic'), key='Right hand side', default_value='Constant', size=(7, 1)),
               sg.Text('Ax^2+By^2+Cx*y+Dx+Ey+F', size=(20, 1), auto_size_text=True, justification='left')],
              [sg.Text('Domain', size=(12, 1), auto_size_text=False, justification='left'), sg.InputText("[0,1] x [0,1]", key='Domain', size=(12, 1)),
               sg.InputText("1", key='rhs_coeff', size=(20, 1))],
              [sg.Text('Number of nodes', size=(12, 1), auto_size_text=False, justification='left'), sg.InputText("400", key='Number of nodes', size=(12, 1))],
              [sg.Text('Solution method', size=(12, 1), auto_size_text=False, justification='left'),
               sg.InputCombo(('FDM', 'FEM'),key='Method of solution', size=(5, 1),default_value='FDM')],
              [sg.Text('Boundary conditions', size=(20, 1), auto_size_text=False, justification='left')],
              [sg.Text('Lower boundary', size=(12, 1), auto_size_text=False, justification='left'),
               sg.InputText("0", key='bclow', size=(12, 1))],
              [sg.Text('Upper boundary', size=(12, 1), auto_size_text=False, justification='left'),
               sg.InputText("0", key='bcupp', size=(12, 1))],
              [sg.Text('Left boundary', size=(12, 1), auto_size_text=False, justification='left'),
               sg.InputText("0", key='bcleft', size=(12, 1))],
              [sg.Text('Right boundary', size=(12, 1), auto_size_text=False, justification='left'),
               sg.InputText("0", key='bcright', size=(12, 1))],
              [sg.Button('Run'), sg.Button('Cancel'), sg.Button("Stop")]]


    window = sg.Window('Poisson equation parameters', layout)
    # Event Loop to process "events" and get the "values" of the inputs
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Cancel':  # if user closes window or clicks cancel
            break

        if event == 'Read input file':
            read_from_file(values)

        if event == 'Save input':
            filename = sg.PopupGetFile('Write input to file', save_as=True, no_window=True)
            validinputs = ('right hand side', 'number of nodes', 'domain', 'method of solution', 'rhs_coeff', 'bclow', 'bcupp', 'bcleft', 'bcright')
            try:
                file = open(filename, 'w')
                for i in values:
                    if i.lower() in validinputs:
                        text = i + ':'
                        text += values[i] + '\n'
                        file.write(text)
                file.close()

            except IOError:
                print("Cannot write to file: ", filename)


        if event == 'Run':
            print('Running simulation')

            (x_min, x_max, y_min, y_max) = parse_domain(values['Domain'])

            nnos = values['Number of nodes']
            if not nnos.isnumeric():
                print("Number of nodes is not an integer")
                continue
            else:
                nno = int(nnos)

            sol_met = values['Method of solution']
            if (not sol_met):
                sol_met = 'FDM'

            # for now, assume a  box
            lx = (x_max - x_min)
            ly = (y_max - y_min)
            
            nno_x = int((nno*lx/ly)**0.5) 
            nno_y = int(nno/nno_x)

            if sol_met == 'FEM':
                print("Creating FEM solver")
                print("NOTE: There is  bug in the FEM solver. \n Works if rhs=0")
                parameters = dsc.DiscPrms(nnx=nno_x, nny=nno_y, dt=1000, t_max=2000)
                grid = dsc.Grid2d(parameters)
                # For now, assume only Dirichlet BC's 
                boind_list = [1,2,3,4]
                grid.setBoindWithEssBC(boind_list)                

                sim = FEMPoissonSub(parameters, grid)
                sim.rhs.attachRHS(parse_rhs(values['Right hand side'], values['rhs_coeff'].strip()))
                sim.bc.attachBC(parse_bc(values))
                if not sim.rhs.rhs:
                    print("Illegal right hand side ", sim.rhs.rhs)
                    continue
                
                sim.solve(-1)

            else: # Finite differences
                print("Creating FDM solver")
                solver = PoissonSub(nno_x, nno_y, x_max, y_max, x_min, y_min)
                solver.rhs.attachRHS(parse_rhs(values['Right hand side'], values['rhs_coeff'].strip()))
                solver.bc.attachBC(parse_bc(values))

                print(parse_bc(values))

                if not solver.rhs.rhs:
                    print("Illegal right hand side ", solver.rhs.rhs)
                    continue

                solver.solve()

        if event == 'Stop':
            print('Stopping simulation')
            break

    window.close()

