"""
SPICE 2 - A program that takes in a netlist file and returns the voltages at each node and current in each branch.
"""

import sys 
import numpy as np

infinity = 1e50
epsilon = 1 / infinity
start = -1; end = -1
elements = {'R' : [], 'C' : [],'L' : [], 'V' : [], 'I' : []}
nodes_mapping = {"GND" : 0}
omega = epsilon
ac_flag = False

class Element:
    
    def __init__(self,line):
        tokens = line.split()
        self.name = tokens[0]
        self.value = complex(tokens[-1])
        self.type = tokens[0][0]
        
        # checking for erroneous elements
        if self.type not in 'RLCIV':
            print("Improper element type")
            sys.exit(0)
        
        self.n1 = tokens[1]
        self.n2 = tokens[2]
        
        # checking for the AC part of the circuit
        if len(tokens) == 6 and tokens[3] == 'ac' and self.type in 'IV':
            self.value = complex(tokens[-2]) / 2
            self.phase = float(tokens[-1])
            self.value = self.value*(np.cos(self.phase) + 1j*np.sin(self.phase))
        elif len(tokens) == 4 and self.type in 'RLC':
            if self.type == 'L':
                self.value *=omega*1j
            elif self.type == 'C':
                self.value = -1j/(self.value*omega)


def get_data():

    if len(sys.argv) != 2:                            
        sys.exit('Usage: assignment2.py <filename>.netlist')       
            
        # Opening the file and saving its content

    if not sys.argv[1].endswith('.netlist'):
        print("Please enter a .netlist only")
        sys.exit(0)

    try:
        name = sys.argv[1]
        f = open(name)
    except FileNotFoundError:
        print("File not found")
        sys.exit(0)

    # checking for .circuit, .end and .ac
    lines = []
    for line in f.readlines():
        lines.append(line.split('#')[0].split('\n')[0])
    
    try:
        start = lines.index(".circuit")
    except:
        print("File must have a .circuit line")
        sys.exit(0)
    
    try:
        end = lines.index(".end")
    except:
        print("File must have a .end line")
        sys.exit(0)

    circuit = lines[start+1:end]
    for i in range(end+1, len(lines)):
        if lines[i].startswith('.ac'):
            global ac_flag, omega
            ac_flag = True
            line = lines[i].split()
            omega = float(line[-1])
            print(omega)

    
    for line in circuit:
        e = Element(line)
        elements[e.type].append(e)
        if e.n1 not in nodes_mapping.keys():
            nodes_mapping[e.n1] = len(nodes_mapping)
        if e.n2 not in nodes_mapping.keys():
            nodes_mapping[e.n2] = len(nodes_mapping)


def make_matrix():
    n = len(nodes_mapping) - 1
    m = len(elements['V'])
    
    G = np.zeros((n,n),dtype='complex')
    B = np.zeros((n,m),dtype='complex')
    D = np.zeros((m,m),dtype='complex')

    for x in 'RLC':
        for i in range(len(elements[x])):
            e = elements[x][i]
            admittance = 1/e.value
            n1 = nodes_mapping[e.n1]
            n2 = nodes_mapping[e.n2]
            if n1 == 0 or n2 == 0:
                if n1 < n2:
                    n1,n2 = n2,n1
                n1 -= 1
                G[n1,n1] += admittance
            else:
                n1 -= 1
                n2 -= 1
                G[n1,n1] += admittance
                G[n2,n2] += admittance
                G[n1,n2] -= admittance
                G[n2,n1] -= admittance

    for i in range(m):
        e = elements['V'][i]
        n1 = nodes_mapping[e.n1]
        n2 = nodes_mapping[e.n2]
        if n1 != 0:
            B[n1-1,i] = 1
        if n2 != 0:
            B[n2-1,i] = -1
    
    C= np.transpose(B)

    GB = np.concatenate((G,B),axis=1)
    CD = np.concatenate((C,D),axis=1)
    A = np.concatenate((GB,CD),axis=0)

    b = np.zeros((m+n,1),dtype='complex')

    for i in range(len(elements['I'])):
        e = elements['I'][i]
        n1 = nodes_mapping[e.n1]
        n2 = nodes_mapping[e.n2]
        if n1 != 0:
            b[n1-1] += e.value
        if n2 != 0:
            b[n2-1] =- e.value

    for i in range(m):
        e = elements['V'][i]
        b[i+n] = e.value

    solve(A,b)

def solve(A,b):
    x=np.linalg.solve(A,b)
    display(x)

def get_key(d,val):
    for key, value in d.items():
        if d[key] == val:
            return str(key)
    return None

def display(x):
    np.set_printoptions(precision = 4,suppress = True)

    for i in range(len(nodes_mapping)-1):
        print("V_" + get_key(nodes_mapping,i+1) + " = ", np.format_float_scientific(np.real(x[i]),precision=7), "+", np.format_float_scientific(np.imag(x[i]),precision=7), "\bj V")
    for i in range(len(elements['V'])):
        print("I_" + elements['V'][i].n1 + "_"+elements['V'][i].n2 + " = ", np.format_float_scientific(np.real(x[i+len(nodes_mapping)-1]), precision=4), "+", np.format_float_scientific(np.imag(x[i+len(nodes_mapping)-1]), precision=4), "\bj A")
    if ac_flag:
        print("The Angluar Frequency of these values are", omega,"rad/s")

    
if __name__ == "__main__":
    get_data()
    make_matrix()