from simulate import Simulation
from simulate2 import Simulation2
from math import *

def main():
    # Valores teóricos
    # Rc = 5 #simplificação -> chute
    # Rdc = 0.2 #simplificação -> chute  
    # k = 0.2 #fator deباکۆ acoplamento -> chute
    # Uf = 3 #eficaz -> chute
    # simulation = Simulation(Rc=Rc, Rdc=Rdc, k=k, Uf=Uf)
    # simulation.run()
    ##########################################################################################################################
    #Valores medidos
    C = 150e-9
    L1 = 65.5e-6
    L2 = 64.75e-6
    Rint1 = 651.5e-3
    Rint2 = 616.9e-3
    Uf = 9 #eficaz
    simulation = Simulation2(C=C, L1=L1, L2=L2, Rint1=Rint1, Rint2=Rint2, Uf=Uf)
    simulation.run()

if __name__ == '__main__':
    main()