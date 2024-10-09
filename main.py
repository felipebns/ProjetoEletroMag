from simulate import Simulation

if __name__ == '__main__':
    Rc = 5 #simplificação -> chute
    Rdc = 0.2 #simplificação -> chute
    k = 0.2 #fator de acoplamento -> chute
    Uf = 10 #eficaz -> chute
    simulation = Simulation(Rc=Rc, Rdc=Rdc, k=k, Uf=Uf)
    simulation.run()