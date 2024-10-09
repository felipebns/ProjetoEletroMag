from cmath import *
from numpy import linalg
import numpy as np
import matplotlib.pyplot as plt

class Simulation():

    def __init__(self, Rc, Rdc, k, Uf):
        self.lista_freq = np.arange(0,100e3,10)
        self.lista_capac = [150e-9, 0.1e-6, 0.47e-6, 1e-6, 4.7e-9]
        self.colors = ['b','g','r','c','m']
        self.Rc = Rc
        self.Rdc = Rdc
        self.k = k
        self.Uf = Uf

    def calculate_current_serie(self, C, L, R, w) -> tuple:
        XRc = self.Rc
        XR = R
        XC = 1/((w*C)*1j)
        XL = (w*L)*1j
        XM = self.k*XL

        Z=np.array([[XR+XL+XC, -XM],[-XM, XL+XR+XRc+XC]])
        V=np.array([self.Uf,0])
        i=np.dot(linalg.inv(Z),V)
        return i[0], i[1]

    def calculate_current_parallel(self, C, L, R, w) -> tuple:
        XRc = self.Rc
        XR = R
        XC = 1/((w*C)*1j)
        XL = (w*L)*1j
        XM = self.k*XL
        Z=np.array([[XR+XL+XC, -XM],[-XM, XL+XR+((XRc*XC)/(XRc+XC))]])
        V=np.array([self.Uf,0])
        i=np.dot(linalg.inv(Z),V)
        return i[0], i[1]
    
    def calculate(self):
        c_para_rendimento = {}
        c_para_saida = {}

        c_para_rendimento_p = {}
        c_para_saida_p = {}

        for C in self.lista_capac:
            rend_aux = []
            saida_aux = []
            rend_aux_p = []
            saida_aux_p = []
            for fr in self.lista_freq:
                w = 2*pi*fr
                L = ((1/(2*pi*fr))**2)*(1/C)
                R = self.Rdc + (self.Rdc/(100e3))*fr
                i1, i2 = self.calculate_current_serie(C, L, R, w)
                i1_p, i2_p = self.calculate_current_parallel(C, L, R, w)

                v2 = self.Rc*i2
                v1 = self.Uf
                Psaida = (v2*np.conj(i2))/2
                Pentrada = (v1*np.conj(i1))/2
                rendimento = Psaida/Pentrada
                rend_aux.append(rendimento*100)
                saida_aux.append(abs(v2)) #PICO

                v2_p = self.Rc*i2_p
                Psaida_p = (v2_p*np.conj(i2_p))/2
                Pentrada_p = (v1*np.conj(i1_p))/2
                rendimento_p = Psaida_p/Pentrada_p
                rend_aux_p.append(rendimento_p*100)
                saida_aux_p.append(abs(v2_p)) #PICO
                
            c_para_rendimento[C] = rend_aux
            c_para_saida[C] = saida_aux

            c_para_rendimento_p[C] = rend_aux_p
            c_para_saida_p[C] = saida_aux_p

        return c_para_rendimento, c_para_saida, c_para_rendimento_p, c_para_saida_p
    
    def plot(self):
        c_para_rendimento, c_para_saida, c_para_rendimento_p, c_para_saida_p = self.calculate()
        self.plot_serie(c_para_rendimento, c_para_saida)
        self.plot_parallel(c_para_rendimento_p, c_para_saida_p)

    def plot_serie(self, c_para_rendimento, c_para_saida) -> None:
        # Create the figure and the first axis outside the loop
        fig, ax1 = plt.subplots(figsize=(20,14))

        # Set up the first y-axis for 'Rendimento'
        ax1.set_xlabel('Frequências')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        # Create a second y-axis sharing the same x-axis
        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2 pico')
        ax2.tick_params(axis='y')

        # Loop over the items and plot each one on the respective axis
        for i, (c, rendimento) in enumerate(c_para_rendimento.items()):
            # Plot on the first y-axis
            ax1.plot(self.lista_freq, rendimento, color=self.colors[i], label=f'Rendimento para capacitância de: {c}F', linestyle=':')
            
            # Plot on the second y-axis
            ax2.plot(self.lista_freq, list(c_para_saida.values())[i], color=self.colors[i], label=f'Saída para capacitância de: {c}F', linestyle='-')

        # Add a title
        plt.title('Frequência VS Rendimento & Saída ---- Em série')

        # Combine the legends from both axes
        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

        # Display the plot
        plt.show()

    def plot_parallel(self, c_para_rendimento_p, c_para_saida_p) -> None:
        # Create the figure and the first axis outside the loop
        fig, ax1 = plt.subplots(figsize=(20,14))

        # Set up the first y-axis for 'Rendimento'
        ax1.set_xlabel('Frequências')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        # Create a second y-axis sharing the same x-axis
        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2 pico')
        ax2.tick_params(axis='y')

        # Loop over the items and plot each one on the respective axis
        for i, (c, rendimento) in enumerate(c_para_rendimento_p.items()):
            # Plot on the first y-axis
            ax1.plot(self.lista_freq, rendimento, color=self.colors[i], label=f'Rendimento para capacitância de: {c}F', linestyle=':')
            
            # Plot on the second y-axis
            ax2.plot(self.lista_freq, list(c_para_saida_p.values())[i], color=self.colors[i], label=f'Saída para capacitância de: {c}F', linestyle='-')

        # Add a title
        plt.title('Frequência VS Rendimento & Saída ---- Em paralelo')

        # Combine the legends from both axes
        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

        # Display the plot
        plt.show()



    def run(self) -> None:
        self.plot()