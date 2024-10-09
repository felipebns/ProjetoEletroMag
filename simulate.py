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
    
    def calculate(self) -> tuple:
        cap_freq_vs_rend_s = {}
        cap_freq_vs_saida_s = {}

        cap_freq_vs_rend_p = {}
        cap_freq_vs_saida_p = {}

        for C in self.lista_capac:
            rend_aux_s = []
            saida_aux_s = []
            rend_aux_p = []
            saida_aux_p = []
            
            for fr in self.lista_freq:
                v1 = self.Uf
                w = 2*pi*fr
                L = ((1/(2*pi*fr))**2)*(1/C)
                R = self.Rdc + (self.Rdc/(100e3))*fr

                i1_s, i2_s = self.calculate_current_serie(C, L, R, w)
                i1_p, i2_p = self.calculate_current_parallel(C, L, R, w)

                v2_s, rendimento_s = self.calculate_serie(i1_s, i2_s, v1)
                rend_aux_s.append(rendimento_s)
                saida_aux_s.append(v2_s)

                v2_p, rendimento_p = self.calculate_parallel(i1_p, i2_p, v1)
                rend_aux_p.append(rendimento_p)
                saida_aux_p.append(v2_p)
                
            cap_freq_vs_rend_s[C] = rend_aux_s
            cap_freq_vs_saida_s[C] = saida_aux_s

            cap_freq_vs_rend_p[C] = rend_aux_p
            cap_freq_vs_saida_p[C] = saida_aux_p

        return cap_freq_vs_rend_s, cap_freq_vs_saida_s, cap_freq_vs_rend_p, cap_freq_vs_saida_p
    
    def calculate_serie(self, i1_s, i2_s, v1):
        v2_s = self.Rc*i2_s
        Psaida = (v2_s*np.conj(i2_s))/2
        Pentrada = (v1*np.conj(i1_s))/2
        rendimento = Psaida/Pentrada

        return abs(v2_s), rendimento*100
    
    def calculate_parallel(self, i1_p, i2_p, v1):
        v2_p = self.Rc*i2_p
        Psaida = (v2_p*np.conj(i2_p))/2
        Pentrada = (v1*np.conj(i1_p))/2
        rendimento = Psaida/Pentrada

        return abs(v2_p), rendimento*100

    def plot(self) -> None:
        cap_freq_vs_rend_s, cap_freq_vs_saida_s, cap_freq_vs_rend_p, cap_freq_vs_saida_p = self.calculate()

        self.plot_serie(cap_freq_vs_rend_s, cap_freq_vs_saida_s)
        self.plot_parallel(cap_freq_vs_rend_p, cap_freq_vs_saida_p)
        
        plt.show()

    def plot_serie(self, cap_freq_vs_rend_s, cap_freq_vs_saida_s) -> None:
        # Create the figure and the first axis outside the loop
        fig, ax1 = plt.subplots(figsize=(20,14))

        ax1.set_xlabel('Frequências')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2_s pico')
        ax2.tick_params(axis='y')

        for i, (c, rendimento) in enumerate(cap_freq_vs_rend_s.items()):
            ax1.plot(self.lista_freq, rendimento, color=self.colors[i], label=f'Rendimento para capacitância de: {c}F', linestyle=':')
            
            ax2.plot(self.lista_freq, list(cap_freq_vs_saida_s.values())[i], color=self.colors[i], label=f'Saída para capacitância de: {c}F', linestyle='-')

        plt.title('Frequência VS Rendimento & Saída ---- Em série')

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

    def plot_parallel(self, cap_freq_vs_rend_p, cap_freq_vs_saida_p) -> None:
        fig, ax1 = plt.subplots(figsize=(20,14))

        ax1.set_xlabel('Frequências')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2_s pico')
        ax2.tick_params(axis='y')

        for i, (c, rendimento) in enumerate(cap_freq_vs_rend_p.items()):
            ax1.plot(self.lista_freq, rendimento, color=self.colors[i], label=f'Rendimento para capacitância de: {c}F', linestyle=':')
            
            ax2.plot(self.lista_freq, list(cap_freq_vs_saida_p.values())[i], color=self.colors[i], label=f'Saída para capacitância de: {c}F', linestyle='-')

        plt.title('Frequência VS Rendimento & Saída ---- Em paralelo')

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

    def run(self) -> None:
        self.plot()