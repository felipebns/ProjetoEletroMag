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
                w = 2*pi*fr
                L = ((1/(2*pi*fr))**2)*(1/C)
                R = self.Rdc + (self.Rdc/(100e3))*fr

                v2_s, rendimento_s = self.calcula_rendimento_saida_serie(C, L, R, w)
                rend_aux_s.append(rendimento_s)
                saida_aux_s.append(v2_s)

                v2_p, rendimento_p = self.calcula_rendimento_saida_paralelo(C, L, R, w)
                rend_aux_p.append(rendimento_p)
                saida_aux_p.append(v2_p)
                
            cap_freq_vs_rend_s[C] = rend_aux_s
            cap_freq_vs_saida_s[C] = saida_aux_s

            cap_freq_vs_rend_p[C] = rend_aux_p
            cap_freq_vs_saida_p[C] = saida_aux_p

        return cap_freq_vs_rend_s, cap_freq_vs_saida_s, cap_freq_vs_rend_p, cap_freq_vs_saida_p

    def plot(self) -> None:
        cap_freq_vs_rend_s, cap_freq_vs_saida_s, cap_freq_vs_rend_p, cap_freq_vs_saida_p = self.calculate()

        self.plot_serie(cap_freq_vs_rend_s, cap_freq_vs_saida_s)
        self.plot_paralelo(cap_freq_vs_rend_p, cap_freq_vs_saida_p)
        
        plt.show()

    def calcula_rendimento_saida_serie(self, C, L, R, w) -> tuple:
        Uf = self.Uf
        XRc = self.Rc
        XR = R
        XC = 1/((w*C)*1j)
        XL = (w*L)*1j
        XM = self.k*XL

        Z=np.array([[XR+XL+XC, -XM],[-XM, XL+XR+XRc+XC]])
        V=np.array([Uf,0])
        i=np.dot(linalg.inv(Z),V)

        i1, i2 = i[0], i[1]

        v2_s = self.Rc*i2
        Psaida = (v2_s*np.conj(i2))/2
        Pentrada = (Uf*np.conj(i1))/2
        rendimento = Psaida/Pentrada

        return abs(v2_s), rendimento*100

    def calcula_rendimento_saida_paralelo(self, C, L, R, w) -> tuple:
        Uf = self.Uf
        XRc = self.Rc
        XR = R
        XC = 1/((w*C)*1j)
        XL = (w*L)*1j
        XM = self.k*XL
        Zcarga = ((XRc*XC)/(XRc+XC))

        Z=np.array([[XR+XL+XC, -XM],[-XM, XL+XR+Zcarga]])
        V=np.array([Uf,0])
        i=np.dot(linalg.inv(Z),V)

        i1, i2 = i[0], i[1]

        v2_p = Zcarga*i2
        Psaida = (v2_p*np.conj(i2))/2
        Pentrada = (Uf*np.conj(i1))/2
        rendimento = Psaida/Pentrada

        return abs(v2_p), rendimento*100

    def plot_serie(self, cap_freq_vs_rend_s, cap_freq_vs_saida_s) -> None:
        fig, ax1 = plt.subplots(figsize=(20,14))

        ax1.set_xlabel('Frequências')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2 pico')
        ax2.tick_params(axis='y')

        for i, (c, rendimento) in enumerate(cap_freq_vs_rend_s.items()):
            ax1.plot(self.lista_freq, rendimento, color=self.colors[i], label=f'Rendimento para capacitância de: {c}F', linestyle=':')
            
            ax2.plot(self.lista_freq, list(cap_freq_vs_saida_s.values())[i], color=self.colors[i], label=f'Saída para capacitância de: {c}F', linestyle='-')

        plt.title('Frequência VS Rendimento & Saída ---- Em série')

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

    def plot_paralelo(self, cap_freq_vs_rend_p, cap_freq_vs_saida_p) -> None:
        fig, ax1 = plt.subplots(figsize=(20,14))

        ax1.set_xlabel('Frequências')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2 pico')
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