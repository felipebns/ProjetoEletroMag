from cmath import *
from numpy import linalg
import numpy as np
import matplotlib.pyplot as plt

class Simulation2():

    def __init__(self, C, L1, L2, Rint1, Rint2, Uf):
        self.colors = ['b','g','r','c','m', 'royalblue', 'slategrey', 'plum']
        self.C = C
        self.L1 = L1
        self.L2 = L2
        self.Rint1 = Rint1
        self.Rint2 = Rint2
        self.Uf = Uf*(2/np.pi) #eficaz quadrada
        self.lista_rc = [5, 10, 15, 20, 25, 30, 35, 40]
        self.lista_k = np.arange(0,0.5,0.01)

    def valida_saida(self) -> dict:

        rc_freq_vs_rend_s = {}
        rc_freq_vs_saida_s = {}
        freq_maior_pico = {}

        freqs = np.arange(0,100e3,10)
        k = 0.2 #Fixando

        for rc in self.lista_rc:
            rend_aux_s = []
            saida_aux_s = []

            for fr in freqs:
                w = 2*pi*fr

                v2_s, rendimento_s = self.calcula_rendimento_saida_serie(w, rc, k)
                rend_aux_s.append(rendimento_s)
                saida_aux_s.append(v2_s)
                
            rc_freq_vs_rend_s[rc] = rend_aux_s
            rc_freq_vs_saida_s[rc] = saida_aux_s

            i = 0
            maior = 0
            for saida in saida_aux_s:
                if saida > maior:
                    maior = saida
                    index = i
                i += 1

            maior_freq = freqs[index]
            freq_maior_pico[rc] = maior_freq

        fig, ax1 = plt.subplots(figsize=(20,14))

        ax1.set_xlabel('freq')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2 pico')
        ax2.tick_params(axis='y')

        for i, (rc, rendimento) in enumerate(rc_freq_vs_rend_s.items()):
            ax1.plot(freqs, rendimento, color=self.colors[i], label=f'Rendimento para Rc de: {rc} Ohms', linestyle=':')
            
            ax2.plot(freqs, list(rc_freq_vs_saida_s.values())[i], color=self.colors[i], label=f'Saída para Rc de: {rc} Ohms', linestyle='-')

        plt.title('Freq VS Pico & Rendimento - Série - K = 0.2')

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')
        plt.show()

        return freq_maior_pico
                

    def calculate(self) -> tuple:
        rc_k_vs_rend_s = {}
        rc_k_vs_saida_s = {}

        rc_k_vs_rend_p = {}
        rc_k_vs_saida_p = {}


        for rc in self.lista_rc:
            rend_aux_s = []
            saida_aux_s = []
            rend_aux_p = []
            saida_aux_p = []
            
            freq_maior_pico = self.valida_saida()
            freq1 = freq_maior_pico[rc]
            w1 = 2*np.pi*freq1

            for k in self.lista_k:
                v2_s, rendimento_s = self.calcula_rendimento_saida_serie(w1, rc, k)
                rend_aux_s.append(rendimento_s)
                saida_aux_s.append(v2_s)

                v2_p, rendimento_p = self.calcula_rendimento_saida_paralelo(w1, rc, k)
                rend_aux_p.append(rendimento_p)
                saida_aux_p.append(v2_p)
                
            rc_k_vs_rend_s[rc] = rend_aux_s
            rc_k_vs_saida_s[rc] = saida_aux_s

            rc_k_vs_rend_p[rc] = rend_aux_p
            rc_k_vs_saida_p[rc] = saida_aux_p

        return rc_k_vs_rend_s, rc_k_vs_saida_s, rc_k_vs_rend_p, rc_k_vs_saida_p

    def plot(self) -> None:
        rc_k_vs_rend_s, rc_k_vs_saida_s, rc_k_vs_rend_p, rc_k_vs_saida_p = self.calculate()

        # plt.plot(self.lista_k, rc_k_vs_saida_s)

        self.plot_serie(rc_k_vs_rend_s, rc_k_vs_saida_s)
        # self.plot_paralelo(rc_k_vs_rend_p, rc_k_vs_saida_p)
        
        plt.show()

    def calcula_rendimento_saida_serie(self, w, rc, k) -> tuple:
        Uf = self.Uf
        XR1 = self.Rint1
        XR2 = self.Rint2
        XC = 1/((w*self.C)*1j)
        XL1 = (w*self.L1)*1j
        XL2 = (w*self.L2)*1j
        XM = k*((XL1*XL2)**0.5)
        XRc = rc

        Z=np.array([[XR1+XL1+XC, -XM],[-XM, XL2+XR2+XRc+XC]])
        V=np.array([Uf,0])
        i=np.dot(linalg.inv(Z),V)

        i1, i2 = i[0], i[1]

        v2_s = XRc*i2
        Psaida = np.real((v2_s*np.conj(i2))/2)
        Pentrada = np.real((Uf*np.conj(i1))/2)
        rendimento = Psaida/Pentrada

        return abs(v2_s), rendimento*100

    def calcula_rendimento_saida_paralelo(self, w, rc, k) -> tuple:
        Uf = self.Uf
        XR1 = self.Rint1
        XR2 = self.Rint2
        XC = 1/((w*self.C)*1j)
        XL1 = (w*self.L1)*1j
        XL2 = (w*self.L2)*1j
        XM = k*((XL1*XL2)**0.5)
        XRc = rc
        Zcarga = ((XRc*XC)/(XRc+XC))

        Z=np.array([[XR1+XL1+XC, -XM],[-XM, XL2+XR2+Zcarga]])
        V=np.array([Uf,0])
        i=np.dot(linalg.inv(Z),V)

        i1, i2 = i[0], i[1]

        v2_p = Zcarga*i2
        ic = (XC*i2)/(XRc+XC)
        Psaida = np.real((v2_p*np.conj(ic))/2)
        Pentrada = np.real((Uf*np.conj(i1))/2)
        rendimento = Psaida/Pentrada

        return abs(v2_p), rendimento*100

    def plot_serie(self, rc_k_vs_rend_s, rc_k_vs_saida_s) -> None:
        fig, ax1 = plt.subplots(figsize=(20,14))

        ax1.set_xlabel('K')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2 pico')
        ax2.tick_params(axis='y')

        for i, (rc, rendimento) in enumerate(rc_k_vs_rend_s.items()):
            ax1.plot(self.lista_k, rendimento, color=self.colors[i], label=f'Rendimento para Rc de: {rc} Ohms', linestyle=':')
            
            ax2.plot(self.lista_k, list(rc_k_vs_saida_s.values())[i], color=self.colors[i], label=f'Saída para Rc de: {rc} Ohms', linestyle='-')

        plt.title('Fator de Acoplamento (K) VS Rendimento & Saída ---- Em série')

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

    def plot_paralelo(self, rc_k_vs_rend_p, rc_k_vs_saida_p) -> None:
        fig, ax1 = plt.subplots(figsize=(20,14))

        ax1.set_xlabel('K')
        ax1.set_ylabel('Rendimento %')
        ax1.tick_params(axis='y')

        ax2 = ax1.twinx()
        ax2.set_ylabel('Saída V2 pico')
        ax2.tick_params(axis='y')

        for i, (rc, rendimento) in enumerate(rc_k_vs_rend_p.items()):
            ax1.plot(self.lista_k, rendimento, color=self.colors[i], label=f'Rendimento para Rc de: {rc} Ohms', linestyle=':')
            
            ax2.plot(self.lista_k, list(rc_k_vs_saida_p.values())[i], color=self.colors[i], label=f'Saída para Rc de: {rc} Ohms', linestyle='-')

        plt.title('Fator de Acoplamento (K) VS Rendimento & Saída ---- Em paralelo')

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

    def run(self) -> None:
        self.valida_saida()
        # self.plot()