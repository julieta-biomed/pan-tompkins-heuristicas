import numpy as np
from scipy import signal as sg

def paso_banda(x, fs, lo=5.0, hi=15.0, orden=3):
    """Etapa 1: aislar la banda donde domina el QRS."""
    sos = sg.butter(orden,[lo,hi],'bandpass',fs=fs,output='sos')
    return sg.sosfiltfilt(sos,x)

def derivada(x, fs):
    """Etapa 2: filtro diferenciador de 5 puntos del Pan-Tompkins original.

    y[n] = (1/8T)(-x[n-2] - 2x[n-1] + 2x[n+1] + x[n+2])
    Resalta pendientes abruptas; aplana ondas P y T, que son lentas.
    """
    h = np.array([-1,-2,0,2,1], dtype=float)*(fs/8.0)
    return np.convolve(x, h[::-1], mode='same')

def cuadrado(x):
    """Etapa 3: no linealidad. Todo positivo y los picos grandes se amplifican."""
    return x**2

def integrar(x, fs, ventana=0.150):
    """Etapa 4: media movil. Convierte cada QRS en una sola joroba."""
    n = max(1,int(round(ventana*fs)))
    return np.convolve(x, np.ones(n)/n, mode='same')

def cadena(x, fs, ventana=0.150, banda=(5.0,15.0)):
    b = paso_banda(x,fs,*banda)
    d = derivada(b,fs)
    c = cuadrado(d)
    i = integrar(c,fs,ventana)
    return b,d,c,i
