import numpy as np
from scipy import signal as sg
FS=500

def ecg(dur=180.0, hr=62.0, amp_r=1.30, amp_t=0.30, sigma_t=0.045,
        modulacion=0.0, desvanecer=False, semilla=7, fs=FS):
    """ECG configurable. modulacion: variacion de amplitud del QRS (0-1).
    desvanecer: la amplitud cae linealmente a la mitad en la segunda mitad."""
    rng=np.random.default_rng(semilla)
    n=int(fs*dur); t=np.arange(n)/fs; x=np.zeros(n); tr=0.5; p=[]
    while tr<dur-0.6:
        p.append(tr)
        ka = 1 + modulacion*np.sin(2*np.pi*0.20*tr)
        if desvanecer: ka *= (1.0 - 0.55*max(0.0,(tr-dur/2))/(dur/2))
        for c,a,s in [(-0.200,0.15,0.025),(-0.035,-0.10,0.008),(0.0,1.10,0.010),
                      (0.035,-0.25,0.010),(0.280,amp_t,sigma_t)]:
            x += a*ka*np.exp(-0.5*((t-(tr+c))/s)**2)
        tr += max(0.35, rng.normal(60/hr,0.045))
    return t, x*(amp_r/1.10), np.array(p)

def ruido(t, fs=FS, emg_amp=0.03, semilla=11):
    rng=np.random.default_rng(semilla)
    dv=0.35*np.sin(2*np.pi*0.28*t)+0.18*np.sin(2*np.pi*0.09*t+1.1)
    rd=0.12*np.sin(2*np.pi*60*t)+0.034*np.sin(2*np.pi*120*t+0.7)
    sos=sg.butter(4,[20,150],'bandpass',fs=fs,output='sos')
    return dv+rd+emg_amp*sg.sosfilt(sos,rng.standard_normal(len(t)))

ESCENARIOS = {
 'normal':            dict(),
 'T picuda':          dict(amp_t=0.70, sigma_t=0.028),
 'amplitud variable': dict(modulacion=0.55),
 'señal que decae':   dict(desvanecer=True),
 'taquicardia':       dict(hr=145),
}
