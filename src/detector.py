import numpy as np
from scipy import signal as sg
from pantom import cadena

def detectar(x, fs, umbral_adaptativo=True, refractario=True, regla_360=True,
             search_back=True, ventana=0.150, k_inicial=0.25):
    """Pan-Tompkins con sus heuristicas conmutables, para ablacion.

    Devuelve las posiciones (en segundos) de los picos R detectados.
    """
    b,d,c,integrada = cadena(x, fs, ventana=ventana)
    n = len(integrada)

    # candidatos: maximos locales de la senal integrada.
    # OJO: el parametro distance de find_peaks YA es un periodo refractario
    # implicito. Para que la ablacion sea honesta, se desactiva junto con el.
    dist = int(0.10*fs) if refractario else 1
    cand,_ = sg.find_peaks(integrada, distance=dist)
    if len(cand)==0: return np.array([]), {}

    # inicializacion sobre los primeros 2 s (fase de aprendizaje del original)
    ini = integrada[:int(2*fs)]
    SPKI = 0.25*ini.max() if len(ini) else integrada.max()*0.25
    NPKI = 0.5*ini.mean() if len(ini) else 0.0
    UMBRAL = NPKI + 0.25*(SPKI-NPKI)
    if not umbral_adaptativo:
        UMBRAL = k_inicial*np.percentile(integrada,99)

    picos=[]; rr=[]; rr_medio=None; ultimo=-np.inf
    # La "pendiente maxima" del original se mide sobre la senal DERIVADA,
    # no sobre la integrada: en el pico de la integrada la pendiente es cero.
    absd = np.abs(d)
    w_p = int(0.05*fs)
    def pend_max(q):
        a,b_ = max(0,q-w_p), min(n,q+w_p)
        return absd[a:b_].max()
    aceptados_idx=[]

    i=0
    while i < len(cand):
        p = cand[i]; tp = p/fs
        if integrada[p] > UMBRAL:
            # --- periodo refractario ---
            if refractario and (tp-ultimo) < 0.200:
                i+=1; continue
            # --- regla de los 360 ms: descartar onda T ---
            if regla_360 and picos and 0.200 <= (tp-ultimo) < 0.360:
                if pend_max(p) < 0.5*pend_max(int(ultimo*fs)):
                    if umbral_adaptativo:            # se trata como ruido
                        NPKI = 0.125*integrada[p] + 0.875*NPKI
                        UMBRAL = NPKI + 0.25*(SPKI-NPKI)
                    i+=1; continue
            picos.append(tp); aceptados_idx.append(p)
            if ultimo>-np.inf: rr.append(tp-ultimo)
            ultimo=tp
            if umbral_adaptativo:
                SPKI = 0.125*integrada[p] + 0.875*SPKI
                UMBRAL = NPKI + 0.25*(SPKI-NPKI)
        else:
            if umbral_adaptativo:
                NPKI = 0.125*integrada[p] + 0.875*NPKI
                UMBRAL = NPKI + 0.25*(SPKI-NPKI)

        # --- search-back ---
        if search_back and len(rr)>=3:
            rr_medio=np.mean(rr[-8:])
            if i+1<len(cand):
                sig_t = cand[i+1]/fs
                if (sig_t-ultimo) > 1.66*rr_medio:
                    # buscar el mayor candidato en el hueco con umbral a la mitad
                    hueco=[q for q in cand if ultimo+0.200 < q/fs < sig_t]
                    if hueco:
                        mejor=max(hueco,key=lambda q: integrada[q])
                        if integrada[mejor] > 0.5*UMBRAL:
                            tm=mejor/fs
                            picos.append(tm); aceptados_idx.append(mejor)
                            rr.append(tm-ultimo); ultimo=tm
                            if umbral_adaptativo:
                                SPKI=0.25*integrada[mejor]+0.75*SPKI
                                UMBRAL=NPKI+0.25*(SPKI-NPKI)
        i+=1

    picos=np.array(sorted(picos))
    return picos, dict(integrada=integrada, umbral_final=UMBRAL)

def evaluar(det, verd, tol=0.15):
    """Sensibilidad y valor predictivo positivo."""
    if len(det)==0: return 0.0,0.0,0,len(verd),0
    usados=set(); VP=0
    for tv in verd:
        j=int(np.argmin(np.abs(det-tv)))
        if abs(det[j]-tv)<tol and j not in usados:
            usados.add(j); VP+=1
    FN=len(verd)-VP; FP=len(det)-VP
    Se = 100*VP/max(VP+FN,1); VPP = 100*VP/max(VP+FP,1)
    return Se,VPP,VP,FN,FP
