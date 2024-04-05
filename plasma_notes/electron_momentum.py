import scipy.constants as sci
import numpy as np
import plotly.graph_objects as go

def electron_momentum(pho_wavelen, scatt_angle):
    en_diff = pho_energy_diff(pho_wavelen, scatt_angle)
    eon_momentum = np.sqrt((en_diff ** 2 + 2 * en_diff * sci.electron_mass * sci.c ** 2) / sci.c ** 2)
    return eon_momentum

def pho_energy_diff(pho_wavelen, scatt_angle):
    pho_en = sci.h * sci.c / pho_wavelen
    common_term_1 = 1 - np.cos(scatt_angle)
    common_term_2 = sci.electron_mass * sci.c ** 2
    en_diff = pho_en ** 2 / (common_term_2) * common_term_1 / (1 + pho_en / common_term_2 * common_term_1)
    return en_diff

def wavelen_from_ev(en_ev):
    wavelen = sci.h * sci.c / sci.e * (1 / en_ev)
    return wavelen

pho_en = np.arange(0, 1000, 0.1)
p_eon_max = np.array([electron_momentum(wavelen_from_ev(i), np.pi / 2) for i in pho_en])
p_eon_mid = np.array([electron_momentum(wavelen_from_ev(i), np.pi / 4) for i in pho_en])
fig = go.Figure()
fig.add_trace(go.Scatter(x=pho_en, y=p_eon_max, mode='lines'))
fig.add_trace(go.Scatter(x=pho_en, y=p_eon_mid, mode='lines'))
fig.show()
