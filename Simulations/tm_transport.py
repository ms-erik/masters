from constants import *
import numpy as np
from dataclasses import dataclass

@dataclass
class ScatterMatrix:
    rL: complex
    rR: complex
    tLR: complex
    tRL: complex


def identity_scatter():
    return ScatterMatrix(rL=0, rR=0, tLR=1, tRL=1)


def interface_scatter(k1, k2, m1, m2):
    q1 = k1 / m1
    q2 = k2 / m2
    den = q1 + q2

    if abs(den) < 1e-30:
        den = 1e-30 + 1e-30j

    rL = (q1 -q2) / den
    rR = (q2 - q1) / den
    tLR = 2 * q1 / den
    tRL = 2 * q2 / den

    return ScatterMatrix(rL=rL, rR=rR, tLR=tLR, tRL=tRL)


def propagation_scatter(k, d_nm):
    d_m = d_nm * nm_to_m
    phi = k * d_m
    p = np.exp(1j * phi)
    return ScatterMatrix(rL=0, rR=0, tLR=p, tRL=p)


def star_product(S1: ScatterMatrix, S2: ScatterMatrix):
    den = 1 - S1.rR * S2.rL
    rL = S1.rL + S1.tRL * S2.rL * S1.tLR / (den)
    rR = S2.rR + S2.tLR * S1.rR * S2.tRL / (den)
    tLR = S1.tLR * S2.tLR / (den)
    tRL = S1.tRL * S2.tRL / (den)
    return ScatterMatrix(rL=rL, rR=rR, tLR=tLR, tRL=tRL)
