"""
experimento_schwarzschild.py
============================

Primer experimento oficial del Motor de Geometría Diferencial (MGD).

Objetivo
--------
Comparar la métrica de Schwarzschild clásica contra la
Métrica Disforme Regularizada propuesta por la teoría de
Atracción Energética.

Se calculan:

    • f(r)
    • Horizonte efectivo
    • Tiempo propio
    • Distancia propia
    • Comparación de ambas métricas

Las geodésicas completas serán incorporadas en la siguiente
iteración utilizando el integrador del MGD.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------
# Constantes físicas
# ----------------------------------------------------------

G = 6.67430e-11
c = 299792458

hbar = 1.054571817e-34

# Planck

l_planck = sp.sqrt(hbar * G / c**3)

m_planck = sp.sqrt(hbar * c / G)


# ----------------------------------------------------------
# Parámetros del cuerpo central
# ----------------------------------------------------------

M = sp.Symbol(
    "M",
    positive=True,
)

r = sp.Symbol(
    "r",
    positive=True,
)

# Schwarzschild

r_s = 2 * G * M / c**2

# Radio mínimo

r_min = l_planck * (M / m_planck) ** sp.Rational(1, 3)


# ----------------------------------------------------------
# Función constitutiva
# ----------------------------------------------------------

sigma = 1 / (
    1
    +
    (r_min / r) ** 4
)


# ----------------------------------------------------------
# Métrica clásica
# ----------------------------------------------------------

f_schwarzschild = 1 - r_s / r


# ----------------------------------------------------------
# Métrica propuesta
# ----------------------------------------------------------

f_disforme = 1 - (r_s / r) * sigma


# ----------------------------------------------------------
# Horizonte efectivo
# ----------------------------------------------------------

try:

    horizonte = sp.solve(

        sp.Eq(
            f_disforme,
            0,
        ),

        r,

    )

except Exception:

    horizonte = []


# ----------------------------------------------------------
# Tiempo propio
# ----------------------------------------------------------

dt = sp.Symbol("dt")

d_tau_s = sp.sqrt(
    f_schwarzschild
) * dt

d_tau_d = sp.sqrt(
    f_disforme
) * dt


# ----------------------------------------------------------
# Distancia radial propia
# ----------------------------------------------------------

dr = sp.Symbol("dr")

dl_s = dr / sp.sqrt(
    f_schwarzschild
)

dl_d = dr / sp.sqrt(
    f_disforme
)


# ----------------------------------------------------------
# Reporte simbólico
# ----------------------------------------------------------

print()

print("=" * 60)

print("VALIDACIÓN EXPERIMENTAL 1.0")

print("=" * 60)

print()

print("Schwarzschild:")

print(sp.simplify(
    f_schwarzschild
))

print()

print("Disforme:")

print(sp.simplify(
    f_disforme
))

print()

print("Radio mínimo:")

print(
    sp.simplify(
        r_min
    )
)

print()

print("Horizonte:")

print(horizonte)

print()

print("Tiempo propio:")

print(d_tau_d)

print()

print("Distancia propia:")

print(dl_d)

print()


# ----------------------------------------------------------
# Comparación numérica
# ----------------------------------------------------------

M_sol = 1.98847e30

rp = float(
    l_planck.subs({})
)

mp = float(
    m_planck.subs({})
)

rs_num = (
    2 * G * M_sol / c**2
)

rmin_num = (
    rp * (M_sol / mp) ** (1 / 3)
)

R = np.linspace(

    rs_num * 1.01,

    rs_num * 10,

    600,

)

sigma_num = 1 / (
    1
    +
    (rmin_num / R) ** 4
)

fsch = 1 - rs_num / R

fdis = 1 - (rs_num / R) * sigma_num


# ----------------------------------------------------------
# Gráfica
# ----------------------------------------------------------

plt.figure(figsize=(8,6))

plt.plot(
    R,
    fsch,
    label="Schwarzschild",
)

plt.plot(
    R,
    fdis,
    "--",
    label="Atracción Energética",
)

plt.xlabel("r (m)")

plt.ylabel("f(r)")

plt.grid(True)

plt.legend()

plt.tight_layout()

plt.savefig(
    "comparacion_metricas.png",
    dpi=300,
)

plt.show()

print()

print("Figura guardada como:")

print("comparacion_metricas.png")

print()