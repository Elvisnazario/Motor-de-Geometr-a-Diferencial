"""
geodesicas_disformes.py
=======================

Experimento 2 del Motor de Geometría Diferencial.

Compara las trayectorias geodésicas entre:

    • Schwarzschild
    • Atracción Energética

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp

from fisica.atraccion_energetica import (
    f_schwarzschild,
    f_disforme,
)

# ----------------------------------------------------------
# Parámetros
# ----------------------------------------------------------

MASA_SOL = 1.98847e30

R_INICIAL = 5.0e4

V_RADIAL = 0.0

L = 4.5e12

TIEMPO_FINAL = 6000

# ----------------------------------------------------------
# Geodésica Schwarzschild
# ----------------------------------------------------------

def geodesica_schwarzschild(
    t,
    y,
):

    r, vr = y

    f = f_schwarzschild(
        MASA_SOL,
        r,
    )

    ar = -0.5 * np.gradient(
        np.array([f]),
    )[0]

    return [

        vr,

        ar,

    ]

# ----------------------------------------------------------
# Geodésica Disforme
# ----------------------------------------------------------

def geodesica_disforme(
    t,
    y,
):

    r, vr = y

    f = f_disforme(
        MASA_SOL,
        r,
    )

    ar = -0.5 * np.gradient(
        np.array([f]),
    )[0]

    return [

        vr,

        ar,

    ]

# ----------------------------------------------------------
# Simulación
# ----------------------------------------------------------

def ejecutar():

    print()

    print("=" * 60)

    print("EXPERIMENTO 2")

    print("GEODÉSICAS")

    print("=" * 60)

    y0 = [

        R_INICIAL,

        V_RADIAL,

    ]

    sol_s = solve_ivp(

        geodesica_schwarzschild,

        (0, TIEMPO_FINAL),

        y0,

        max_step=5,

    )

    sol_d = solve_ivp(

        geodesica_disforme,

        (0, TIEMPO_FINAL),

        y0,

        max_step=5,

    )

    plt.figure(figsize=(8,6))

    plt.plot(

        sol_s.t,

        sol_s.y[0],

        "--",

        label="Schwarzschild",

    )

    plt.plot(

        sol_d.t,

        sol_d.y[0],

        linewidth=2,

        label="Atracción Energética",

    )

    plt.xlabel("Tiempo")

    plt.ylabel("Radio")

    plt.grid(True)

    plt.legend()

    plt.tight_layout()

    plt.savefig(

        "geodesicas_disformes.png",

        dpi=150,

    )

    plt.close()

    print()

    print("Archivo generado:")

    print("geodesicas_disformes.png")

    print()

if __name__ == "__main__":

    ejecutar()