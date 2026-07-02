"""
experimentos/experimento_horizonte.py
=====================================

EXPERIMENTO 4
Horizonte efectivo del Motor de Geometría Diferencial.

Objetivos
---------
1) Localizar raíces de f_MGD(r)=0
2) Comparar con el horizonte clásico de Schwarzschild
3) Calcular el potencial efectivo
4) Visualizar la diferencia entre ambos modelos

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import brentq

from fisica.atraccion_energetica import (

    radio_schwarzschild,
    radio_minimo,

    f_schwarzschild,
    f_disforme,

)

# -------------------------------------------------------
# CONFIGURACIÓN
# -------------------------------------------------------

MASA = 1.98847e30

PUNTOS = 4000

FACTOR_MAX = 20.0

L = 4.5e12

# -------------------------------------------------------
# Potencial efectivo
# -------------------------------------------------------

def potencial_schwarzschild(masa, r):

    f = f_schwarzschild(

        masa,

        r,

    )

    return f * (

        1.0 +

        (L**2)/(r**2)

    )


def potencial_disforme(masa, r):

    f = f_disforme(

        masa,

        r,

    )

    return f * (

        1.0 +

        (L**2)/(r**2)

    )

# -------------------------------------------------------
# Buscador de raíces
# -------------------------------------------------------

def buscar_raices(x, y):

    raices = []

    for i in range(len(x)-1):

        if np.isnan(y[i]) or np.isnan(y[i+1]):
            continue

        if y[i] == 0:

            raices.append(x[i])

            continue

        if y[i]*y[i+1] < 0:

            try:

                raiz = brentq(

                    lambda r: f_disforme(MASA, r),

                    x[i],

                    x[i+1],

                )

                raices.append(raiz)

            except Exception:

                pass

    return raices

# -------------------------------------------------------
# Experimento
# -------------------------------------------------------

def ejecutar():

    print()

    print("="*70)
    print("EXPERIMENTO 4")
    print("HORIZONTE EFECTIVO")
    print("="*70)

    rs = radio_schwarzschild(MASA)

    rmin = radio_minimo(MASA)

    r = np.logspace(

        np.log10(rmin),

        np.log10(rs*FACTOR_MAX),

        PUNTOS,

    )

    print("[INFO] Calculando métricas...")

    fGR = f_schwarzschild(

        MASA,

        r,

    )

    fMGD = f_disforme(

        MASA,

        r,

    )

    print("[INFO] Buscando raíces...")

    raices = buscar_raices(

        r,

        fMGD,

    )

    print("[INFO] Calculando potenciales...")

    Vgr = potencial_schwarzschild(

        MASA,

        r,

    )

    Vmgd = potencial_disforme(

        MASA,

        r,

    )

    print("[INFO] Generando gráficas...")

    fig, ax = plt.subplots(

        figsize=(10,6)

    )

    ax.semilogx(

        r,

        fGR,

        "--",

        linewidth=2,

        label="Schwarzschild",

    )

    ax.semilogx(

        r,

        fMGD,

        linewidth=2,

        label="MGD",

    )

    ax.axhline(

        0,

        color="black",

        linewidth=1,

    )

    ax.axvline(

        rs,

        color="red",

        linestyle=":",

        label="Horizonte clásico",

    )

    for raiz in raices:

        ax.axvline(

            raiz,

            color="green",

            linestyle="--",

            linewidth=2,

        )

    ax.grid(

        True,

        which="both",

        linestyle=":",

    )

    ax.set_xlabel("r (m)")
    ax.set_ylabel("f(r)")
    ax.legend()

    plt.tight_layout()

    plt.savefig(

        "horizonte_mgd.png",

        dpi=150,

    )

    plt.close()

    plt.figure(

        figsize=(10,6)

    )

    plt.semilogx(

        r,

        Vgr,

        "--",

        linewidth=2,

        label="Vef Schwarzschild",

    )

    plt.semilogx(

        r,

        Vmgd,

        linewidth=2,

        label="Vef MGD",

    )

    plt.grid(

        True,

        which="both",

        linestyle=":",

    )

    plt.xlabel("r (m)")
    plt.ylabel("Vef(r)")
    plt.legend()

    plt.tight_layout()

    plt.savefig(

        "potencial_efectivo_mgd.png",

        dpi=150,

    )

    plt.close()

    print()

    print("="*70)
    print("RESULTADOS")
    print("="*70)

    print(f"Radio de Schwarzschild : {rs:.6e} m")
    print(f"Radio mínimo MGD       : {rmin:.6e} m")

    if len(raices)==0:

        print()
        print("No se encontró horizonte efectivo.")
        print("La métrica MGD no cambia de signo.")

    else:

        print()
        print("Horizontes encontrados:")

        for i,raiz in enumerate(raices):

            print(f"{i+1}: {raiz:.6e} m")

    print()
    print("[OK] Archivos generados:")
    print("    horizonte_mgd.png")
    print("    potencial_efectivo_mgd.png")
    print()
    print("="*70)
    print("Experimento finalizado.")
    print("="*70)
    print()


if __name__ == "__main__":

    ejecutar()