"""
experimento_schwarzschild.py
============================

Suite Oficial de Validación Experimental del
Motor de Geometría Diferencial (MGD).

Experimento 1

Comparación entre

    • Schwarzschild
    • Atracción Energética

Proyecto:
    Motor de Geometría Diferencial

Autor:
    Elvis Omar Nazario Espinoza
"""

import numpy as np
import matplotlib.pyplot as plt

from fisica.atraccion_energetica import (

    radio_schwarzschild,

    energia_critica_vacio,

    densidad_energia,

    sigma,

    f_schwarzschild,

    f_disforme,

)


# ---------------------------------------------------------
# Parámetros
# ---------------------------------------------------------

MASA_SOL = 1.98847e30

PUNTOS = 200

FACTOR_MIN = 0.50

FACTOR_MAX = 3.00


# ---------------------------------------------------------
# Experimento
# ---------------------------------------------------------

def ejecutar():

    print()

    print("=" * 65)

    print("MGD EXPERIMENTAL SUITE")

    print("EXPERIMENTO 1")

    print("SCHWARZSCHILD vs ATRACCIÓN ENERGÉTICA")

    print("=" * 65)

    print()

    masa = MASA_SOL

    rs = radio_schwarzschild(masa)

    r = np.linspace(

        rs * FACTOR_MIN,

        rs * FACTOR_MAX,

        PUNTOS,

    )

    print("[INFO] Calculando perfiles métricos...")

    f_s = f_schwarzschild(

        masa,

        r,

    )

    f_d = f_disforme(

        masa,

        r,

    )

    print("[INFO] Calculando energía del campo...")

    energia = densidad_energia(

        masa,

        r,

    )

    k0 = energia_critica_vacio()

    sigma_vacio = sigma(

        masa,

        r,

    )

    print("[INFO] Generando gráficas...")

    fig, ax1 = plt.subplots(

        figsize=(10, 6)

    )

    ax1.plot(

        r,

        f_s,

        "--",

        linewidth=2,

        label="Schwarzschild",

    )

    ax1.plot(

        r,

        f_d,

        linewidth=2,

        label="Atracción Energética",

    )

    ax1.axvline(

        rs,

        color="red",

        linestyle=":",

        label="Horizonte clásico",

    )

    ax1.set_xlabel(

        "Distancia radial (m)"

    )

    ax1.set_ylabel(

        "f(r)"

    )

    ax1.grid(True)

    ax2 = ax1.twinx()

    ax2.plot(

        r,

        sigma_vacio,

        color="orange",

        linestyle="-.",

        label="σ(E)",

    )

    ax2.set_ylabel(

        "Respuesta del vacío"

    )

    lineas1, etiquetas1 = ax1.get_legend_handles_labels()

    lineas2, etiquetas2 = ax2.get_legend_handles_labels()

    ax1.legend(

        lineas1 + lineas2,

        etiquetas1 + etiquetas2,

        loc="lower right",

    )

    plt.title(

        "Validación Experimental MGD"

    )

    plt.tight_layout()

    plt.savefig(

        "comparacion_metricas.png",

        dpi=150,

    )

    plt.close(fig)

    print()

    print("RESULTADOS")

    print("-" * 50)

    print(f"Radio Schwarzschild : {rs:.6e} m")

    print(f"Energía crítica K0 : {k0:.6e}")

    print(f"σ(E) mínimo        : {sigma_vacio.min():.6e}")

    print(f"σ(E) máximo        : {sigma_vacio.max():.6e}")

    print()

    print("[OK] Figura generada:")

    print("comparacion_metricas.png")

    print()

    print("=" * 65)

    print("Experimento finalizado correctamente.")

    print("=" * 65)

    print()


# ---------------------------------------------------------
# Punto de entrada
# ---------------------------------------------------------

if __name__ == "__main__":

    ejecutar()