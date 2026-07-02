"""
experimentos/experimento_transicion.py
======================================

Experimento 3 del Motor de Geometría Diferencial (MGD)

OBJETIVO
--------
Buscar el umbral de transición entre el régimen clásico
(Relatividad General) y el régimen donde la Atracción
Energética comienza a modificar la métrica.

Este experimento NO integra geodésicas.

Su objetivo consiste en barrer masas durante muchos
órdenes de magnitud y medir:

    η(r) = T00 / K0

    σ(E)

    Δ(r)=|f_MGD-f_GR|

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import numpy as np
import matplotlib.pyplot as plt

from fisica.atraccion_energetica import (
    radio_schwarzschild,
    radio_minimo,
    energia_critica_vacio,
    densidad_energia,
    sigma,
    f_schwarzschild,
    f_disforme,
)

# ----------------------------------------------------------
# CONFIGURACIÓN
# ----------------------------------------------------------

MASA_MIN = 1e-12
MASA_MAX = 1e12

MUESTRAS_MASA = 120
MUESTRAS_RADIALES = 250
FACTOR_RADIAL_MAX = 100.0

# ----------------------------------------------------------
# EXPERIMENTO
# ----------------------------------------------------------

def ejecutar():

    print()
    print("=" * 70)
    print("MGD EXPERIMENTAL SUITE")
    print("EXPERIMENTO 3")
    print("TRANSICIÓN DE RÉGIMEN")
    print("=" * 70)
    print()

    masas = np.logspace(
        np.log10(MASA_MIN),
        np.log10(MASA_MAX),
        MUESTRAS_MASA,
    )
    eta_max = []
    sigma_min = []
    delta_max = []
    radios_minimos = []
    radios_schwarzschild = []

    K0 = energia_critica_vacio()

    print("[INFO] Iniciando barrido logarítmico de masas...")

    for i, masa in enumerate(masas):

        rs = radio_schwarzschild(masa)
        rmin = radio_minimo(masa)

        radios_schwarzschild.append(rs)
        radios_minimos.append(rmin)

        # Mallado logarítmico desde el radio mínimo
        r = np.logspace(
            np.log10(rmin),
            np.log10(max(rs * FACTOR_RADIAL_MAX, rmin * 10.0)),
            MUESTRAS_RADIALES,
        )

        T00 = densidad_energia(masa, r)
        eta = T00 / K0
        sigma_vacio = sigma(masa, r)
        f_gr = f_schwarzschild(masa, r)
        f_mgd = f_disforme(masa, r)

        delta = np.abs(f_mgd - f_gr)

        eta_max.append(np.max(eta))
        sigma_min.append(np.min(sigma_vacio))
        delta_max.append(np.max(delta))

        if (i + 1) % 20 == 0:
            print(f"[{i+1:03d}/{MUESTRAS_MASA}] masas procesadas...")

    eta_max = np.asarray(eta_max)
    sigma_min = np.asarray(sigma_min)
    delta_max = np.asarray(delta_max)

    print()
    print("[INFO] Generando gráficas...")

    # --------------------------------------------------
    # Figura 1 - η
    # --------------------------------------------------
    plt.figure(figsize=(9,6))
    plt.loglog(masas, eta_max, linewidth=2, color="tab:blue")
    plt.grid(True, which="both", linestyle=":")
    plt.xlabel("Masa (kg)")
    plt.ylabel("η = T00 / K0")
    plt.title("Transición energética del vacío")
    plt.tight_layout()
    plt.savefig("transicion_eta.png", dpi=150)
    plt.close()

    # --------------------------------------------------
    # Figura 2 - σ(E)
    # --------------------------------------------------
    plt.figure(figsize=(9,6))
    plt.semilogx(masas, sigma_min, linewidth=2, color="tab:orange")
    plt.grid(True, which="both", linestyle=":")
    plt.xlabel("Masa (kg)")
    plt.ylabel("σ(E) mínimo")
    plt.title("Respuesta constitutiva del vacío")
    plt.tight_layout()
    plt.savefig("transicion_sigma.png", dpi=150)
    plt.close()

    # --------------------------------------------------
    # Figura 3 - Diferencia entre métricas
    # --------------------------------------------------
    plt.figure(figsize=(9,6))
    plt.loglog(masas, delta_max, linewidth=2, color="tab:red")
    plt.grid(True, which="both", linestyle=":")
    plt.xlabel("Masa (kg)")
    plt.ylabel("Δ = |fMGD − fGR|")
    plt.title("Separación entre Schwarzschild y Atracción Energética")
    plt.tight_layout()
    plt.savefig("transicion_delta.png", dpi=150)
    plt.close()

    # --------------------------------------------------
    # Exportación CSV
    # --------------------------------------------------
    datos = np.column_stack((
        masas,
        radios_schwarzschild,
        radios_minimos,
        eta_max,
        sigma_min,
        delta_max,
    ))

    np.savetxt(
        "transicion_resultados.csv",
        datos,
        delimiter=",",
        header="masa,radio_schwarzschild,radio_minimo,eta_max,sigma_min,delta_max",
        comments="",
    )

    # --------------------------------------------------
    # Resumen Global
    # --------------------------------------------------
    indice = np.argmax(delta_max)

    print()
    print("=" * 70)
    print("RESULTADOS GLOBALES DEL BARRIDO")
    print("=" * 70)
    print(f"Masa crítica aproximada : {masas[indice]:.6e} kg")
    print(f"Radio mínimo           : {radios_minimos[indice]:.6e} m")
    print(f"η máximo               : {eta_max[indice]:.6e}")
    print(f"σ mínimo               : {sigma_min[indice]:.6e}")
    print(f"Δ máximo               : {delta_max[indice]:.6e}")

    print()
    print("[OK] Archivos generados:")
    print("    transicion_eta.png")
    print("    transicion_sigma.png")
    print("    transicion_delta.png")
    print("    transicion_resultados.csv")
    print()

    # --------------------------------------------------
    # Diagnóstico de Transición (Ahora dentro de ejecutar)
    # --------------------------------------------------
    print("[INFO] Buscando primer punto de transición...")

    umbral_sigma = 0.999999
    indices = np.where(sigma_min < umbral_sigma)[0]

    print()
    if len(indices) == 0:
        print("No se detectó transición física.")
        print("La teoría permanece en el régimen clásico")
        print("para todo el intervalo de masas analizado.")
    else:
        i = indices[0]
        print("=" * 70)
        print("¡¡ PRIMERA TRANSICIÓN DETECTADA !!")
        print("=" * 70)
        print(f"Masa                 : {masas[i]:.6e} kg")
        print(f"Radio Schwarzschild  : {radios_schwarzschild[i]:.6e} m")
        print(f"Radio mínimo         : {radios_minimos[i]:.6e} m")
        print(f"η                    : {eta_max[i]:.6e}")
        print(f"σ(E)                 : {sigma_min[i]:.12f}")
        print(f"Δ                    : {delta_max[i]:.6e}")

    print()
    print("=" * 70)
    print("INTERPRETACIÓN FÍSICA")
    print("=" * 70)
    print(
        "Este experimento determina el régimen donde la\n"
        "rigidez constitutiva del vacío comienza a modificar\n"
        "la solución clásica de Schwarzschild."
    )
    print()
    print(
        "Si η << 1 entonces el MGD reproduce exactamente\n"
        "la Relatividad General (Principio de Correspondencia)."
    )
    print()
    print(
        "Cuando η≈1 la función σ(E) empieza a disminuir\n"
        "y aparecen las primeras diferencias geométricas."
    )
    print()
    print(
        "La gráfica Δ(r) constituye el indicador principal\n"
        "para localizar el umbral de transición."
    )
    print()
    print("=" * 70)
    print("EXPERIMENTO FINALIZADO")
    print("=" * 70)
    print()

# ----------------------------------------------------------
# Punto de entrada
# ----------------------------------------------------------
if __name__ == "__main__":
    ejecutar()
