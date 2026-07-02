"""
experimentos/experimento_kretschmann.py
=======================================
EXPERIMENTO 5: Evaluación del Escalar de Kretschmann Esférico.
Compara la curvatura clásica contra el modelo MGD incorporando limpieza nan_to_num.

Proyecto: Motor de Geometría Diferencial (MGD)
Autor: Elvis Omar Nazario Espinoza
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

from fisica.atraccion_energetica import radio_schwarzschild, radio_minimo
from fisica.invariantes import kretschmann_schwarzschild, kretschmann_mgd

# CONFIGURACIÓN DEL EXPERIMENTO
MASA = 1.98847e30  # Una masa solar
PUNTOS = 4000
FACTOR_MAX = 5.0

def ejecutar():
    print()
    print("="*70)
    print("EXPERIMENTO 5: ESCALAR DE KRETSCHMANN (CÁLCULO BLINDADO DEFINITIVO)")
    print("="*70)
    
    rs = radio_schwarzschild(MASA)
    rmin = radio_minimo(MASA)
    
    # Malla logarítmica profunda sin tocar regiones prohibidas
    r = np.logspace(
        np.log10(rmin),
        np.log10(rs * FACTOR_MAX),
        PUNTOS
    )
    
    print("[INFO] Evaluando curvatura de Schwarzschild (K_GR)...")
    K_gr = kretschmann_schwarzschild(MASA, r)
    
    print("[INFO] Evaluando curvatura del Motor MGD (K_MGD)...")
    K_mgd = kretschmann_mgd(MASA, r)
    
    # Analizar el cociente con protección épsilon y limpieza de singularidades artificiales
    print("[INFO] Analizando cociente de curvaturas K_MGD / K_GR...")
    eps = np.finfo(float).eps
    cociente = K_mgd / np.maximum(K_gr, eps)
    cociente = np.nan_to_num(cociente, nan=0.0, posinf=np.inf, neginf=0.0)
    
    # Verificar límites numéricos cerca de r_min
    val_K_gr_min = K_gr[0]
    val_K_mgd_min = K_mgd[0]
    
    print("\n" + "-"*50)
    print("ANÁLISIS DE CERCANÍA AL LÍMITE CENTRAL (r -> r_min):")
    print(f"Radio mínimo evaluado : {rmin:.6e} m")
    print(f"Kretschmann GR (Clásico) en r_min : {val_K_gr_min:.6e} m^-4")
    print(f"Kretschmann MGD (Tu Motor) en r_min : {val_K_mgd_min:.6e} m^-4")
    print("-"*50 + "\n")
    
    # GENERAR GRÁFICA DEL ESCALAR DE KRETSCHMANN
    print("[INFO] Generando gráfica de curvaturas...")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.loglog(r, K_gr, "--", linewidth=2, label="K_GR (Schwarzschild)")
    ax.loglog(r, K_mgd, linewidth=2, label="K_MGD (Motor)")
    
    ax.axvline(rs, color="red", linestyle=":", label="Horizonte clásico")
    ax.axvline(rmin, color="green", linestyle="--", label="Radio Mínimo MGD")
    
    ax.grid(True, which="both", linestyle=":")
    ax.set_xlabel("Radio r (m)")
    ax.set_ylabel("Escalar de Kretschmann K (m^-4)")
    ax.set_title("Comparativa de Curvatura de Kretschmann: Regularidad vs Singularidad")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig("curvatura_kretschmann.png", dpi=150)
    plt.close()
    
    # GENERAR GRÁFICA DEL COCIENTE K_MGD / K_GR
    print("[INFO] Generando gráfica del cociente de atenuación...")
    plt.figure(figsize=(10, 6))
    plt.semilogx(r, cociente, color="purple", linewidth=2, label="K_MGD / K_GR")
    plt.axvline(rs, color="red", linestyle=":")
    plt.grid(True, which="both", linestyle=":")
    plt.xlabel("Radio r (m)")
    plt.ylabel("Cociente K_MGD / K_GR")
    plt.title("Factor de Amortiguación de Curvatura del MGD respecto a Relatividad General")
    plt.legend()
    plt.tight_layout()
    plt.savefig("cociente_curvatura.png", dpi=150)
    plt.close()
    
    # EXPORTAR RESULTADOS A CSV
    print("[INFO] Exportando datos analíticos a CSV...")
    csv_path = "resultados_kretschmann.csv"
    with open(csv_path, mode="w", newline="") as f_csv:
        writer = csv.writer(f_csv)
        writer.writerow(["radio", "K_GR", "K_MGD", "cociente_MGD_GR"])
        for i in range(len(r)):
            writer.writerow([r[i], K_gr[i], K_mgd[i], cociente[i]])
            
    print("\n" + "="*70)
    print("¡EXPERIMENTO 5 COMPLETADO CON ÉXITO!")
    print("Archivos generados en tu suite con corrección de sesgo en borde.")
    print("  -> curvatura_kretschmann.png")
    print("  -> cociente_curvatura.png")
    print("  -> resultados_kretschmann.csv")
    print("="*70 + "\n")

if __name__ == "__main__":
    ejecutar()
