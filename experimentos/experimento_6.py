import numpy as np
import json
import os

# Corrección analítica de la ruta: importamos desde el nombre real de tu repositorio
from experimentos.colapso_detenido import ejecutar_pipeline_colapso

# ==============================================================================
# EXPERIMENTO 6: TEST DE ROBUSTEZ GLOBAL Y ESTABILIDAD PARAMÉTRICA (MGD)
# ==============================================================================

def ejecutar_barrido_robustez():
    """
    Lanza un barrido automático multiescala combinando diferentes masas y
    órdenes de magnitud de alpha para auditar la estabilidad de R_eq / r_s.
    """
    # 1. Definición del espacio de parámetros a explorar
    masas_prueba = [
        1.98847e30,  # 1 Masa Solar (Estrella estándar)
        1.0e31,       # ~5 Masas Solares (Estructura colapsada ligera)
        1.0e33        # ~500 Masas Solares (Estructura colapsada intermedia)
    ]
    
    alphas_prueba = [1.0e-46, 5.0e-46, 1.0e-45, 5.0e-45, 1.0e-44]
    
    # Estructura para almacenar la matriz de resultados científicos
    tabla_resultados = []
    
    # 2. Bucle automatizado de simulación masiva
    for m in masas_prueba:
        for a in alphas_prueba:
            try:
                # El pipeline colapso_detenido calcula y exporta silenciosamente
                res = ejecutar_pipeline_colapso(masa=m, alpha=a, resolucion=4000)
                
                # Extraemos los observables clave para el análisis de estabilidad
                tabla_resultados.append({
                    "Masa (kg)": f"{m:.2e}",
                    "Alpha": f"{a:.2e}",
                    "R_eq (m)": f"{res['r_eq']:.4e}",
                    "R_eq / r_s": f"{res['r_eq_rs_ratio']:.6e}",
                    "Pendiente dF/dr": f"{res['pendiente']:.4e}"
                })
            except Exception as e:
                # Si una combinación extrema diverge, se reporta sin detener el bucle
                tabla_resultados.append({
                    "Masa (kg)": f"{m:.2e}",
                    "Alpha": f"{a:.2e}",
                    "R_eq (m)": "DIVERGENCIA",
                    "R_eq / r_s": "NaN",
                    "Pendiente dF/dr": "NaN"
                })

    # 3. Persistencia de la tabla de robustez en disco
    archivo_reporte = "experimento_6_robustez.json"
    with open(archivo_reporte, "w") as f_out:
        json.dump(tabla_resultados, f_out, indent=4)
        
    return archivo_reporte


if __name__ == "__main__":
    print("[MGD ACCIÓN] Iniciando Experimento 6: Barrido automatizado de robustez paramétrica...")
    reporte = ejecutar_barrido_robustez()
    print(f"[MGD ACCIÓN] Matriz de robustez finalizada con éxito. Datos guardados en: '{reporte}'")
