import numpy as np
import matplotlib.pyplot as plt
import json

# Importación nativa del motor para garantizar coherencia analítica absoluta
from fisica.atraccion_energetica import f_disforme, radio_schwarzschild, radio_minimo

# ==============================================================================
# MÓDULO PREDICTIVO MGD: COLAPSO GRAVITACIONAL DETENIDO (V1.9 ARQUITECTURA FINAL)
# ==============================================================================

def generar_perfil_vacio_modelo(r, r_min, m_nucleo, alpha):
    """
    Calcula el andamiaje fenomenológico provisional para el vacío disforme.
    NOTA: Bloque temporal a sustituir por soluciones exactas de las ecuaciones MGD.
    """
    rho_vac_modelo = (m_nucleo / (4 * np.pi * r**3)) * np.exp(-r_min / r)
    P_eff_modelo = alpha * (rho_vac_modelo**2)
    return rho_vac_modelo, P_eff_modelo


def calcular_fuerzas_dinamicas(r, f, rho_vac_modelo, P_eff_modelo):
    """
    Calcula el balance de fuerzas por unidad de volumen usando edge_order=2.
    """
    df_dr = np.gradient(f, r, edge_order=2)
    a_grav = -0.5 * df_dr  
    F_grav = rho_vac_modelo * a_grav

    grad_P = np.gradient(P_eff_modelo, r, edge_order=2)
    eps = np.finfo(float).eps
    F_presion = - grad_P / np.maximum(rho_vac_modelo, eps) 

    F_neta = F_grav + F_presion
    
    if not np.all(np.isfinite(F_neta)):
        raise ValueError("[MGD ERROR CRÍTICO] F_neta contiene valores NaN o Inf no finitos.")
        
    return F_grav, F_presion, F_neta


def buscar_punto_equilibrio(r, F_neta):
    """
    Identifica el índice de minimización de fuerza y extrae sus propiedades locales.
    """
    idx_eq = np.nanargmin(np.abs(F_neta))
    r_equilibrio = r[idx_eq]
    residuo_eq = np.abs(F_neta[idx_eq])

    dF_dr = np.gradient(F_neta, r, edge_order=2)
    pendiente_eq = dF_dr[idx_eq]
    
    return idx_eq, r_equilibrio, residuo_eq, pendiente_eq


def exportar_datos_simulacion(archivo_csv, r, rho, p, f_g, f_p, f_n):
    """
    Guarda la matriz numérica cruda para garantizar reproducibilidad científica.
    """
    datos = np.column_stack((r, rho, p, f_g, f_p, f_n))
    np.savetxt(
        archivo_csv,
        datos,
        delimiter=",",
        header="r,rho_vac_modelo,P_eff_modelo,F_grav,F_presion,F_neta",
        comments=""
    )


def generar_grafica_diagnostico(archivo_png, r, F_grav, F_presion, r_equilibrio, idx_eq):
    """
    Genera el entregable visual del perfil de fuerzas marcando explícitamente el punto de equilibrio.
    """
    plt.figure(figsize=(10, 6))
    plt.loglog(r, np.abs(F_grav), 'r--', label='Gravedad Relatividad MGD (-0.5*df/dr)')
    plt.loglog(r, np.abs(F_presion), 'b-', label='Presión Constitutiva Vacío Disforme')
    
    # Marcador exacto del punto donde se intersectan/minimizan las fuerzas
    plt.scatter(r_equilibrio, np.abs(F_grav[idx_eq]), s=50, color='black', zorder=5, label='Punto de Equilibrio')
    plt.axvline(r_equilibrio, color='g', linestyle=':', label=f'R_eq ({r_equilibrio:.2e} m)')
    
    plt.title('Colapso Gravitacional Detenido (Módulo Predictivo MGD v1.9)')
    plt.xlabel('Radio r (m)')
    plt.ylabel('Fuerza por Unidad de Volumen')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.savefig(archivo_png, dpi=150)
    plt.close()


def ejecutar_pipeline_colapso(masa=1.98847e30, alpha=1.0e-45, resolucion=4000):
    """
    Pipeline silencioso de cálculo y exportación de datos. No escribe en stdout.
    """
    r_min = radio_minimo(masa)
    r_s = radio_schwarzschild(masa)
    
    r = np.logspace(np.log10(r_min), np.log10(r_s * 1.5), resolucion)
    
    f = f_disforme(masa, r)
    rho_vac_m, P_eff_m = generar_perfil_vacio_modelo(r, r_min, masa, alpha)
    F_g, F_p, F_n = calcular_fuerzas_dinamicas(r, f, rho_vac_m, P_eff_m)
    idx_eq, r_eq, res_eq, pend_eq = buscar_punto_equilibrio(r, F_n)
    
    # Generación dinámica de rutas basadas en los parámetros de entrada
    nombre_base = f"colapso_m_{masa:.1e}_alpha_{alpha:.1e}".replace("+", "")
    csv_path = f"{nombre_base}.csv"
    png_path = f"{nombre_base}.png"
    json_path = f"{nombre_base}.json"
    
    # Persistencia en disco
    exportar_datos_simulacion(csv_path, r, rho_vac_m, P_eff_m, F_g, F_p, F_n)
    generar_grafica_diagnostico(png_path, r, F_g, F_p, r_eq, idx_eq)
    
    config = {
        "masa": masa,
        "alpha": alpha,
        "resolucion": resolucion,
        "r_min": r_min,
        "r_s": r_s,
        "archivos_generados": [csv_path, png_path, json_path]
    }
    with open(json_path, "w") as f_json:
        json.dump(config, f_json, indent=4)
    
    return {
        "r_eq": r_eq,
        "r_eq_rs_ratio": r_eq / r_s,
        "residuo": res_eq,
        "pendiente": pend_eq,
        "r_s": r_s,
        "r_min": r_min,
        "alpha": alpha,
        "masa": masa,
        "csv_path": csv_path,
        "png_path": png_path,
        "json_path": json_path
    }


if __name__ == "__main__":
    # La interfaz maneja de forma exclusiva la salida en pantalla y el control de excepciones
    try:
        M_SOLAR = 1.98847e30
        ALPHA_TEST = 1.0e-45
        
        res = ejecutar_pipeline_colapso(masa=M_SOLAR, alpha=ALPHA_TEST)
        
        print(f"[MGD ACCIÓN] Pipeline v1.9 finalizado con éxito.")
        print(f" -> Parámetros: Masa = {res['masa']:.2e} kg | Alpha = {res['alpha']:.2e}")
        print(f" -> Radio de equilibrio (R_eq): {res['r_eq']:.6e} m")
        print(f" -> Relación con el Horizonte (R_eq / r_s): {res['r_eq_rs_ratio']:.6e}")
        print(f" -> Residuo de fuerza en R_eq: {res['residuo']:.6e}")
        print(f" -> Estabilidad dF/dr en R_eq: {res['pendiente']:.6e}")
        print(f" -> Archivos exportados con éxito:\n    - {res['csv_path']}\n    - {res['png_path']}\n    - {res['json_path']}")

    except Exception as e:
        print(f"\n[MGD ERROR DE EJECUCIÓN] Ocurrió una anomalía en el pipeline: {e}")
        raise
