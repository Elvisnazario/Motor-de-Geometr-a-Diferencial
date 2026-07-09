import sympy as sp
import time
import os
# Importamos el núcleo matemático que ya homologamos en la terminal
from motor_tensorial import crear_motor_tensorial, matriz_es_cero

def calcular_tensores_atraccion_energetica():
    print("=======================================================================")
    print(" EXTRACCIÓN DE TENSORES: TEORÍA DE LA ATRACCIÓN ENERGÉTICA")
    print(" Autor: Elvis Omar Nazario Espinoza")
    print("=======================================================================\n")
    
    # 1. DEFINICIÓN DE COORDENADAS Y CONSTANTES FUNDAMENTALES (G, c, M)
    t, r, theta, phi = sp.symbols('t r theta phi')
    G, c, M = sp.symbols('G c M', positive=True)
    K_0 = sp.symbols('K_0', positive=True) # Rigidez crítica del vacío (Escala de Planck)
    
    # 2. DEFINICIÓN DE FUNCIONES SIMBÓLICAS ABSTRACTAS DE LA TEORÍA
    # Para obtener las ecuaciones maestras sin asumir perfiles arbitrarios
    A_int = sp.Function('A_int')(r)  # Geometría modificada interior (sector temporal)
    B_int = sp.Function('B_int')(r)  # Geometría modificada interior (sector radial)
    T_00 = sp.Function('T_00')(r)    # Densidad de energía local efectiva
    
    # Componentes estándar de la solución exterior de Schwarzschild
    A_sch = 1 - (2 * G * M) / (c**2 * r)
    B_sch = 1 / A_sch
    
    # Postulado P2: Parámetro de acoplamiento y Función de respuesta constitutiva sigma(E)
    eta_r = T_00 / K_0
    sigma_E = 1 / (1 + eta_r)
    
    # Construcción Matemática de la Métrica Unificada (Ecuaciones 3 y 4 del artículo)
    g_00_teoria = -(sigma_E * A_sch + (1 - sigma_E) * A_int)
    g_11_teoria = sigma_E * B_sch + (1 - sigma_E) * B_int
    g_22_teoria = r**2
    g_33_teoria = r**2 * sp.sin(theta)**2
    
    metric_components = [g_00_teoria, g_11_teoria, g_22_teoria, g_33_teoria]
    
    print("[PROCESO] Enviando la métrica elástica unificada al Motor MGD...")
    print(f"-> g_00 = {g_00_teoria}")
    print(f"-> g_11 = {g_11_teoria}\n")
    
    # 3. COMPUTACIÓN SIMBÓLICA MEDIANTE EL MOTOR HOMOLOGADO
    t_inicio = time.time()
    G_tensor, R_escalar = crear_motor_tensorial(metric_components, nombre_cache="cache_atraccion_energetica")
    t_fin = time.time()
    
    print(f"✔ [ÉXITO] Ecuaciones analíticas calculadas en {t_fin - t_inicio:.2f} segundos.")
    print("=======================================================================")
    print(" TENSORES RESULTANTES DE LA ATRACCIÓN ENERGÉTICA (G_mu_nu):")
    print("=======================================================================")
    
    # Limpieza algebraica estricta elemento a elemento
    G_limpio = G_tensor.applyfunc(lambda x: sp.together(sp.cancel(x)))
    R_limpio = sp.together(sp.cancel(R_escalar))
    
    print("\n> ESCALAR DE RICCI ACTIVO (R):")
    sp.pprint(R_limpio)
    print("\n" + "-"*70)
    
    print("\n> COMPONENTE TEMPORAL G[0,0] (Respuesta a la Densidad de Energía):")
    sp.pprint(G_limpio[0,0])
    print("\n" + "-"*70)
    
    print("\n> COMPONENTE RADIAL G[1,1] (Tensión Elástica del Vacío):")
    sp.pprint(G_limpio[1,1])
    print("\n" + "-"*70)
    
    print("\n> COMPONENTE ANGULAR G[2,2]:")
    sp.pprint(G_limpio[2,2])
    print("\n=======================================================================")
    
    # 4. PRUEBA DE LÍMITE ASINTÓTICO (Sección 3.1 del artículo)
    # Evaluamos el límite cuando T_00 -> 0, lo que implica eta -> 0 y sigma -> 1
    print("\n[PROBANDO LÍMITE MACROSCÓPICO ASINTÓTICO T_00 -> 0]:")
    G_limite_vacio = G_limpio.subs(T_00, 0)
    
    if matriz_es_cero(G_limite_vacio):
        print("✔ EXCELENTE: Cuando T_00 -> 0, la función constitutiva satura en 1")
        print("   y su teoría recupera matemáticamente el vacío puro de Schwarzschild (G_mu_nu = 0).")
    else:
        print("⚠ ATVERTENCIA: Se detectaron residuos algebraicos en el límite macroscópico.")
    print("=======================================================================")

if __name__ == "__main__":
    calcular_tensores_atraccion_energetica()
