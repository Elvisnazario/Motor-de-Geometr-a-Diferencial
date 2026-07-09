import sympy as sp
import time
from motor_tensorial import crear_motor_tensorial, matriz_es_cero

def calcular_tensores_atraccion_energetica_optimizado():
    print("=======================================================================")
    print(" EXTRACCIÓN DE TENSORES (OPTIMIZADA) - ATRACCIÓN ENERGÉTICA")
    print(" Autor: Elvis Omar Nazario Espinoza")
    print("=======================================================================\n")
    
    t, r, theta, phi = sp.symbols('t r theta phi')
    G, c, M = sp.symbols('G c M', positive=True)
    K_0 = sp.symbols('K_0', positive=True)
    
    # Declaramos las funciones base de su modelo
    A_int = sp.Function('A_int')(r)
    B_int = sp.Function('B_int')(r)
    T_00 = sp.Function('T_00')(r)
    
    # Optimizador analítico: Tratamos a sigma como una función directa de r para evitar
    # que SymPy expanda fracciones masivas con denominadores anidados durante las derivadas.
    sigma = sp.Function('sigma')(r)
    
    A_sch = 1 - (2 * G * M) / (c**2 * r)
    B_sch = 1 / A_sch
    
    # Métrica parametrizada elásticamente
    g_00_calc = -(sigma * A_sch + (1 - sigma) * A_int)
    g_11_calc = sigma * B_sch + (1 - sigma) * B_int
    g_22_calc = r**2
    g_33_calc = r**2 * sp.sin(theta)**2
    
    metric_components = [g_00_calc, g_11_calc, g_22_calc, g_33_calc]
    
    print("[PROCESO] Calculando tensores analíticos con desacoplamiento funcional...")
    t_inicio = time.time()
    G_tensor, R_escalar = crear_motor_tensorial(metric_components, nombre_cache=None)
    
    # Postulado P2: Definición real de la función constitutiva de su artículo
    eta_real = T_00 / K_0
    sigma_real = 1 / (1 + eta_real)
    
    # SUSTITUCIÓN DE RETORNO: Inyectamos la definición exacta de su teoría (P2)
    # Reemplazamos la función sigma(r) y su derivada por las expresiones explícitas de su ecuación.
    G_final = G_tensor.subs(sigma, sigma_real)
    R_final = R_escalar.subs(sigma, sigma_real)
    
    # Reemplazamos también la derivada formal de sigma(r) respecto a r
    derivada_sigma_real = sp.diff(sigma_real, r)
    G_final = G_final.subs(sp.diff(sigma, r), derivada_sigma_real)
    R_final = R_final.subs(sp.diff(sigma, r), derivada_sigma_real)
    
    t_fin = time.time()
    print(f"✔ [ÉXITO] Ecuaciones analíticas extraídas en {t_fin - t_inicio:.2f} segundos.")
    print("=======================================================================")
    print(" TENSORES FINALES INYECTADOS CON LA FUNCIÓN CONSTITUTIVA POSTULADA:")
    print("=======================================================================")
    
    # Aplicamos simplificación directa y limpia
    G_limpio = G_final.applyfunc(lambda x: sp.together(sp.cancel(x)))
    R_limpio = sp.together(sp.cancel(R_final))
    
    print("\n> ESCALAR DE RICCI ACTIVO (R):")
    sp.pprint(R_limpio)
    print("\n" + "-"*70)
    
    print("\n> COMPONENTE TEMPORAL G[0,0] (Densidad de Energía Efectiva):")
    sp.pprint(G_limpio[0,0])
    print("\n" + "-"*70)
    
    print("\n> COMPONENTE RADIAL G[1,1] (Tensión Elástica):")
    sp.pprint(G_limpio[1,1])
    print("\n=======================================================================")
    
    print("\n[PROBANDO LÍMITE ASINTÓTICO T_00 -> 0]:")
    G_limite_vacio = G_limpio.subs(T_00, 0)
    if matriz_es_cero(G_limite_vacio):
        print("✔ EXCELENTE: El límite T_00 -> 0 recupera perfectamente Schwarzschild.")
    else:
        print("⚠ ADVERTENCIA: Se encontraron residuos en el límite asintótico.")

if __name__ == "__main__":
    calcular_tensores_atraccion_energetica_optimizado()
