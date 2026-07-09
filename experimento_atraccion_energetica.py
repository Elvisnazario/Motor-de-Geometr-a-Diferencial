import sympy as sp
import time
from motor_tensorial import crear_motor_tensorial, matriz_es_cero

def calcular_tensores_rigor_analitico():
    print("=======================================================================")
    print(" EXTRACCIÓN RIGUROSA DE TENSORES - TEORÍA DE LA ATRACCIÓN ENERGÉTICA")
    print(" Autor: Elvis Omar Nazario Espinoza")
    print("=======================================================================\n")
    
    t, r, theta, phi = sp.symbols('t r theta phi')
    G, c, M = sp.symbols('G c M', positive=True)
    K_0 = sp.symbols('K_0', positive=True)
    
    # 1. FUNCIONES ABSTRACTAS GENERALES PARA EL MOTOR (Máximo Rigor Geométrico)
    A = sp.Function('A')(r)
    B = sp.Function('B')(r)
    
    # Métrica exacta en su forma covariante pura
    metric_components = [-A, B, r**2, r**2 * sp.sin(theta)**2]
    
    print("[PROCESO] Calculando el Tensor de Einstein formal exacto...")
    t_inicio = time.time()
    G_tensor, R_escalar = crear_motor_tensorial(metric_components, nombre_cache=None)
    t_fin = time.time()
    print(f"✔ [ÉXITO] Estructura geométrica calculada en {t_fin - t_inicio:.2f} segundos.")
    
    # 2. DEFINICIÓN DE LOS POSTULADOS DE TU ARTÍCULO
    A_int = sp.Function('A_int')(r)
    B_int = sp.Function('B_int')(r)
    T_00 = sp.Function('T_00')(r)
    
    A_sch = 1 - (2 * G * M) / (c**2 * r)
    B_sch = 1 / A_sch
    
    sigma = 1 / (1 + (T_00 / K_0))
    
    # Ecuaciones constitutivas unificadas de tu modelo (Ec. 3 y 4)
    A_teoria = sigma * A_sch + (1 - sigma) * A_int
    B_teoria = sigma * B_sch + (1 - sigma) * B_int
    
    # 3. CÁLCULO DE DERIVADAS ANALÍTICAS EXACTAS MEDIANTE REGLA DE LA CADENA
    # Calculamos rigurosamente las derivadas primeras y segundas de tu modelo
    A_p1 = sp.diff(A_teoria, r)
    A_p2 = sp.diff(A_p1, r)
    B_p1 = sp.diff(B_teoria, r)
    
    # 4. INYECCIÓN COVARIANTE DE TU TEORÍA EN EL TENSOR DE EINSTEIN
    print("\n[PROCESO] Acoplando las funciones constitutivas elásticas de la teoría...")
    
    # Sustituimos A(r), B(r) y sus derivadas por las expresiones exactas de tu modelo
    sustituciones = [
        (sp.diff(A, (r, 2)), A_p2),
        (sp.diff(A, r), A_p1),
        (sp.diff(B, r), B_p1),
        (A, A_teoria),
        (B, B_teoria)
    ]
    
    G_00_final = G_tensor[0,0].subs(sustituciones)
    G_11_final = G_tensor[1,1].subs(sustituciones)
    R_final = R_escalar.subs(sustituciones)
    
    print("=======================================================================")
    print(" COMPONENTES ANALÍTICAS EXACTAS DE TU MODELO:")
    print("=======================================================================")
    
    print("\n> ESCALAR DE RICCI (R):")
    sp.pprint(R_final)
    print("\n" + "-"*70)
    
    print("\n> COMPONENTE TEMPORAL G[0,0] (Densidad de Energía Efectiva):")
    sp.pprint(G_00_final)
    print("\n" + "-"*70)
    
    print("\n> COMPONENTE RADIAL G[1,1] (Tensión Elástica del Vacío):")
    sp.pprint(G_11_final)
    print("=======================================================================")

if __name__ == "__main__":
    calcular_tensores_rigor_analitico()
