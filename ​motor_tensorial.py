import sympy as sp
import os
import pickle
import time

# Identificador único de versión matemática del motor geométrico
VERSION_MOTOR = "RG_MGD_V1.3_PRODUCCION"

def matriz_es_cero(M):
    """
    Evalúa elemento por elemento si una matriz de SymPy es matemáticamente cero.
    Evita que .is_zero_matrix devuelva 'None' ante expresiones algebraicas complejas.
    """
    return all(sp.simplify(sp.cancel(x)) == 0 for x in M)

def crear_motor_tensorial(metric_components, nombre_cache=None):
    """
    Motor tensorial puro de alto rendimiento.
    Usa sp.cancel() en pasos intermedios para evitar la explosión en RAM
    y reserva el esfuerzo de simplificación para el tensor final.
    """
    if nombre_cache and os.path.exists(f"{nombre_cache}.pkl"):
        try:
            print(f"-> [CACHÉ] Intentando cargar tensores precalculados desde '{nombre_cache}.pkl'...")
            with open(f"{nombre_cache}.pkl", "rb") as f:
                datos = pickle.load(f)
            
            if datos.get("version") == VERSION_MOTOR:
                print("✔ [CACHÉ] Datos válidos cargados de forma instantánea.")
                return datos["G_matriz"], datos["R_escalar"]
            else:
                print("⚠ [CACHÉ] Versión desactualizada detectada. Forzando recalculo analítico...")
        except Exception:
            print("⚠ [CACHÉ] Archivo corrupto o ilegible detectado. Forzando recalculo analítico...")

    t, r, theta, phi = sp.symbols('t r theta phi')
    variables = [t, r, theta, phi]
    
    # 1. MATRIZ MÉTRICA COVARIANTE E INVERSA
    g = sp.Matrix([
        [metric_components[0], 0, 0, 0],
        [0, metric_components[1], 0, 0],
        [0, 0, metric_components[2], 0],
        [0, 0, 0, metric_components[3]]
    ])
    g_inv = g.inv()
    
    # 2. SÍMBOLOS DE CHRISTOFFEL
    # Optimización: Se usa sp.cancel() para mantener las fracciones limpias sin el costo de simplify()
    christoffel = {}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                suma = 0
                for l in range(4):
                    dg_il_dj = sp.diff(g[i, l], variables[j])
                    dg_jl_di = sp.diff(g[j, l], variables[i])
                    dg_ij_dl = sp.diff(g[i, j], variables[l])
                    suma += g_inv[k, l] * (dg_il_dj + dg_jl_di - dg_ij_dl)
                christoffel[(k, i, j)] = sp.cancel(0.5 * suma)

    # 3. TENSOR DE RICCI (Convenio estándar de contracción indexada)
    ricci_matriz = sp.zeros(4, 4)
    for mu in range(4):
        for nu in range(4):
            valor = 0
            for lam in range(4):
                term1 = sp.diff(christoffel[(lam, mu, nu)], variables[lam])
                term2 = sp.diff(christoffel[(lam, mu, lam)], variables[nu])
                
                term3 = sum(christoffel[(lam, lam, sig)] * christoffel[(sig, mu, nu)] for sig in range(4))
                term4 = sum(christoffel[(lam, nu, sig)] * christoffel[(sig, mu, lam)] for sig in range(4))
                
                valor += (term1 - term2 + term3 - term4)
            ricci_matriz[mu, nu] = valor

    # 4. TENSOR DE EINSTEIN Y ESCALAR DE RICCI (Simplificación única al final)
    ricci_matriz = ricci_matriz.applyfunc(sp.cancel)
    
    R_escalar = sp.cancel(
        sum(g_inv[i, j] * ricci_matriz[i, j] for i in range(4) for j in range(4))
    )

    G_matriz = sp.zeros(4, 4)
    for mu in range(4):
        for nu in range(4):
            # Combinamos cancel y together para preparar el terreno analítico eficientemente
            G_matriz[mu, nu] = sp.together(sp.cancel(ricci_matriz[mu, nu] - 0.5 * g[mu, nu] * R_escalar))
            
    # Guardar en caché con el protocolo más alto de serialización binaria de Python
    if nombre_cache:
        try:
            print(f"-> [CACHÉ] Guardando resultados (Versión: {VERSION_MOTOR}) en '{nombre_cache}.pkl'...")
            with open(f"{nombre_cache}.pkl", "wb") as f:
                pickle.dump(
                    {"version": VERSION_MOTOR, "G_matriz": G_matriz, "R_escalar": R_escalar}, 
                    f, 
                    protocol=pickle.HIGHEST_PROTOCOL
                )
        except Exception as e:
            print(f"⚠ No se pudo escribir el caché en el disco: {e}")

    return G_matriz, R_escalar

def ejecutar_bateria_pruebas_avanzadas():
    r, M, Q, Lambda = sp.symbols('r M Q Lambda')
    theta = sp.symbols('theta')
    
    print("=======================================================================")
    print(" INICIANDO VALIDACIÓN DE LA ECUACIÓN DE CAMPO MULTIESCALA")
    print("=======================================================================\n")
    
    # -------------------------------------------------------------------------
    # TEST 1: MINKOWSKI
    # -------------------------------------------------------------------------
    print("[TEST 1] Evaluando Espaciotiempo de Minkowski (Plano)...")
    m_minkowski = [-1, 1, r**2, r**2 * sp.sin(theta)**2]
    G_mink, R_mink = crear_motor_tensorial(m_minkowski, nombre_cache="cache_minkowski")
    
    if matriz_es_cero(G_mink) and sp.simplify(sp.cancel(R_mink)) == 0:
        print("✔ PASÓ: Minkowski es plano (G = 0, R = 0).")
    else:
        print(f"❌ FALLÓ Minkowski.")
        return

    # -------------------------------------------------------------------------
    # TEST 2: SCHWARZSCHILD
    # -------------------------------------------------------------------------
    print("\n[TEST 2] Evaluando Métrica de Schwarzschild (Vacío)...")
    f_schwarz = 1 - 2*M/r
    m_schwarz = [-f_schwarz, 1/f_schwarz, r**2, r**2 * sp.sin(theta)**2]
    G_schw, _ = crear_motor_tensorial(m_schwarz, nombre_cache="cache_schwarzschild")
    
    if matriz_es_cero(G_schw):
        print("✔ PASÓ: Schwarzschild en el vacío es consistente (G = 0).")
    else:
        print(f"❌ FALLÓ Schwarzschild.")
        return

    # -------------------------------------------------------------------------
    # TEST 3: REISSNER-NORDSTRÖM (Electromagnetismo)
    # -------------------------------------------------------------------------
    print("\n[TEST 3] Evaluando Reissner-Nordström (Carga Q)...")
    t_i3 = time.time()
    f_reissner = 1 - 2*M/r + Q**2/r**2
    m_reissner = [-f_reissner, 1/f_reissner, r**2, r**2 * sp.sin(theta)**2]
    G_reis, R_reis = crear_motor_tensorial(m_reissner, nombre_cache="cache_reissner")
    t_f3 = time.time()
    
    if sp.simplify(sp.cancel(R_reis)) == 0:
        print(f"✔ PASÓ: Traza electromagnética nula verificada (R = 0). Tiempo: {t_f3 - t_i3:.2f}s")
    else:
        print(f"❌ FALLÓ Reissner-Nordström en contracción escalar.")
        return

    # -------------------------------------------------------------------------
    # TEST 4: DE SITTER (Vacío Cuántico / Lambda)
    # -------------------------------------------------------------------------
    print("\n[TEST 4] Evaluando de Sitter (Constante Cosmológica Lambda)...")
    f_desitter = 1 - (Lambda/3)*r**2
    m_desitter = [-f_desitter, 1/f_desitter, r**2, r**2 * sp.sin(theta)**2]
    G_ds, R_ds = crear_motor_tensorial(m_desitter, nombre_cache="cache_desitter")
    
    g_ds_formal = sp.Matrix([
        [-f_desitter, 0, 0, 0],
        [0, 1/f_desitter, 0, 0],
        [0, 0, r**2, 0],
        [0, 0, 0, r**2 * sp.sin(theta)**2]
    ])
    
    tensor_balance_total = G_ds + Lambda * g_ds_formal
    
    if matriz_es_cero(tensor_balance_total):
        print("✔ PASÓ: La ecuación cuántica de vacío cierra en CERO absoluto.")
    else:
        print("❌ FALLÓ de Sitter.")
        return

    print("\n=======================================================================")
    print(" MOTOR TOTALMENTE HOMOLOGADO Y CONFIGURADO AL 100%")
    print("=======================================================================")

if __name__ == "__main__":
    ejecutar_bateria_pruebas_avanzadas()
