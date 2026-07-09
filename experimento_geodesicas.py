import sympy as sp
import time

def calcular_geodesicas_atraccion_energetica():
    print("=======================================================================")
    print(" CÁLCULO RIGUROSO DE GEODÉSICAS - TEORÍA DE LA ATRACCIÓN ENERGÉTICA")
    print(" Autor: Elvis Omar Nazario Espinoza")
    print("=======================================================================\n")
    
    # 1. DEFINICIÓN DE COORDENADAS Y VARIABLES
    t, r, theta, phi = sp.symbols('t r theta phi')
    tau = sp.symbols('tau') # Tiempo propio para la parametrización de la línea de universo
    
    G, c, M = sp.symbols('G c M', positive=True)
    K_0 = sp.symbols('K_0', positive=True)
    
    # Coordenadas como funciones del tiempo propio tau para las geodésicas
    t_tau = sp.Function('t')(tau)
    r_tau = sp.Function('r')(tau)
    theta_tau = sp.Function('theta')(tau)
    phi_tau = sp.Function('phi')(tau)
    
    # Velocidades de la partícula (derivadas respecto a tau)
    dt = sp.diff(t_tau, tau)
    dr = sp.diff(r_tau, tau)
    dtheta = sp.diff(theta_tau, tau)
    dphi = sp.diff(phi_tau, tau)
    
    # Aceleraciones de la partícula (segundas derivadas respecto a tau)
    d2t = sp.diff(dt, tau)
    d2r = sp.diff(dr, tau)
    d2theta = sp.diff(dtheta, tau)
    d2phi = sp.diff(dphi, tau)
    
    # 2. DEFINICIÓN DE LOS POSTULADOS DE TU ARTÍCULO (Inyección exacta)
    A_int = sp.Function('A_int')(r)
    B_int = sp.Function('B_int')(r)
    T_00 = sp.Function('T_00')(r)
    
    A_sch = 1 - (2 * G * M) / (c**2 * r)
    B_sch = 1 / A_sch
    
    sigma = 1 / (1 + (T_00 / K_0))
    A_teoria = sigma * A_sch + (1 - sigma) * A_int
    B_teoria = sigma * B_sch + (1 - sigma) * B_int
    
    # Derivadas de la métrica respecto a la coordenada radial 'r'
    A_p1 = sp.diff(A_teoria, r)
    B_p1 = sp.diff(B_teoria, r)
    
    # 3. CONSTRUCCIÓN DIRECTA DE LAS ECUACIONES GEODÉSICAS
    # Usamos la formulación Lagrangiana de la métrica estática esférica, la cual es matemáticamente
    # equivalente a calcular los Símbolos de Christoffel (Gamma), pero instantánea en CPU.
    print("[PROCESO] Derivando analíticamente los coeficientes de movimiento...")
    t_inicio = time.time()
    
    # Ecuaciones diferenciales crudas para simetría esférica general:
    # Ecuación Temporal (t)
    eq_t_cruda = d2t + (A_p1 / A_teoria) * dt * dr
    
    # Ecuación Radial (r)
    eq_r_cruda = d2r + (A_p1 / (2 * B_teoria)) * dt**2 + (B_p1 / (2 * B_teoria)) * dr**2 - (r / B_teoria) * dtheta**2 - (r * sp.sin(theta)**2 / B_teoria) * dphi**2
    
    # Ecuación Angular Theta
    eq_theta_cruda = d2theta + (2 / r) * dr * dtheta - sp.sin(theta) * sp.cos(theta) * dphi**2
    
    # Ecuación Angular Phi
    eq_phi_cruda = d2phi + (2 / r) * dr * dphi + 2 * (sp.cos(theta) / sp.sin(theta)) * dtheta * dphi
    
    # 4. ACOPLAMIENTO DE LA COORDENADA RADIAL DINÁMICA r(tau)
    # Reemplazamos la variable estática 'r' por la función dinámica 'r(tau)' para que las 
    # derivadas tengan sentido físico a lo largo de la trayectoria.
    print("[PROCESO] Mapeando la trayectoria sobre el radio dinámico r(tau)...")
    sustituciones_dinamicas = {
        r: r_tau,
        theta: theta_tau
    }
    
    eq_t_final = sp.together(sp.cancel(eq_t_cruda.subs(sustituciones_dinamicas)))
    eq_r_final = sp.together(sp.cancel(eq_r_cruda.subs(sustituciones_dinamicas)))
    eq_theta_final = sp.together(sp.cancel(eq_theta_cruda.subs(sustituciones_dinamicas)))
    eq_phi_final = sp.together(sp.cancel(eq_phi_cruda.subs(sustituciones_dinamicas)))
    
    t_fin = time.time()
    print(f"✔ [ÉXITO] Ecuaciones de las Geodésicas calculadas en {t_fin - t_inicio:.2f} segundos.")
    
    # 5. DESPLIEGUE ELEGANTE EN PANTALLA
    print("=======================================================================")
    print(" ECUACIONES DIFERENCIALES DE LAS GEODÉSICAS EXACTAS: d²x^mu / d_tau² = ...")
    print("=======================================================================")
    
    print("\n> [DIMENSIÓN 0] COMPONENTE TEMPORAL (Conservación de la Energía Orbital):")
    print("d²t/dtau² = ")
    sp.pprint(sp.solve(eq_t_final, d2t)[0])
    print("\n" + "-"*70)
    
    print("\n> [DIMENSIÓN 1] COMPONENTE RADIAL (Dinámica de Atracción Elástica Efectiva):")
    print("d²r/dtau² = ")
    sp.pprient_expr = sp.solve(eq_r_final, d2r)[0]
    sp.pprint(sp.together(sp.cancel(sp.solve(eq_r_final, d2r)[0])))
    print("\n" + "-"*70)
    
    print("\n> [DIMENSIÓN 2 Y 3] GEODÉSICAS ANGULARES (Sectores de Momento Angular):")
    print("d²theta/dtau² = ")
    sp.pprint(sp.solve(eq_theta_final, d2theta)[0])
    print("d²phi/dtau² = ")
    sp.pprint(sp.solve(eq_phi_final, d2phi)[0])
    print("=======================================================================")

if __name__ == "__main__":
    calcular_geodesicas_atraccion_energetica()
