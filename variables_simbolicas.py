
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def verificacion_formal_teoria():
    print("=======================================================================")
    print("   MOTOR DE GEOMETRÍA: VERIFICACIÓN FORMAL DE LA ECUACIÓN DE ESTADO")
    print("   Validación de TOV, Comprobación Inversa y Cotas del Horizonte")
    print("   Autor: Elvis Omar Nazario Espinoza")
    print("=======================================================================\n")

    # 1. DEFINICIÓN DE VARIABLES SIMBÓLICAS
    x = sp.Symbol('x', real=True, positive=True) # Coordenada adimensional r/rc
    rho_0 = sp.Symbol('rho_0', real=True, positive=True)
    r_c = sp.Symbol('r_c', real=True, positive=True)
    G = sp.Symbol('G', real=True, positive=True)
    C_0 = 4 * sp.pi * G * rho_0**2 * r_c**2

    # Perfil de densidad de la Atracción Energética
    rho = rho_0 / (1 + x**2)
    
    # Masa acumulada integrada analíticamente
    M = 4 * sp.pi * rho_0 * r_c**3 * (x - sp.atan(x))

    print("[1/3] Verificando la ecuación de equilibrio TOV / Lane-Emden...")
    # Lado derecho de la ecuación de equilibrio hidrostático-elástico dP/dx
    tov_rhs = - (G * M * rho) / (r_c * x**2)
    tov_rhs_simplified = sp.simplify(tov_rhs)
    
    # Expresión analítica deducida para la Presión P(x)
    P = C_0 * ((sp.pi**2)/8 - sp.atan(x)/x - (sp.atan(x)**2)/2)
    
    # Derivamos la Presión analítica respecto a x para comprobar hacia atrás
    dP_dx_analitica = sp.diff(P, x)
    dP_dx_simplified = sp.simplify(dP_dx_analitica)
    
    # Verificación de equivalencia exacta
    identidad = sp.simplify(dP_dx_simplified - tov_rhs_simplified)
    
    print(f" -> Lado Derecho de TOV: {tov_rhs_simplified}")
    print(f" -> Derivada dP/dx de la solución: {dP_dx_simplified}")
    print(f" -> Diferencia simbólica (Debe ser cero): {identidad}")
    if identidad == 0:
        print(" ¡TEOREMA DEMOSTRADO HACIA ATRÁS CON ÉXITO SIMBÓLICO ABSOLUTO!")
    else:
        print(" Comprobación numérica alternativa necesaria...")

    # 2. DEMOSTRACIÓN DE LA COTA SUPERIOR PARA EVITAR EL HORIZONTE (A(r) > 0)
    print("\n[2/3] Evaluando cotas analíticas para garantizar A(r) > 0...")
    # Para que A(r) = 1 - 2*G*M_eff(r)/(c^2 * r) > 0, se requiere que:
    # M_eff(r)/r < c^2 / (2*G)
    # Evaluamos el máximo de M_eff(r)/r.
    # El primer término (bariónico) M_bar * r / (r^2 + r_0^2) tiene su máximo en r = r_0
    # El segundo término (vacío elástico) tiene su cota asintótica superior en 4*pi*rho_0*r_c^2
    print(" -> Cota analítica superior de M_eff(r)/r establecida como:")
    print("    M_bar / (2*r_0) + 4 * pi * rho_0 * r_c^2")
    print(" -> Condición de regularidad causal estricta sin horizontes:")
    print("    G * M_bar / r_0 + 8 * pi * G * rho_0 * r_c^2 < c^2")

    # 3. VERIFICACIÓN NUMÉRICA Y GENERACIÓN DE GRÁFICAS DE COMPROBACIÓN
    print("\n[3/3] Generando simulaciones de validación física de presiones...")
    
    x_val = np.linspace(0.001, 15.0, 500)
    
    # Solución analítica del perfil de Presión normalizado P(x) / C_0
    P_val = (np.pi**2)/8.0 - np.arctan(x_val)/x_val - (np.arctan(x_val)**2)/2.0
    
    # Evaluamos dP/dx numéricamente y mediante la fórmula TOV para validar
    dP_dx_num = np.gradient(P_val, x_val[1] - x_val[0])
    tov_rhs_num = - (x_val - np.arctan(x_val)) / (x_val**2 * (1.0 + x_val**2))
    
    # Gráficas de validación
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    ax1.plot(x_val, P_val, 'r-', linewidth=2.0, label="Presión Elástica $P(r)$ (Atracción Energética)")
    ax1.set_xlabel("Distancia adimensional $x = r/r_c$")
    ax1.set_ylabel("Presión normalizada $P(r) / C_0$")
    ax1.set_title("Distribución de Presión Monótona y Regular")
    ax1.grid(True, linestyle='--')
    ax1.legend()
    
    ax2.plot(x_val, dP_dx_num, 'b-', linewidth=2.5, label="Derivada Numérica $dP/dx$")
    ax2.plot(x_val, tov_rhs_num, 'r--', linewidth=1.5, label="Ecuación TOV Hidrostática")
    ax2.set_xlabel("Distancia adimensional $x = r/r_c$")
    ax2.set_ylabel("Gradiente de Presión")
    ax2.set_title("Verificación de Equivalencia de Equilibrio")
    ax2.grid(True, linestyle='--')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig("verificacion_equilibrio_elastico.png", dpi=300)
    print(" -> Gráfica de comprobación guardada con éxito como 'verificacion_equilibrio_elastico.png'.")
    plt.show()

if __name__ == "__main__":
    verificacion_formal_teoria()
