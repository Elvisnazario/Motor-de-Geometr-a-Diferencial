import sympy as sp
import time
from motor_tensorial import crear_motor_tensorial, matriz_es_cero

def calcular_y_exportar_reporte_latex():
    print("=======================================================================")
    print(" EXTRACCIÓN RIGUROSA Y EXPORTACIÓN A LATEX - ATRACCIÓN ENERGÉTICA")
    print(" Autor: Elvis Omar Nazario Espinoza")
    print("=======================================================================\n")
    
    t, r, theta, phi = sp.symbols('t r theta phi')
    G, c, M = sp.symbols('G c M', positive=True)
    K_0 = sp.symbols('K_0', positive=True)
    
    # 1. GEOMETRÍA ABSTRACTA BASE
    A = sp.Function('A')(r)
    B = sp.Function('B')(r)
    metric_components = [-A, B, r**2, r**2 * sp.sin(theta)**2]
    
    print("[PROCESO] Ejecutando cálculo tensorial subyacente...")
    G_tensor, R_escalar = crear_motor_tensorial(metric_components, nombre_cache=None)
    
    # 2. POSTULADOS DE TU ARTÍCULO
    A_int = sp.Function('A_int')(r)
    B_int = sp.Function('B_int')(r)
    T_00 = sp.Function('T_00')(r)
    
    A_sch = 1 - (2 * G * M) / (c**2 * r)
    B_sch = 1 / A_sch
    
    sigma = 1 / (1 + (T_00 / K_0))
    A_teoria = sigma * A_sch + (1 - sigma) * A_int
    B_teoria = sigma * B_sch + (1 - sigma) * B_int
    
    # Derivadas analíticas exactas por regla de la cadena
    A_p1 = sp.diff(A_teoria, r)
    A_p2 = sp.diff(A_p1, r)
    B_p1 = sp.diff(B_teoria, r)
    
    # Sustitución de covariantes
    sustituciones = [
        (sp.diff(A, (r, 2)), A_p2),
        (sp.diff(A, r), A_p1),
        (sp.diff(B, r), B_p1),
        (A, A_teoria),
        (B, B_teoria)
    ]
    
    print("[PROCESO] Acoplando funciones constitutivas elásticas...")
    G_00_final = G_tensor[0,0].subs(sustituciones)
    G_11_final = G_tensor[1,1].subs(sustituciones)
    R_final = R_escalar.subs(sustituciones)
    
    # 3. GENERACIÓN AUTOMÁTICA DEL REPORTE LATEX (.tex)
    nombre_archivo = "reporte_tensores.tex"
    print(f"\n[PROCESO] Exportando ecuaciones analíticas completas a '{nombre_archivo}'...")
    
    with open(nombre_archivo, "w", encoding="utf-8") as f:
        f.write("% REPORTE GENERADO AUTOMÁTICAMENTE POR EL MOTOR DE GEOMETRÍA DIFERENCIAL\n")
        f.write("\\section{Resultados Analíticos de los Tensores de Campo}\n\n")
        f.write("A continuación se presentan las expresiones algebraicas exactas obtenidas ")
        f.write("mediante el cálculo tensorial computacional para el modelo de Atracción Energética.\n\n")
        
        f.write("\\subsection{Escalar de Ricci ($R$)}\n")
        f.write("\\begin{equation}\n")
        f.write(sp.latex(R_final) + "\n")
        f.write("\\end{equation}\n\n")
        
        f.write("\\subsection{Componente Temporal del Tensor de Einstein ($G_{00}$)}\n")
        f.write("\\begin{style}{\\small}\n") # Añadimos entorno para prevenir desbordes de página en LaTeX
        f.write("\\begin{equation}\n")
        f.write(sp.latex(G_00_final) + "\n")
        f.write("\\end{equation}\n")
        f.write("\\end{style}\n\n")
        
        f.write("\\subsection{Componente Radial del Tensor de Einstein ($G_{11}$)}\n")
        f.write("\\begin{equation}\n")
        f.write(sp.latex(G_11_final) + "\n")
        f.write("\\end{equation}\n")

    print("=======================================================================")
    print(f"✔ [ÉXITO] Reporte científico guardado con éxito en: {nombre_archivo}")
    print("=======================================================================")

if __name__ == "__main__":
    calcular_y_exportar_reporte_latex()
