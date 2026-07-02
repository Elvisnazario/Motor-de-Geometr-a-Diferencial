"""
fisica/invariantes.py
=====================
Módulo numérico de alta precisión para el cálculo de invariantes de curvatura.
Implementa diferencias finitas híbridas de segundo orden homogéneo tanto en el
espacio libre (centradas) como en la frontera cuántica r_min (hacia adelante).

Proyecto: Motor de Geometría Diferencial (MGD)
Autor: Elvis Omar Nazario Espinoza
"""

import numpy as np
from fisica.atraccion_energetica import (
    f_schwarzschild,
    f_disforme,
    radio_schwarzschild,
    radio_minimo
)

def calcular_derivadas_num_schwarzschild(masa, r):
    """
    Calcula analíticamente f, f' y f'' para la solución clásica de Schwarzschild.
    Nota: Para r < rs esta carta de coordenadas cambia de signo, lo cual es el
    comportamiento divergente clásico analizado en el experimento.
    """
    rs = radio_schwarzschild(masa)
    f = 1.0 - (rs / r)
    f_primo = rs / (r**2)
    f_segundo = -2.0 * rs / (r**3)
    return f, f_primo, f_segundo

def calcular_derivadas_num_mgd(masa, r):
    """
    Calcula f y sus derivadas primera y segunda para el MGD de forma adaptativa.
    Aplica diferencias finitas hacia adelante de SEGUNDO ORDEN en el borde físico r_min,
    y diferencias centradas puras en el espacio libre para garantizar precisión uniforme.
    """
    f = f_disforme(masa, r)
    rmin = radio_minimo(masa)
    
    # dx adaptativo absoluto óptimo para evitar errores de cancelación numéricos
    dx = np.maximum(1e-35, 1e-4 * np.abs(r))
    
    # Inicialización de arreglos para derivadas vectorizadas
    f_primo = np.zeros_like(r)
    f_segundo = np.zeros_like(r)
    
    # Máscara booleana para detectar puntos críticamente pegados al borde r_min
    en_borde = (r <= rmin + dx)
    fuera_borde = ~en_borde
    
    # CASO A: Puntos en el borde físico (Diferencias finitas hacia adelante de SEGUNDO ORDEN)
    if np.any(en_borde):
        r_b = r[en_borde]
        dx_b = dx[en_borde]
        f_b = f[en_borde]
        f_mas1 = f_disforme(masa, r_b + dx_b)
        f_mas2 = f_disforme(masa, r_b + 2.0 * dx_b)
        f_mas3 = f_disforme(masa, r_b + 3.0 * dx_b) # Punto de control para segundo orden estable
        
        # Fórmulas progresivas de alta precisión (orden 2)
        f_primo[en_borde] = (-3.0 * f_b + 4.0 * f_mas1 - f_mas2) / (2.0 * dx_b)
        f_segundo[en_borde] = (2.0 * f_b - 5.0 * f_mas1 + 4.0 * f_mas2 - f_mas3) / (dx_b**2)
        
    # CASO B: Puntos en el espacio libre (Diferencias finitas centradas estándar de orden 2)
    if np.any(fuera_borde):
        r_f = r[fuera_borde]
        dx_f = dx[fuera_borde]
        f_f = f[fuera_borde]
        f_mas = f_disforme(masa, r_f + dx_f)
        f_menos = f_disforme(masa, r_f - dx_f)
        
        f_primo[fuera_borde] = (f_mas - f_menos) / (2.0 * dx_f)
        f_segundo[fuera_borde] = (f_mas - 2.0 * f_f + f_menos) / (dx_f**2)
        
    return f, f_primo, f_segundo

def kretschmann_esferico(f, f_primo, f_segundo, r):
    """
    Implementación del escalar de Kretschmann para métricas estáticas,
    diagonales y esféricamente simétricas generales:
    K = (f'')^2 + 4*(f')^2 / r^2 + 4*(1-f)^2 / r^4
    """
    termino1 = f_segundo**2
    termino2 = 4.0 * (f_primo**2) / (r**2)
    termino3 = 4.0 * ((1.0 - f)**2) / (r**4)
    return termino1 + termino2 + termino3

def kretschmann_schwarzschild(masa, r):
    """Calcula el escalar de Kretschmann para Schwarzschild usando NumPy."""
    f, f_primo, f_segundo = calcular_derivadas_num_schwarzschild(masa, r)
    return kretschmann_esferico(f, f_primo, f_segundo, r)

def kretschmann_mgd(masa, r):
    """Calcula el escalar de Kretschmann para el MGD usando NumPy."""
    f, f_primo, f_segundo = calcular_derivadas_num_mgd(masa, r)
    return kretschmann_esferico(f, f_primo, f_segundo, r)
