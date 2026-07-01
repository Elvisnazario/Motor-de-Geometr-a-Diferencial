"""
fisica/atraccion_energetica.py
==============================

Implementación oficial de la teoría de Atracción Energética.

Toda la formulación física de la hipótesis debe vivir aquí.

El resto del MGD únicamente consumirá estas funciones.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import numpy as np

# 1. IMPORTACIÓN DEL ADMINISTRADOR EXCLUSIVO DE TU ARQUITECTURA
from constantes.parametros_fundamentales import ParametrosFundamentales

# 2. INICIALIZACIÓN Y REGISTRO DE CONSTANTES UNIVERSALES
# Usamos tu clase gestora para centralizar la física del universo
registro = ParametrosFundamentales()

# Constantes Fundamentales (Valores CODATA de alta precisión)
registro.registrar("G", valor=6.67430e-11, descripcion="Constante de gravitación universal")
registro.registrar("c", valor=299792458, descripcion="Velocidad de la luz en el vacío")
registro.registrar("hbar", valor=1.054571817e-34, descripcion="Constante de Planck reducida")

# Extracción de valores base mediante tus métodos nativos para el cálculo numérico
G = registro.valor("G")
c = registro.valor("c")
hbar = registro.valor("hbar")

# Parámetros derivados en la escala cuántica de Planck
l_planck_num = (hbar * G / c**3)**0.5
m_planck_num = (hbar * c / G)**0.5

# Registro oficial de los parámetros derivados dentro de tu administrador
registro.registrar("l_planck", valor=l_planck_num, descripcion="Longitud de Planck")
registro.registrar("m_planck", valor=m_planck_num, descripcion="Masa de Planck")

l_planck = registro.valor("l_planck")
m_planck = registro.valor("m_planck")


# ----------------------------------------------------------
# Radio de Schwarzschild
# ----------------------------------------------------------

def radio_schwarzschild(masa):
    """
    Radio de Schwarzschild clásico.
    """

    return (2 * G * masa) / (c ** 2)


# ----------------------------------------------------------
# Radio mínimo de la teoría
# ----------------------------------------------------------

def radio_minimo(masa):
    """
    Radio mínimo absoluto de la teoría.

        r_min = l_p (M/m_p)^(1/3)
    """

    return l_planck * (masa / m_planck) ** (1.0 / 3.0)


# ----------------------------------------------------------
# Energía crítica del vacío
# ----------------------------------------------------------

def energia_critica_vacio():
    """
    Densidad crítica del vacío.
    """

    return (
        m_planck * c ** 2
    ) / (
        (4.0 / 3.0)
        * np.pi
        * l_planck ** 3
    )


# ----------------------------------------------------------
# Densidad gravitatoria efectiva
# ----------------------------------------------------------

def densidad_energia(masa, r):
    """
    T00 efectiva del campo gravitatorio.
    """

    return (
        G * masa ** 2
    ) / (
        8.0
        * np.pi
        * r ** 4
    )


# ----------------------------------------------------------
# Función constitutiva
# ----------------------------------------------------------

def sigma(masa, r):
    """
    Respuesta constitutiva del vacío.

        σ(E)=1/(1+T00/K0)
    """

    T00 = densidad_energia(
        masa,
        r,
    )

    K0 = energia_critica_vacio()

    return 1.0 / (
        1.0
        + (T00 / K0)
    )


# ----------------------------------------------------------
# Función métrica clásica
# ----------------------------------------------------------

def f_schwarzschild(masa, r):
    """
    Función métrica clásica.
    """

    rs = radio_schwarzschild(masa)

    return 1.0 - (rs / r)


# ----------------------------------------------------------
# Función métrica disforme
# ----------------------------------------------------------

def f_disforme(masa, r):
    """
    Función métrica oficial del MGD.
    """

    rs = radio_schwarzschild(masa)

    return 1.0 - (
        (rs / r)
        * sigma(
            masa,
            r,
        )
    )
