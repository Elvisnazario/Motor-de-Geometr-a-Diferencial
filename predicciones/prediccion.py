"""
prediccion.py
=============

Infraestructura para generar predicciones físicas
a partir del Motor de Geometría Diferencial (MGD).

Una predicción puede ser cualquier magnitud obtenida
a partir de objetos geométricos del motor.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp


class Prediccion:
    """
    Representa una predicción física.

    Ejemplos
    --------

    • radio mínimo

    • constante elástica

    • desplazamiento Lamb

    • sección eficaz

    • potencial efectivo

    • densidad de energía

    • etc.
    """

    def __init__(
        self,
        nombre,
        expresion,
        descripcion="",
        unidades="",
    ):

        self.nombre = nombre

        self.expresion = sp.simplify(
            expresion
        )

        self.descripcion = descripcion

        self.unidades = unidades

    # -------------------------------------------------

    def evaluar(
        self,
        sustituciones=None,
    ):

        if sustituciones is None:
            sustituciones = {}

        return sp.N(
            self.expresion.subs(
                sustituciones
            )
        )

    # -------------------------------------------------

    def simbolica(self):

        return self.expresion

    # -------------------------------------------------

    def latex(self):

        return sp.latex(
            self.expresion
        )

    # -------------------------------------------------

    def __repr__(self):

        return (
            f"Prediccion("
            f"{self.nombre})"
        )

    __str__ = __repr__