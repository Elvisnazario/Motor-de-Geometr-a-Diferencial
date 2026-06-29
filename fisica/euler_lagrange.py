"""
euler_lagrange.py
=================

Operador general de Euler-Lagrange.

Permite obtener las ecuaciones de movimiento a partir de una
densidad lagrangiana.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp


class EulerLagrange:
    """
    Operador de Euler-Lagrange.

    Calcula

        ∂L/∂φ
        -
        ∂μ(∂L/∂(∂μφ))
        = 0
    """

    def __init__(
        self,
        lagrangiano,
        campo,
        coordenadas,
    ):

        self.L = lagrangiano
        self.campo = campo
        self.coords = coordenadas

    # --------------------------------------------------

    def ecuacion(self):

        phi = self.campo

        # Primer término
        termino1 = sp.diff(
            self.L,
            phi,
        )

        # Segundo término
        termino2 = 0

        for coord in self.coords:

            derivada_campo = sp.diff(
                phi,
                coord,
            )

            parcial = sp.diff(
                self.L,
                derivada_campo,
            )

            termino2 += sp.diff(
                parcial,
                coord,
            )

        return sp.simplify(
            termino1 - termino2
        )

    # --------------------------------------------------

    def __call__(self):

        return self.ecuacion()

    # --------------------------------------------------

    def __repr__(self):

        return (
            "EulerLagrange()"
        )

    __str__ = __repr__