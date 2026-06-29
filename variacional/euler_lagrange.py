"""
euler_lagrange.py
=================

Operador de Euler-Lagrange del
Motor de Geometría Diferencial (MGD).

Calcula las ecuaciones de movimiento obtenidas
a partir de un lagrangiano.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp


class EulerLagrange:
    """
    Calculadora de ecuaciones de Euler-Lagrange.

    Implementa

        ∂L/∂φ
        -
        ∂_μ(
            ∂L/∂(∂_μφ)
        )
        = 0
    """

    def __init__(
        self,
        lagrangiano,
        coordenadas,
    ):

        self.L = lagrangiano
        self.coords = coordenadas

    # --------------------------------------------------

    def ecuacion(
        self,
        campo,
    ):
        """
        Calcula la ecuación de Euler-Lagrange
        para un campo escalar.
        """

        termino1 = sp.diff(
            self.L,
            campo,
        )

        termino2 = 0

        for coord in self.coords:

            derivada = sp.diff(
                campo,
                coord,
            )

            parcial = sp.diff(
                self.L,
                derivada,
            )

            termino2 += sp.diff(
                parcial,
                coord,
            )

        return sp.simplify(
            termino1
            -
            termino2
        )

    # --------------------------------------------------

    def verificar(
        self,
        campo,
    ):
        """
        Devuelve True si la ecuación se anula.
        """

        return (
            self.ecuacion(campo)
            == 0
        )

    # --------------------------------------------------

    def __repr__(self):

        return (
            "EulerLagrange("
            f"dim={len(self.coords)})"
        )

    __str__ = __repr__