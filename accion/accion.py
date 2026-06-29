"""
accion.py
=========

Representación de la Acción de una teoría de campos.

La acción constituye el objeto variacional fundamental del
Motor de Geometría Diferencial (MGD).

        S = ∫ L √|g| dⁿx

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from lagrangianos.lagrangiano_elastico import LagrangianoElastico


class Accion:
    """
    Acción variacional.

    Permite construir la acción simbólica de cualquier
    teoría implementada sobre el MGD.
    """

    def __init__(
        self,
        variedad,
        metrica,
        lagrangiano,
    ):

        if not isinstance(
            lagrangiano,
            LagrangianoElastico,
        ):
            raise TypeError(
                "Se esperaba un LagrangianoElastico."
            )

        self.variedad = variedad
        self.metrica = metrica
        self.lagrangiano = lagrangiano

    # --------------------------------------------------

    @property
    def determinante(self):
        """
        Determinante de la métrica.
        """

        return sp.simplify(
            self.metrica.g_cov.det()
        )

    # --------------------------------------------------

    @property
    def volumen(self):
        """
        Elemento de volumen:

            √|g|
        """

        return sp.sqrt(
            sp.Abs(
                self.determinante
            )
        )

    # --------------------------------------------------

    @property
    def densidad(self):
        """
        Densidad lagrangiana.

            √|g| L
        """

        return sp.simplify(

            self.volumen

            *

            self.lagrangiano.expresion

        )

    # --------------------------------------------------

    def integrando(self):
        """
        Devuelve el integrando completo.
        """

        return self.densidad

    # --------------------------------------------------

    def accion_formal(self):
        """
        Devuelve la acción escrita formalmente
        como integral múltiple.
        """

        expr = self.densidad

        for coord in reversed(
            self.variedad.coordenadas
        ):

            expr = sp.Integral(
                expr,
                coord,
            )

        return expr

    # --------------------------------------------------

    def latex(self):

        return sp.latex(
            self.accion_formal()
        )

    # --------------------------------------------------

    def __repr__(self):

        return (
            "Accion("
            f"dimension={self.variedad.dimension})"
        )

    __str__ = __repr__