"""
escalar_ricci.py
================

Construcción del Escalar de Ricci.

    R = g^{μν} R_{μν}

Autor:
    Elvis Omar Nazario Espinoza

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

import sympy as sp

from ricci import Ricci
from metrica import Metrica


class EscalarRicci:
    """
    Escalar de Ricci.

    Se obtiene mediante la contracción

        R = g^{μν} R_{μν}
    """

    def __init__(self, ricci, metrica):

        if not isinstance(ricci, Ricci):
            raise TypeError(
                "Se esperaba un objeto Ricci."
            )

        if not isinstance(metrica, Metrica):
            raise TypeError(
                "Se esperaba un objeto Metrica."
            )

        self.ricci = ricci
        self.metrica = metrica
        self.variedad = ricci.variedad

        self.valor = self._calcular()

    # --------------------------------------------------
    # Cálculo
    # --------------------------------------------------

    def _calcular(self):

        dim = self.variedad.dimension

        g_inv = self.metrica.inversa()

        escalar = 0

        for mu in range(dim):
            for nu in range(dim):

                escalar += (
                    g_inv[mu, nu]
                    * self.ricci[(mu, nu)]
                )

        return sp.simplify(escalar)

    # --------------------------------------------------
    # Conversión
    # --------------------------------------------------

    def __sympy__(self):
        return self.valor

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (
            f"EscalarRicci({self.valor})"
        )

    def __str__(self):

        return str(self.valor)