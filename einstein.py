"""
einstein.py
===========

Construcción del Tensor de Einstein.

    G_{μν} = R_{μν} - 1/2 R g_{μν}

Autor:
    Elvis Omar Nazario Espinoza

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

import sympy as sp

from nucleo.tensor import Tensor
from ricci import Ricci
from escalar_ricci import EscalarRicci
from metrica import Metrica


class Einstein(Tensor):
    """
    Tensor de Einstein.

    G_{μν} = R_{μν} - (1/2) R g_{μν}
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

        escalar = EscalarRicci(
            ricci,
            metrica
        ).valor

        componentes = {}

        dim = self.variedad.dimension

        for mu in range(dim):
            for nu in range(dim):

                valor = sp.simplify(

                    ricci[(mu, nu)]

                    - sp.Rational(1, 2)

                    * escalar

                    * metrica.g_cov[mu, nu]

                )

                if valor != 0:
                    componentes[(mu, nu)] = valor

        super().__init__(
            nombre="Einstein",
            variedad=self.variedad,
            componentes=componentes,
            indices=[
                ("μ", "abajo"),
                ("ν", "abajo"),
            ],
        )

    # --------------------------------------------------
    # Utilidad
    # --------------------------------------------------

    def es_nulo(self):
        """
        Devuelve True si todas las componentes
        del tensor son cero.
        """

        return len(self.componentes) == 0

    # --------------------------------------------------
    # Acceso
    # --------------------------------------------------

    def __getitem__(self, clave):

        return super().__getitem__(clave)

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (
            f"Einstein("
            f"dimension={self.dimension}, "
            f"componentes={len(self.componentes)})"
        )

    def __str__(self):

        return self.__repr__()