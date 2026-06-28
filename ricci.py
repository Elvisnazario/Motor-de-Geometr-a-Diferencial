"""
ricci.py
========

Construcción del Tensor de Ricci a partir del Tensor de Riemann.

R_{μν} = R^α_{μαν}

Autor:
    Elvis Omar Nazario Espinoza

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

import sympy as sp

from nucleo.tensor import Tensor
from riemann import Riemann


class Ricci(Tensor):
    """
    Tensor de Ricci obtenido por contracción del
    tensor de Riemann.

        R_{μν} = R^α_{μαν}
    """

    def __init__(self, riemann):

        if not isinstance(riemann, Riemann):
            raise TypeError(
                "Se esperaba un objeto Riemann."
            )

        self.riemann = riemann

        variedad = riemann.variedad
        dim = variedad.dimension

        componentes = {}

        #
        # Contracción:
        #
        # R_{μν} = Σ_α R^α_{μαν}
        #
        for mu in range(dim):
            for nu in range(dim):

                suma = 0

                for alpha in range(dim):

                    suma += riemann[
                        alpha,
                        mu,
                        alpha,
                        nu
                    ]

                suma = sp.simplify(suma)

                if suma != 0:
                    componentes[(mu, nu)] = suma

        super().__init__(
            nombre="Ricci",
            variedad=variedad,
            componentes=componentes,
            indices=[
                ("μ", "abajo"),
                ("ν", "abajo"),
            ],
        )

    # --------------------------------------------------
    # Traza
    # --------------------------------------------------

    def es_nulo(self):
        """
        Devuelve True si todas las componentes son cero.
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
            f"Ricci("
            f"dimension={self.dimension}, "
            f"componentes={len(self.componentes)})"
        )