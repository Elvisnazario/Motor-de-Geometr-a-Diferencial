"""
tensor_energia_maxwell.py
=========================

Tensor energía-momento del campo electromagnético.

              1
Tμν = FμαFν α - ─ gμν FαβFαβ
              4

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from nucleo.tensor import Tensor
from fisica.electromagnetismo.tensor_faraday import (
    TensorFaraday,
)
from metrica import Metrica


class TensorEnergiaMaxwell(Tensor):
    """
    Tensor energía-momento del campo electromagnético.
    """

    def __init__(
        self,
        faraday,
        metrica,
    ):

        if not isinstance(
            faraday,
            TensorFaraday,
        ):
            raise TypeError(
                "Se esperaba un TensorFaraday."
            )

        if not isinstance(
            metrica,
            Metrica,
        ):
            raise TypeError(
                "Se esperaba una Metrica."
            )

        self.F = faraday
        self.g = metrica
        self.variedad = faraday.variedad

        dim = self.variedad.dimension

        ginv = metrica.g_con

        componentes = {}

        # --------------------------------------------------
        # Invariante FαβFαβ
        # --------------------------------------------------

        invariante = 0

        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):

                        invariante += (

                            ginv[a, c]

                            * ginv[b, d]

                            * faraday[(a, b)]

                            * faraday[(c, d)]

                        )

        invariante = sp.simplify(invariante)

        # --------------------------------------------------
        # Tensor energía
        # --------------------------------------------------

        for mu in range(dim):

            for nu in range(dim):

                termino1 = 0

                for alpha in range(dim):

                    for beta in range(dim):

                        termino1 += (

                            faraday[(mu, alpha)]

                            * ginv[alpha, beta]

                            * faraday[(nu, beta)]

                        )

                termino2 = (

                    sp.Rational(1, 4)

                    * metrica.g_cov[mu, nu]

                    * invariante

                )

                valor = sp.simplify(

                    termino1 - termino2

                )

                if valor != 0:

                    componentes[(mu, nu)] = valor

        super().__init__(

            nombre="TensorEnergiaMaxwell",

            variedad=self.variedad,

            componentes=componentes,

            indices=[

                ("μ", "abajo"),

                ("ν", "abajo"),

            ],

        )

    # --------------------------------------------------

    def traza(self):

        """
        Calcula

            Tμμ

        En vacío debe ser cero.
        """

        dim = self.dimension

        resultado = 0

        for mu in range(dim):

            for nu in range(dim):

                resultado += (

                    self.g.g_con[mu, nu]

                    * self[(mu, nu)]

                )

        return sp.simplify(resultado)

    # --------------------------------------------------

    def es_sin_traza(self):

        return self.traza() == 0

    # --------------------------------------------------

    def __repr__(self):

        return (

            "TensorEnergiaMaxwell("

            f"componentes={len(self.componentes)}"

            ")"

        )

    __str__ = __repr__