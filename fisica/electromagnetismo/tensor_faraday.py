"""
tensor_faraday.py
=================

Construcción del tensor electromagnético de Faraday

    F_{μν} = ∂μAν − ∂νAμ

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from nucleo.tensor import Tensor
from fisica.electromagnetismo.potencial import (
    PotencialElectromagnetico,
)


class TensorFaraday(Tensor):
    """
    Tensor electromagnético de Faraday.

    Es un tensor covariante antisimétrico de rango 2.
    """

    def __init__(self, potencial):

        if not isinstance(
            potencial,
            PotencialElectromagnetico,
        ):
            raise TypeError(
                "Se esperaba un PotencialElectromagnetico."
            )

        self.potencial = potencial
        self.variedad = potencial.variedad

        coords = self.variedad.coordenadas
        dim = self.variedad.dimension

        componentes = {}

        for mu in range(dim):

            for nu in range(dim):

                valor = sp.simplify(

                    sp.diff(
                        potencial[nu],
                        coords[mu],
                    )

                    -

                    sp.diff(
                        potencial[mu],
                        coords[nu],
                    )

                )

                if valor != 0:

                    componentes[(mu, nu)] = valor

        super().__init__(

            nombre="Faraday",

            variedad=self.variedad,

            componentes=componentes,

            indices=[

                ("μ", "abajo"),

                ("ν", "abajo"),

            ],

        )

    # --------------------------------------------------

    def es_antisimetrico(self):
        """
        Verifica

            Fμν = −Fνμ
        """

        dim = self.dimension

        for mu in range(dim):

            for nu in range(dim):

                a = self[(mu, nu)]
                b = self[(nu, mu)]

                if sp.simplify(a + b) != 0:
                    return False

        return True

    # --------------------------------------------------

    def dual(self):
        """
        Reservado para la futura implementación
        del tensor dual

            *Fμν
        """

        raise NotImplementedError(
            "Tensor dual pendiente."
        )

    # --------------------------------------------------

    def __repr__(self):

        return (

            f"TensorFaraday("
            f"dimension={self.dimension})"

        )

    __str__ = __repr__