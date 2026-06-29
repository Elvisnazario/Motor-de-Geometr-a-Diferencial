"""
ecuacion_campo_mgd.py
=====================

Primera ecuación de campo del
Motor de Geometría Diferencial (MGD).

Esta clase une el Tensor de Einstein con el
Tensor de Energía Recursiva del medio geométrico.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from einstein import Einstein
from tensor_energia_recursiva import TensorEnergiaRecursiva
from nucleo.tensor import Tensor


class EcuacionCampoMGD:
    """
    Ecuación fundamental del MGD.

        G_{μν} = κ T^(ER)_{μν}

    donde

        G_{μν}   Tensor de Einstein

        T^(ER)   Tensor de Energía Recursiva

        κ        Constante de acoplamiento
    """

    def __init__(

        self,

        einstein,

        tensor_energia,

        kappa,

    ):

        if not isinstance(einstein, Einstein):
            raise TypeError(
                "Se esperaba un Tensor de Einstein."
            )

        if not isinstance(
            tensor_energia,
            TensorEnergiaRecursiva,
        ):
            raise TypeError(
                "Se esperaba un TensorEnergiaRecursiva."
            )

        self.einstein = einstein
        self.tensor_energia = tensor_energia

        self.kappa = sp.simplify(kappa)

        self.variedad = einstein.variedad

    # --------------------------------------------------
    # Residuo de la ecuación
    # --------------------------------------------------

    def residuo(self):

        componentes = {}

        dim = self.variedad.dimension

        for mu in range(dim):
            for nu in range(dim):

                valor = sp.simplify(

                    self.einstein[(mu, nu)]

                    -

                    self.kappa

                    *

                    self.tensor_energia[(mu, nu)]

                )

                if valor != 0:

                    componentes[(mu, nu)] = valor

        return Tensor(

            nombre="ResiduoMGD",

            variedad=self.variedad,

            componentes=componentes,

            indices=[
                ("μ", "abajo"),
                ("ν", "abajo"),
            ],
        )

    # --------------------------------------------------
    # Verificación
    # --------------------------------------------------

    def satisfecha(self):

        return len(
            self.residuo().componentes
        ) == 0

    # --------------------------------------------------
    # Norma del error
    # --------------------------------------------------

    def error_maximo(self):

        r = self.residuo()

        if len(r.componentes) == 0:
            return sp.Integer(0)

        return max(

            sp.Abs(v)

            for v in r.componentes.values()

        )

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (

            "EcuacionCampoMGD("

            f"kappa={self.kappa})"

        )