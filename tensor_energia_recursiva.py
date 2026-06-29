"""
tensor_energia_recursiva.py
===========================

Tensor de Energía Recursiva (TER).

Este tensor representa la respuesta constitutiva del
medio geométrico del Motor de Geometría Diferencial (MGD).

A diferencia del tensor energía-momento clásico de la
Relatividad General, el TER describe la energía almacenada
por la deformación del propio medio geométrico.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from nucleo.tensor import Tensor


class TensorEnergiaRecursiva(Tensor):
    """
    Tensor constitutivo del medio geométrico.

    T^μ_ν =
        diag(
            -ρ,
             p_r,
             p_t,
             p_t
        )

    donde

        ρ   = densidad efectiva de energía.

        p_r = presión radial.

        p_t = presión tangencial.

    En implementaciones posteriores estas cantidades
    podrán depender de:

        ψ
        ω
        K0
        σ
        o cualquier otra variable del modelo MGD.
    """

    def __init__(
        self,
        variedad,
        rho,
        p_r,
        p_t,
    ):

        dim = variedad.dimension

        if dim != 4:
            raise ValueError(
                "El Tensor de Energía Recursiva "
                "está definido actualmente "
                "para espacio-tiempo 4D."
            )

        componentes = {

            (0, 0): -sp.simplify(rho),

            (1, 1): sp.simplify(p_r),

            (2, 2): sp.simplify(p_t),

            (3, 3): sp.simplify(p_t),

        }

        super().__init__(

            nombre="TensorEnergiaRecursiva",

            variedad=variedad,

            componentes=componentes,

            indices=[
                ("μ", "abajo"),
                ("ν", "abajo"),
            ],
        )

        self.rho = sp.simplify(rho)
        self.p_r = sp.simplify(p_r)
        self.p_t = sp.simplify(p_t)

    # --------------------------------------------------
    # Magnitudes físicas
    # --------------------------------------------------

    @property
    def ecuacion_estado_radial(self):

        if self.rho == 0:
            return sp.nan

        return sp.simplify(
            self.p_r / self.rho
        )

    @property
    def ecuacion_estado_tangencial(self):

        if self.rho == 0:
            return sp.nan

        return sp.simplify(
            self.p_t / self.rho
        )

    @property
    def es_isotropo(self):

        return sp.simplify(
            self.p_r - self.p_t
        ) == 0

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (
            "TensorEnergiaRecursiva("
            f"rho={self.rho}, "
            f"pr={self.p_r}, "
            f"pt={self.p_t})"
        )