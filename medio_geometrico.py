"""
medio_geometrico.py
==================

Modelo del Medio Geométrico del MGD.

Integra la Ley Constitutiva con el Tensor de Energía
Recursiva.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

from ley_constitutiva import LeyConstitutiva
from tensor_energia_recursiva import TensorEnergiaRecursiva


class MedioGeometrico:
    """
    Representa el estado físico del medio.

    A partir de los parámetros fundamentales del MGD
    construye automáticamente el Tensor de Energía
    Recursiva.
    """

    def __init__(
        self,
        variedad,
        psi,
        omega,
        sigma,
        K0,
    ):

        self.variedad = variedad

        self.ley = LeyConstitutiva(
            psi=psi,
            omega=omega,
            sigma=sigma,
            K0=K0,
        )

        self.tensor = TensorEnergiaRecursiva(
            variedad=variedad,
            rho=self.ley.rho,
            p_r=self.ley.p_r,
            p_t=self.ley.p_t,
        )

    @property
    def rho(self):
        return self.ley.rho

    @property
    def presion_radial(self):
        return self.ley.p_r

    @property
    def presion_tangencial(self):
        return self.ley.p_t

    @property
    def tensor_energia(self):
        return self.tensor

    @property
    def energia_total(self):
        return self.ley.energia_total

    def __repr__(self):
        return (
            "MedioGeometrico("
            f"rho={self.rho}, "
            f"pr={self.presion_radial}, "
            f"pt={self.presion_tangencial})"
        )