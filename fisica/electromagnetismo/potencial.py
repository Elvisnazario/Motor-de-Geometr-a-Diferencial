"""
potencial.py
============

Implementación del potencial electromagnético de cuatro componentes.

A partir de este objeto se construirá posteriormente el tensor
electromagnético de Faraday.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from fisica.campo import CampoFisico


class PotencialElectromagnetico(CampoFisico):
    """
    Cuatripotencial electromagnético

        A_μ

    formado por cuatro funciones escalares sobre la variedad.
    """

    def __init__(
        self,
        variedad,
        simbolo="A",
    ):

        super().__init__(
            nombre="PotencialElectromagnetico",
            variedad=variedad,
        )

        coords = variedad.coordenadas

        self.componentes = [

            sp.Function(f"{simbolo}{i}")(*coords)

            for i in range(variedad.dimension)

        ]

    # --------------------------------------------------

    def __getitem__(self, indice):

        return self.componentes[indice]

    # --------------------------------------------------

    def gradiente(self):

        """
        Devuelve

            ∂ν Aμ

        como matriz.
        """

        dim = self.variedad.dimension

        coords = self.variedad.coordenadas

        G = sp.MutableDenseMatrix.zeros(
            dim,
            dim,
        )

        for mu in range(dim):

            for nu in range(dim):

                G[mu, nu] = sp.diff(

                    self.componentes[mu],

                    coords[nu],

                )

        return G

    # --------------------------------------------------

    def lagrangiano(self):

        raise NotImplementedError(
            "Será implementado junto con Maxwell."
        )

    # --------------------------------------------------

    def tensor_energia(self):

        raise NotImplementedError

    # --------------------------------------------------

    def ecuaciones(self):

        raise NotImplementedError

    # --------------------------------------------------

    def __repr__(self):

        return (

            "PotencialElectromagnetico("

            f"dimension={self.variedad.dimension}"

            ")"

        )

    __str__ = __repr__