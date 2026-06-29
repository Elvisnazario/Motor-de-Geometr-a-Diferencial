"""
tensor_energia_momento.py
=========================

Implementación del Tensor Energía-Momento.

Este objeto representa el tensor T_{μν} que actúa como
fuente de las ecuaciones de campo.

En el Motor de Geometría Diferencial (MGD) esta clase
servirá como clase base para las distintas fuentes
(material, electromagnética, medio elástico, etc.).

Autor:
    Elvis Omar Nazario Espinoza

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

import sympy as sp

from nucleo.tensor import Tensor


class TensorEnergiaMomento(Tensor):
    """
    Tensor Energía-Momento.

        T_{μν}

    Clase base para todas las fuentes físicas.
    """

    def __init__(
        self,
        variedad,
        componentes=None,
        nombre="TensorEnergiaMomento",
    ):

        if componentes is None:
            componentes = {}

        super().__init__(
            nombre=nombre,
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

    def traza(self, metrica):

        g_inv = metrica.inversa()

        tr = sp.Integer(0)

        for mu in range(self.dimension):
            for nu in range(self.dimension):

                tr += (
                    g_inv[mu, nu]
                    * self[(mu, nu)]
                )

        return sp.simplify(tr)

    # --------------------------------------------------
    # Conservación
    # --------------------------------------------------

    def es_nulo(self):

        return len(self.componentes) == 0

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (
            "TensorEnergiaMomento("
            f"dimension={self.dimension}, "
            f"componentes={len(self.componentes)})"
        )

    def __str__(self):

        return self.__repr__()