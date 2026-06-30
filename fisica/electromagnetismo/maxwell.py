"""
maxwell.py
==========

Implementación de las ecuaciones de Maxwell
en forma covariante.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from fisica.electromagnetismo.tensor_faraday import (
    TensorFaraday,
)


class Maxwell:
    """
    Ecuaciones de Maxwell.

    Implementa las dos familias:

        ∇μFμν = Jν

    y

        ∇[αFβγ] = 0
    """

    def __init__(
        self,
        faraday,
        conexion=None,
    ):

        if not isinstance(
            faraday,
            TensorFaraday,
        ):
            raise TypeError(
                "Se esperaba un TensorFaraday."
            )

        self.F = faraday
        self.variedad = faraday.variedad
        self.conexion = conexion

    # --------------------------------------------------
    # Corriente
    # --------------------------------------------------

    def corriente(self):

        """
        Construye simbólicamente

            Jμ

        como funciones arbitrarias.
        """

        coords = self.variedad.coordenadas

        return [

            sp.Function(f"J{i}")(*coords)

            for i in range(self.variedad.dimension)

        ]

    # --------------------------------------------------
    # Maxwell inhomogéneo
    # --------------------------------------------------

    def ecuacion_inhomogenea(self):

        """
        Devuelve formalmente

            ∇μFμν

        La derivada covariante completa se
        incorporará cuando el operador tensorial
        esté completamente estabilizado.
        """

        raise NotImplementedError(
            "Pendiente de conectar con "
            "derivada_covariante tensorial."
        )

    # --------------------------------------------------
    # Maxwell homogéneo
    # --------------------------------------------------

    def ecuacion_homogenea(self):

        """
        Calcula

            ∂αFβγ
            +
            ∂βFγα
            +
            ∂γFαβ

        que debe anularse.
        """

        coords = self.variedad.coordenadas
        dim = self.variedad.dimension

        resultado = {}

        for a in range(dim):

            for b in range(dim):

                for c in range(dim):

                    valor = sp.simplify(

                        sp.diff(
                            self.F[(b, c)],
                            coords[a],
                        )

                        +

                        sp.diff(
                            self.F[(c, a)],
                            coords[b],
                        )

                        +

                        sp.diff(
                            self.F[(a, b)],
                            coords[c],
                        )

                    )

                    if valor != 0:

                        resultado[(a, b, c)] = valor

        return resultado

    # --------------------------------------------------

    def satisface_bianchi(self):

        """
        Verifica la identidad de Bianchi
        del tensor de Faraday.
        """

        return len(
            self.ecuacion_homogenea()
        ) == 0

    # --------------------------------------------------

    def __repr__(self):

        return (
            "Maxwell("
            f"dimension={self.variedad.dimension}"
            ")"
        )

    __str__ = __repr__