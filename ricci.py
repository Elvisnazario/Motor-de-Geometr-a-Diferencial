"""
ricci.py
=========

Construcción del tensor de Ricci a partir del tensor de Riemann.

Autor:
    Elvis Omar Nazario Espinoza

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

from tensor import Tensor
from riemann import Riemann


class Ricci(Tensor):
    """
    Tensor de Ricci.

    Se obtiene mediante la contracción

        R_{μν} = R^α_{μαν}

    del tensor de Riemann.
    """

    def __init__(self, variedad):
        super().__init__(
            nombre="Ricci",
            variedad=variedad,
            indices=[
                ("μ", "abajo"),
                ("ν", "abajo"),
            ],
        )

    @classmethod
    def desde_riemann(cls, riemann):
        """
        Construye el tensor de Ricci contrayendo
        el primer índice superior con el segundo
        índice inferior del tensor de Riemann.

        Parámetros
        ----------
        riemann : Riemann

        Retorna
        -------
        Ricci
        """

        if not isinstance(riemann, Riemann):
            raise TypeError(
                "Se esperaba un objeto de tipo Riemann."
            )

        ricci = cls(riemann.variedad)

        # Contracción:
        #
        # R_{μν} = R^α_{μαν}
        #
        # El método contraer() será el responsable
        # de realizar la suma de Einstein.
        return riemann.contraer(
            posicion_propia=0,
            posicion_otra=1
        )
