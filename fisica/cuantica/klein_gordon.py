"""
klein_gordon.py
===============

Implementación de la ecuación de Klein-Gordon
sobre una variedad diferenciable.

Esta constituye la puerta de entrada del MGD
al régimen cuántico relativista.

Ecuación:

    □φ + m²φ = 0

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sys
from pathlib import Path
import sympy as sp

# Navega desde la ubicación actual hasta la carpeta donde está 'experimentos'
DIRECTORIO_ACTUAL = Path(__file__).resolve().parent
RAIZ_FISICA = DIRECTORIO_ACTUAL.parent

if str(RAIZ_FISICA) not in sys.path:
    sys.path.insert(0, str(RAIZ_FISICA))

# Importa directamente desde el directorio 'experimentos' estando dentro de 'fisica'
from experimentos.campo_escalar import CampoEscalar


class KleinGordon(CampoEscalar):
    """
    Campo escalar relativista de Klein-Gordon.
    """

    def __init__(
        self,
        variedad,
        masa,
        simbolo="phi",
    ):

        super().__init__(
            nombre="Klein-Gordon",
            variedad=variedad,
            simbolo=simbolo,
        )

        self.masa = sp.sympify(masa)

    # --------------------------------------------------

    def operador_dalembert(self, metrica):

        return self.laplaciano(metrica)

    # --------------------------------------------------

    def ecuaciones(self, metrica):

        """
        Devuelve

            □φ + m²φ
        """

        return sp.simplify(

            self.operador_dalembert(metrica)

            +

            self.masa**2

            * self.phi

        )

    # --------------------------------------------------

    def lagrangiano(self):

        """
        Densidad lagrangiana formal.

        La construcción completa mediante la
        métrica será implementada posteriormente.
        """

        return sp.Symbol("L_KG")

    # --------------------------------------------------

    def tensor_energia(self):

        raise NotImplementedError(
            "Pendiente de implementar."
        )

    # --------------------------------------------------

    def __repr__(self):

        return (

            "KleinGordon("

            f"masa={self.masa}"

            ")"

        )

    __str__ = __repr__
