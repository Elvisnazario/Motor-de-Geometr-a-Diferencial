"""
ecuacion_campo.py
=================

Representación general de una ecuación tensorial de campo.

Permite comparar un miembro geométrico con un miembro físico
sin asumir una teoría particular.

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

import sympy as sp


class EcuacionCampo:
    """
    Representa una ecuación tensorial

        L_{μν} = R_{μν}

    donde ambos miembros pueden pertenecer a cualquier teoría.
    """

    def __init__(
        self,
        miembro_izquierdo,
        miembro_derecho,
        nombre="EcuacionCampo",
    ):

        self.nombre = nombre
        self.izquierdo = miembro_izquierdo
        self.derecho = miembro_derecho

    # --------------------------------------------------
    # Residuo
    # --------------------------------------------------

    def residuo(self):

        residuo = {}

        claves = (
            set(self.izquierdo.componentes.keys())
            |
            set(self.derecho.componentes.keys())
        )

        for clave in claves:

            valor = sp.simplify(

                self.izquierdo.componentes.get(clave, 0)

                -

                self.derecho.componentes.get(clave, 0)

            )

            if valor != 0:
                residuo[clave] = valor

        return residuo

    # --------------------------------------------------
    # Verificación
    # --------------------------------------------------

    def satisfecha(self):

        return len(self.residuo()) == 0

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (
            f"{self.nombre}"
            f"(satisfecha={self.satisfecha()})"
        )

    def __str__(self):

        return self.__repr__()