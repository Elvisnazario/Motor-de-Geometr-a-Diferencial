"""
parametros_fundamentales.py
===========================

Administrador de parámetros fundamentales para teorías
construidas sobre el Motor de Geometría Diferencial (MGD).

El MGD no asume el significado físico de ningún parámetro.
Cada teoría es libre de interpretarlos.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp


class ParametrosFundamentales:
    """
    Registro simbólico de constantes y parámetros.

    Ejemplos
    --------
    r_min
    K0
    lambda
    alpha
    beta
    """

    def __init__(self):

        self._parametros = {}

    # --------------------------------------------------

    def registrar(
        self,
        nombre,
        valor=None,
        descripcion="",
    ):

        simbolo = sp.Symbol(
            nombre,
            real=True,
        )

        self._parametros[nombre] = {
            "simbolo": simbolo,
            "valor": valor,
            "descripcion": descripcion,
        }

        return simbolo

    # --------------------------------------------------

    def simbolo(
        self,
        nombre,
    ):

        return self._parametros[nombre]["simbolo"]

    # --------------------------------------------------

    def valor(
        self,
        nombre,
    ):

        return self._parametros[nombre]["valor"]

    # --------------------------------------------------

    def fijar_valor(
        self,
        nombre,
        valor,
    ):

        self._parametros[nombre]["valor"] = valor

    # --------------------------------------------------

    def descripcion(
        self,
        nombre,
    ):

        return self._parametros[nombre]["descripcion"]

    # --------------------------------------------------

    def sustituciones(self):

        reglas = {}

        for parametro in self._parametros.values():

            if parametro["valor"] is not None:

                reglas[
                    parametro["simbolo"]
                ] = parametro["valor"]

        return reglas

    # --------------------------------------------------

    def nombres(self):

        return list(
            self._parametros.keys()
        )

    # --------------------------------------------------

    def __repr__(self):

        return (
            f"ParametrosFundamentales("
            f"{self.nombres()})"
        )

    __str__ = __repr__