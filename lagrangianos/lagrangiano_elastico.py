"""
lagrangiano_elastico.py
=======================

Lagrangiano fundamental del Medio Geométrico Elástico.

Este archivo constituye la base variacional de todas las
teorías construidas sobre el Motor de Geometría Diferencial.

No implementa una teoría particular.

Describe únicamente la estructura general

    L = L_gravedad
      + L_medio
      + L_materia
      + L_acoplamiento

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp


class LagrangianoElastico:
    """
    Lagrangiano general del Medio Geométrico.

    Cada término puede activarse o desactivarse
    independientemente.
    """

    def __init__(self):

        self.terminos = {}

    # --------------------------------------------------

    def agregar(
        self,
        nombre,
        expresion,
    ):

        self.terminos[nombre] = sp.expand(expresion)

    # --------------------------------------------------

    def eliminar(
        self,
        nombre,
    ):

        if nombre in self.terminos:
            del self.terminos[nombre]

    # --------------------------------------------------

    def obtener(
        self,
        nombre,
    ):

        return self.terminos.get(
            nombre,
            sp.Integer(0),
        )

    # --------------------------------------------------

    @property
    def expresion(self):

        total = sp.Integer(0)

        for termino in self.terminos.values():
            total += termino

        return sp.simplify(total)

    # --------------------------------------------------

    @property
    def gravedad(self):

        return self.obtener("gravedad")

    @property
    def medio(self):

        return self.obtener("medio")

    @property
    def materia(self):

        return self.obtener("materia")

    @property
    def acoplamiento(self):

        return self.obtener("acoplamiento")

    # --------------------------------------------------

    def activar_gravedad(
        self,
        escalar_ricci,
        constante=1,
    ):
        """
        L_gravedad

            constante * R
        """

        self.agregar(
            "gravedad",
            constante * escalar_ricci,
        )

    # --------------------------------------------------

    def activar_medio(
        self,
        densidad_energia,
    ):
        """
        L_medio

        Energía interna del medio geométrico.
        """

        self.agregar(
            "medio",
            densidad_energia,
        )

    # --------------------------------------------------

    def activar_materia(
        self,
        lagrangiano_materia,
    ):

        self.agregar(
            "materia",
            lagrangiano_materia,
        )

    # --------------------------------------------------

    def activar_acoplamiento(
        self,
        termino,
    ):

        self.agregar(
            "acoplamiento",
            termino,
        )

    # --------------------------------------------------

    def __repr__(self):

        return (
            f"LagrangianoElastico("
            f"terminos={list(self.terminos.keys())})"
        )

    __str__ = __repr__