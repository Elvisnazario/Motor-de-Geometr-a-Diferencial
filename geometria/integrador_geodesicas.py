"""
integrador_geodesicas.py
========================

Integrador numérico de ecuaciones geodésicas.

Convierte las ecuaciones simbólicas obtenidas por la clase
Geodesica en un sistema de ecuaciones diferenciales de primer
orden listo para integración numérica.

No depende de una teoría gravitacional particular.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp
from scipy.integrate import solve_ivp

from geodesica import Geodesica


class IntegradorGeodesicas:

    def __init__(self, geodesica):

        if not isinstance(geodesica, Geodesica):
            raise TypeError(
                "Se esperaba un objeto Geodesica."
            )

        self.geodesica = geodesica
        self.dimension = geodesica.dimension
        self.lambda_ = geodesica.lambda_

        self.x = geodesica.x

        self._preparar()

    # --------------------------------------------------

    def _preparar(self):

        velocidades = [

            sp.Symbol(f"v{i}")

            for i in range(self.dimension)

        ]

        sustituciones = {}

        for i in range(self.dimension):

            sustituciones[
                sp.diff(
                    self.x[i],
                    self.lambda_
                )
            ] = velocidades[i]

        self.funciones = []

        for ecuacion in self.geodesica.obtener():

            aceleracion = sp.diff(
                self.x[
                    len(self.funciones)
                ],
                self.lambda_,
                2
            )

            rhs = sp.solve(
                ecuacion,
                aceleracion
            )[0]

            rhs = rhs.subs(sustituciones)

            argumentos = (
                self.x +
                velocidades
            )

            self.funciones.append(

                sp.lambdify(
                    argumentos,
                    rhs,
                    "numpy"
                )

            )

    # --------------------------------------------------

    def integrar(

        self,

        intervalo,

        condiciones_iniciales,

        pasos=500,

    ):

        dim = self.dimension

        if len(condiciones_iniciales) != 2 * dim:

            raise ValueError(
                "Número incorrecto de condiciones iniciales."
            )

        def sistema(lmbda, Y):

            posiciones = list(Y[:dim])

            velocidades = list(Y[dim:])

            derivadas = velocidades

            aceleraciones = []

            argumentos = posiciones + velocidades

            for funcion in self.funciones:

                aceleraciones.append(

                    funcion(*argumentos)

                )

            return derivadas + aceleraciones

        solucion = solve_ivp(

            sistema,

            intervalo,

            condiciones_iniciales,

            dense_output=True,

            max_step=(
                intervalo[1] - intervalo[0]
            ) / pasos,

        )

        return solucion

    # --------------------------------------------------

    def __repr__(self):

        return (
            "IntegradorGeodesicas("
            f"dimension={self.dimension})"
        )

    __str__ = __repr__