"""
accion_geodesica.py
===================

Cálculo de la acción asociada a una trayectoria.

Este módulo permite evaluar funcionales de acción sobre
cualquier trayectoria parametrizada.

No supone una teoría gravitacional específica.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import numpy as np


class AccionGeodesica:
    """
    Calcula la acción de una trayectoria.

    La acción se aproxima mediante integración numérica
    sobre los puntos entregados por el integrador.
    """

    def __init__(self, metrica):

        self.metrica = metrica

        self.dimension = metrica.dimension

    # --------------------------------------------------

    def longitud_propia(self, solucion):

        t = solucion.t
        y = solucion.y

        accion = 0.0

        for k in range(len(t) - 1):

            dt = t[k + 1] - t[k]

            posicion = y[:self.dimension, k]
            velocidad = y[self.dimension:, k]

            g = self.metrica.g_cov

            ds2 = 0

            for mu in range(self.dimension):
                for nu in range(self.dimension):

                    ds2 += (
                        g[mu, nu]
                        * velocidad[mu]
                        * velocidad[nu]
                    )

            accion += np.sqrt(abs(float(ds2))) * dt

        return accion

    # --------------------------------------------------

    def accion(self, solucion):

        return self.longitud_propia(solucion)

    # --------------------------------------------------

    def comparar(self, solucion_a, solucion_b):

        S1 = self.accion(solucion_a)
        S2 = self.accion(solucion_b)

        return {

            "accion_A": S1,
            "accion_B": S2,
            "diferencia": S1 - S2,
            "minima": "A" if S1 < S2 else "B"

        }

    # --------------------------------------------------

    def __repr__(self):

        return (
            "AccionGeodesica("
            f"dimension={self.dimension})"
        )

    __str__ = __repr__