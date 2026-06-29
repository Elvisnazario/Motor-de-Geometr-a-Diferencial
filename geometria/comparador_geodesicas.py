"""
comparador_geodesicas.py
========================

Herramientas para comparar trayectorias obtenidas mediante
diferentes métricas o diferentes conexiones.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import numpy as np


class ComparadorGeodesicas:
    """
    Compara dos soluciones numéricas de geodésicas.
    """

    def __init__(self, solucion_a, solucion_b):

        self.A = solucion_a
        self.B = solucion_b

    # --------------------------------------------------

    def distancia_maxima(self):

        if len(self.A.t) != len(self.B.t):
            raise ValueError(
                "Las soluciones deben tener el mismo mallado temporal."
            )

        diferencia = self.A.y - self.B.y

        norma = np.linalg.norm(
            diferencia,
            axis=0,
        )

        return np.max(norma)

    # --------------------------------------------------

    def distancia_media(self):

        diferencia = self.A.y - self.B.y

        norma = np.linalg.norm(
            diferencia,
            axis=0,
        )

        return np.mean(norma)

    # --------------------------------------------------

    def error_final(self):

        diferencia = (

            self.A.y[:, -1]
            -
            self.B.y[:, -1]

        )

        return np.linalg.norm(diferencia)

    # --------------------------------------------------

    def punto_mayor_divergencia(self):

        diferencia = self.A.y - self.B.y

        norma = np.linalg.norm(
            diferencia,
            axis=0,
        )

        indice = np.argmax(norma)

        return (

            self.A.t[indice],

            norma[indice],

        )

    # --------------------------------------------------

    def resumen(self):

        instante, divergencia = (
            self.punto_mayor_divergencia()
        )

        return {

            "distancia_maxima":
                self.distancia_maxima(),

            "distancia_media":
                self.distancia_media(),

            "error_final":
                self.error_final(),

            "instante_maximo":
                instante,

            "divergencia_maxima":
                divergencia,

        }

    # --------------------------------------------------

    def __repr__(self):

        return (
            "ComparadorGeodesicas()"
        )

    __str__ = __repr__