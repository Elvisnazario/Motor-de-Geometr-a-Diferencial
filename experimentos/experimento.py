"""
experimento.py
==============

Clase base para todos los experimentos del
Motor de Geometría Diferencial (MGD).

Todo experimento debe ser completamente reproducible,
independiente de la teoría gravitacional utilizada.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

from abc import ABC, abstractmethod
import time


class Experimento(ABC):
    """
    Clase base para cualquier experimento del MGD.
    """

    def __init__(self, nombre):

        self.nombre = nombre

        self.resultados = {}

        self.tiempo_inicio = None
        self.tiempo_fin = None

    # --------------------------------------------------

    def ejecutar(self):

        """
        Ejecuta el experimento completo.
        """

        self.tiempo_inicio = time.perf_counter()

        self.preparar()

        self.simular()

        self.analizar()

        self.tiempo_fin = time.perf_counter()

        self.resultados["tiempo_ejecucion"] = (
            self.tiempo_fin
            - self.tiempo_inicio
        )

        return self.resultados

    # --------------------------------------------------

    @abstractmethod
    def preparar(self):
        """
        Configuración inicial.
        """
        pass

    # --------------------------------------------------

    @abstractmethod
    def simular(self):
        """
        Ejecuta la simulación.
        """
        pass

    # --------------------------------------------------

    @abstractmethod
    def analizar(self):
        """
        Analiza los resultados.
        """
        pass

    # --------------------------------------------------

    def resumen(self):

        return self.resultados

    # --------------------------------------------------

    def __repr__(self):

        return (
            f"Experimento("
            f"{self.nombre})"
        )

    __str__ = __repr__