"""
resultado_experimento.py
========================

Contenedor oficial de resultados del
Motor de Geometría Diferencial (MGD).

Todos los experimentos deben devolver un objeto de esta clase
para garantizar un formato uniforme y facilitar comparaciones.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

from datetime import datetime


class ResultadoExperimento:
    """
    Resultado uniforme de cualquier experimento del MGD.
    """

    def __init__(self, nombre):

        self.nombre = nombre

        self.fecha = datetime.now()

        self.metricas = {}

        self.observaciones = []

        self.conclusion = ""

    # --------------------------------------------------

    def agregar_metrica(self, nombre, valor):

        self.metricas[nombre] = valor

    # --------------------------------------------------

    def agregar_observacion(self, texto):

        self.observaciones.append(texto)

    # --------------------------------------------------

    def establecer_conclusion(self, texto):

        self.conclusion = texto

    # --------------------------------------------------

    def obtener(self, nombre):

        return self.metricas.get(nombre)

    # --------------------------------------------------

    def resumen(self):

        return {

            "nombre": self.nombre,

            "fecha": self.fecha,

            "metricas": self.metricas,

            "observaciones": self.observaciones,

            "conclusion": self.conclusion,

        }

    # --------------------------------------------------

    def __repr__(self):

        return (

            f"ResultadoExperimento("
            f"{self.nombre}, "
            f"{len(self.metricas)} métricas)"

        )

    __str__ = __repr__