"""
nucleo/cache.py
===============

Sistema de caché del Motor de Geometría Diferencial (MGD).

Su objetivo es evitar el recálculo de cantidades geométricas
costosas como:

    • métricas inversas
    • Christoffel
    • Riemann
    • Ricci
    • Escalar de Ricci
    • Einstein
    • derivadas covariantes
    • identidades de Bianchi

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import hashlib


class CacheMGD:
    """
    Caché central del motor.
    """

    def __init__(self):

        self._datos = {}

    # -------------------------------------------------

    def limpiar(self):
        """
        Vacía completamente el caché.
        """

        self._datos.clear()

    # -------------------------------------------------

    def existe(self, clave):

        return clave in self._datos

    # -------------------------------------------------

    def guardar(self, clave, objeto):

        self._datos[clave] = objeto

    # -------------------------------------------------

    def obtener(self, clave):

        return self._datos.get(clave)

    # -------------------------------------------------

    def eliminar(self, clave):

        if clave in self._datos:
            del self._datos[clave]

    # -------------------------------------------------

    def total_objetos(self):

        return len(self._datos)

    # -------------------------------------------------

    @staticmethod
    def generar_clave(*objetos):
        """
        Genera una clave única a partir
        de cualquier combinación de objetos.
        """

        texto = ""

        for obj in objetos:

            texto += repr(obj)

            texto += "|"

        return hashlib.sha256(
            texto.encode("utf8")
        ).hexdigest()

    # -------------------------------------------------

    def __repr__(self):

        return (
            f"CacheMGD("
            f"objetos={len(self._datos)})"
        )

    __str__ = __repr__


# ---------------------------------------------------------
# Caché global del motor
# ---------------------------------------------------------

CACHE = CacheMGD()