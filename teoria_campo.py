"""
teoria_campo.py
===============

Clase base para cualquier teoría de campos construida
sobre el Motor de Geometría Diferencial (MGD).

Una teoría física queda definida por:

    Geometría = Fuente

sin imponer una teoría concreta.

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

from abc import ABC, abstractmethod

from ecuacion_campo import EcuacionCampo


class TeoriaCampo(ABC):
    """
    Clase base abstracta para cualquier teoría.

    Ejemplos:

        - Relatividad General
        - MGD
        - Teorías f(R)
        - Einstein-Cartan
        - Yang-Mills geométrico
        - etc.
    """

    def __init__(self, variedad):

        self.variedad = variedad

    # --------------------------------------------------
    # Miembro geométrico
    # --------------------------------------------------

    @abstractmethod
    def miembro_geometrico(self):
        """
        Devuelve el tensor geométrico de la teoría.
        """
        pass

    # --------------------------------------------------
    # Miembro físico
    # --------------------------------------------------

    @abstractmethod
    def miembro_fisico(self):
        """
        Devuelve el tensor fuente.
        """
        pass

    # --------------------------------------------------
    # Construcción automática
    # --------------------------------------------------

    def ecuacion(self):

        return EcuacionCampo(

            self.miembro_geometrico(),

            self.miembro_fisico(),

            nombre=self.__class__.__name__,

        )

    # --------------------------------------------------
    # Verificación
    # --------------------------------------------------

    def satisfecha(self):

        return self.ecuacion().satisfecha()

    # --------------------------------------------------
    # Residuo
    # --------------------------------------------------

    def residuo(self):

        return self.ecuacion().residuo()

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (

            f"{self.__class__.__name__}"

            f"(dimension={self.variedad.dimension})"

        )

    def __str__(self):

        return self.__repr__()