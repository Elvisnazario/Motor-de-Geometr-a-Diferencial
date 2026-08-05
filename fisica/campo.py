"""
campo.py
========

Clase base para todos los campos físicos del
Motor de Geometría Diferencial (MGD).

Autor:
    Elvis Omar Nazario Espinoza
"""

from abc import ABC, abstractmethod


class CampoFisico(ABC):
    """
    Clase base abstracta para cualquier campo físico
    definido sobre una variedad diferenciable.
    """

    def __init__(self, nombre, variedad):

        self.nombre = nombre
        self.variedad = variedad

    # --------------------------------------------------

    @abstractmethod
    def lagrangiano(self):
        """
        Devuelve la densidad lagrangiana del campo.
        """
        pass

    # --------------------------------------------------

    @abstractmethod
    def tensor_energia(self):
        """
        Devuelve el tensor energía-momento.
        """
        pass

    # --------------------------------------------------

    @abstractmethod
    def ecuaciones(self):
        """
        Devuelve las ecuaciones de campo.
        """
        pass

    # --------------------------------------------------

    def __repr__(self):

        return f"{self.__class__.__name__}(nombre={self.nombre})"

    __str__ = __repr__