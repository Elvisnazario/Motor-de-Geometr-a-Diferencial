"""
nucleo/identidades.py
=====================

Verificación automática de identidades tensoriales.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from nucleo.tensor import Tensor


class Identidades:
    """
    Rutinas para comprobar identidades entre tensores.
    """

    @staticmethod
    def iguales(tensor1, tensor2):
        """
        Comprueba si dos tensores son idénticos.
        """

        if not isinstance(tensor1, Tensor):
            raise TypeError(
                "tensor1 debe ser Tensor."
            )

        if not isinstance(tensor2, Tensor):
            raise TypeError(
                "tensor2 debe ser Tensor."
            )

        if tensor1.indices != tensor2.indices:
            return False

        claves = (
            set(tensor1.componentes.keys())
            |
            set(tensor2.componentes.keys())
        )

        for clave in claves:

            v1 = tensor1.componentes.get(clave, 0)
            v2 = tensor2.componentes.get(clave, 0)

            if sp.simplify(v1 - v2) != 0:
                return False

        return True

    @staticmethod
    def diferencia(tensor1, tensor2):
        """
        Devuelve únicamente las componentes
        donde ambos tensores difieren.
        """

        diferencias = {}

        claves = (
            set(tensor1.componentes.keys())
            |
            set(tensor2.componentes.keys())
        )

        for clave in claves:

            valor = sp.simplify(

                tensor1.componentes.get(clave, 0)

                -

                tensor2.componentes.get(clave, 0)

            )

            if valor != 0:
                diferencias[clave] = valor

        return diferencias

    @staticmethod
    def es_nulo(tensor):
        """
        Comprueba si un tensor es exactamente cero.
        """

        if not isinstance(tensor, Tensor):
            raise TypeError(
                "Se esperaba un Tensor."
            )

        for valor in tensor.componentes.values():

            if sp.simplify(valor) != 0:
                return False

        return True

    @staticmethod
    def imprimir_diferencias(tensor1, tensor2):
        """
        Imprime únicamente las diferencias.
        """

        diferencias = Identidades.diferencia(
            tensor1,
            tensor2,
        )

        if len(diferencias) == 0:

            print(
                "Las identidades coinciden exactamente."
            )

            return

        print("Diferencias encontradas:")

        for clave, valor in diferencias.items():

            print(
                f"{clave} : {valor}"
            )