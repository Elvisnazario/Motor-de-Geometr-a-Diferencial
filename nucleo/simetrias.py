"""
nucleo/simetrias.py
===================

Verificación de simetrías tensoriales.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from nucleo.tensor import Tensor


class Simetrias:
    """
    Rutinas para verificar propiedades de simetría
    de cualquier tensor.
    """

    @staticmethod
    def simetrico(tensor, indice1, indice2):
        """
        Verifica

            T(...i...j...) = T(...j...i...)
        """

        if not isinstance(tensor, Tensor):
            raise TypeError(
                "Se esperaba un Tensor."
            )

        for clave, valor in tensor.componentes.items():

            clave2 = list(clave)
            clave2[indice1], clave2[indice2] = (
                clave2[indice2],
                clave2[indice1],
            )

            valor2 = tensor[tuple(clave2)]

            if sp.simplify(valor - valor2) != 0:
                return False

        return True

    @staticmethod
    def antisimetrico(tensor, indice1, indice2):
        """
        Verifica

            T(...i...j...) = -T(...j...i...)
        """

        if not isinstance(tensor, Tensor):
            raise TypeError(
                "Se esperaba un Tensor."
            )

        for clave, valor in tensor.componentes.items():

            clave2 = list(clave)
            clave2[indice1], clave2[indice2] = (
                clave2[indice2],
                clave2[indice1],
            )

            valor2 = tensor[tuple(clave2)]

            if sp.simplify(valor + valor2) != 0:
                return False

        return True

    @staticmethod
    def intercambio_doble(
        tensor,
        par1,
        par2,
    ):
        """
        Verifica

            T(abcd)=T(cdab)

        Muy utilizado para Riemann.
        """

        if not isinstance(tensor, Tensor):
            raise TypeError(
                "Se esperaba un Tensor."
            )

        i, j = par1
        k, l = par2

        for clave, valor in tensor.componentes.items():

            nueva = list(clave)

            nueva[i], nueva[k] = nueva[k], nueva[i]
            nueva[j], nueva[l] = nueva[l], nueva[j]

            valor2 = tensor[tuple(nueva)]

            if sp.simplify(valor - valor2) != 0:
                return False

        return True

    @staticmethod
    def resumen(tensor):
        """
        Devuelve un pequeño informe de simetrías.
        """

        resultado = {}

        r = tensor.rango

        for i in range(r):

            for j in range(i + 1, r):

                resultado[f"{i}-{j}"] = {

                    "simetrico":
                    Simetrias.simetrico(
                        tensor,
                        i,
                        j,
                    ),

                    "antisimetrico":
                    Simetrias.antisimetrico(
                        tensor,
                        i,
                        j,
                    )

                }

        return resultado