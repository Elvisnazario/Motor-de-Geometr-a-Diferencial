"""
bianchi.py
==========

Verificación de la identidad contraída de Bianchi.

    ∇^μ G_{μν} = 0

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from einstein import Einstein
from operadores.derivada_covariante import derivada_covariante


class Bianchi:
    """
    Verificador de la identidad contraída de Bianchi.

    Calcula

        ∇^μ G_{μν}

    utilizando la derivada covariante del tensor de Einstein.
    """

    def __init__(self, einstein, conexion):

        if not isinstance(einstein, Einstein):
            raise TypeError(
                "Se esperaba un objeto Einstein."
            )

        self.einstein = einstein
        self.conexion = conexion
        self.variedad = einstein.variedad
        self.dimension = self.variedad.dimension

        self._calcular()

    # --------------------------------------------------

    def _calcular(self):

        self.divergencia = {}

        for mu in range(self.dimension):

            derivada = derivada_covariante(
                self.einstein,
                self.conexion,
                mu,
            )

            for clave, valor in derivada.componentes.items():

                #
                # clave esperada:
                #
                # (μ, ν, λ)
                #
                # Contraemos λ con μ
                #

                if len(clave) != 3:
                    continue

                indice_mu = clave[0]
                indice_nu = clave[1]
                indice_lambda = clave[2]

                if indice_mu != indice_lambda:
                    continue

                self.divergencia[indice_nu] = (
                    self.divergencia.get(indice_nu, 0)
                    + valor
                )

        for clave in list(self.divergencia.keys()):
            self.divergencia[clave] = sp.simplify(
                self.divergencia[clave]
            )

    # --------------------------------------------------

    def satisfecha(self):

        for valor in self.divergencia.values():

            if sp.simplify(valor) != 0:
                return False

        return True

    # --------------------------------------------------

    def residuo(self):

        return self.divergencia

    # --------------------------------------------------

    def __repr__(self):

        return (
            f"Bianchi("
            f"satisfecha={self.satisfecha()})"
        )

    def __str__(self):

        return self.__repr__()