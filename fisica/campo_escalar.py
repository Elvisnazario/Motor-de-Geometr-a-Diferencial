"""
campo_escalar.py
================

Implementación de un campo escalar clásico sobre una variedad.

Esta clase constituye la base para:

- Klein-Gordon
- Inflación cosmológica
- Quintesencia
- Campos escalares modificados
- Campo escalar recursivo del MGD

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from fisica.campo import CampoFisico


class CampoEscalar(CampoFisico):
    """
    Campo escalar φ definido sobre una variedad.
    """

    def __init__(
        self,
        nombre,
        variedad,
        simbolo="phi",
    ):

        super().__init__(nombre, variedad)

        self.funcion = sp.Function(simbolo)(
            *variedad.coordenadas
        )

    # --------------------------------------------------

    @property
    def phi(self):

        return self.funcion

    # --------------------------------------------------

    def gradiente(self):

        """
        Devuelve

            ∂μφ
        """

        return [

            sp.diff(
                self.funcion,
                coord,
            )

            for coord in self.variedad.coordenadas

        ]

    # --------------------------------------------------

    def hessiano(self):

        """
        Devuelve

            ∂μ∂νφ
        """

        dim = self.variedad.dimension

        H = sp.MutableDenseMatrix.zeros(
            dim,
            dim,
        )

        coords = self.variedad.coordenadas

        for mu in range(dim):

            for nu in range(dim):

                H[mu, nu] = sp.diff(

                    self.funcion,

                    coords[mu],

                    coords[nu],

                )

        return H

    # --------------------------------------------------

    def laplaciano(self, metrica):

        """
        Calcula

            □φ = g^{μν} ∂μ∂νφ

        usando la métrica inversa.
        """

        ginv = metrica.inversa()

        H = self.hessiano()

        dim = self.variedad.dimension

        resultado = 0

        for mu in range(dim):

            for nu in range(dim):

                resultado += (

                    ginv[mu, nu]

                    * H[mu, nu]

                )

        return sp.simplify(resultado)

    # --------------------------------------------------

    def lagrangiano(self):

        """
        La densidad lagrangiana será implementada
        por cada teoría particular.
        """

        raise NotImplementedError

    # --------------------------------------------------

    def tensor_energia(self):

        raise NotImplementedError

    # --------------------------------------------------

    def ecuaciones(self):

        raise NotImplementedError

    # --------------------------------------------------

    def __repr__(self):

        return (

            f"CampoEscalar("

            f"{self.nombre})"

        )

    __str__ = __repr__