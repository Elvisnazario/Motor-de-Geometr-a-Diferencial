"""
geodesica.py
============

Construcción simbólica de las ecuaciones geodésicas.

A partir de una conexión afín cualquiera se generan
automáticamente las ecuaciones

    d²x^μ/dλ²
    + Γ^μ_{αβ}
      dx^α/dλ
      dx^β/dλ
    = 0

Este módulo no asume Relatividad General.
Puede utilizarse con cualquier conexión compatible
con el Motor de Geometría Diferencial.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from conexion import Conexion


class Geodesica:
    """
    Constructor simbólico de geodésicas.
    """

    def __init__(self, conexion):

        if not isinstance(conexion, Conexion):
            raise TypeError(
                "Se esperaba un objeto Conexion."
            )

        self.conexion = conexion
        self.variedad = conexion.variedad

        self.dimension = self.variedad.dimension

        self.coords = self.variedad.coordenadas

        # Parámetro afín
        self.lambda_ = sp.symbols("lambda")

        # Coordenadas como funciones del parámetro
        self.x = [

            sp.Function(str(coord))(self.lambda_)

            for coord in self.coords

        ]

        self._sustituciones()

        self._construir()

    # --------------------------------------------------

    def _sustituciones(self):

        """
        Sustituye las coordenadas por x(λ)
        para que las Γ puedan escribirse
        sobre la trayectoria.
        """

        self.subs = {

            self.coords[i]: self.x[i]

            for i in range(self.dimension)

        }

    # --------------------------------------------------

    def _construir(self):

        self.ecuaciones = []

        for mu in range(self.dimension):

            aceleracion = sp.diff(

                self.x[mu],

                self.lambda_,

                2,

            )

            suma = 0

            for alpha in range(self.dimension):

                for beta in range(self.dimension):

                    gamma = (

                        self.conexion[

                            mu,

                            alpha,

                            beta,

                        ]

                    ).subs(self.subs)

                    velocidad_alpha = sp.diff(

                        self.x[alpha],

                        self.lambda_,

                    )

                    velocidad_beta = sp.diff(

                        self.x[beta],

                        self.lambda_,

                    )

                    suma += (

                        gamma

                        * velocidad_alpha

                        * velocidad_beta

                    )

            ecuacion = sp.Eq(

                aceleracion + sp.simplify(suma),

                0,

            )

            self.ecuaciones.append(ecuacion)

    # --------------------------------------------------

    def obtener(self):

        """
        Devuelve la lista de ecuaciones geodésicas.
        """

        return self.ecuaciones

    # --------------------------------------------------

    def latex(self):

        return [

            sp.latex(eq)

            for eq in self.ecuaciones

        ]

    # --------------------------------------------------

    def __getitem__(self, indice):

        return self.ecuaciones[indice]

    # --------------------------------------------------

    def __len__(self):

        return len(self.ecuaciones)

    # --------------------------------------------------

    def __repr__(self):

        return (

            "Geodesica("

            f"dimension={self.dimension})"

        )

    __str__ = __repr__