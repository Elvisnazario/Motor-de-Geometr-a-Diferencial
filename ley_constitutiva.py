"""
ley_constitutiva.py
===================

Ley constitutiva del Medio Geométrico.

Este módulo define la respuesta constitutiva del
espacio-tiempo dentro del Motor de Geometría Diferencial
(MGD).

Toda interacción física modifica el medio mediante una
densidad de energía efectiva.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp


class LeyConstitutiva:
    """
    Modelo constitutivo del medio geométrico.

    Variables fundamentales

        ψ      Campo escalar.

        ω      Norma del gradiente.

        σ      Acoplamiento elástico.

        K0     Rigidez fundamental.

    A partir de ellas se calculan las cantidades
    físicas que alimentarán el Tensor de Energía
    Recursiva.
    """

    def __init__(

        self,

        psi,

        omega,

        sigma,

        K0,

    ):

        self.psi = sp.simplify(psi)

        self.omega = sp.simplify(omega)

        self.sigma = sp.simplify(sigma)

        self.K0 = sp.simplify(K0)

    # --------------------------------------------------
    # Energía almacenada
    # --------------------------------------------------

    @property
    def rho(self):
        """
        Densidad efectiva de energía.

        Primera aproximación constitutiva del MGD.
        """

        return sp.simplify(

            self.K0

            *

            (1 - self.sigma)

            *

            self.omega

        )

    # --------------------------------------------------
    # Presión radial
    # --------------------------------------------------

    @property
    def p_r(self):

        return sp.simplify(

            self.sigma

            *

            self.rho

        )

    # --------------------------------------------------
    # Presión tangencial
    # --------------------------------------------------

    @property
    def p_t(self):

        return sp.simplify(

            sp.Rational(1, 2)

            *

            (1 + self.sigma)

            *

            self.rho

        )

    # --------------------------------------------------
    # Ecuaciones de estado
    # --------------------------------------------------

    @property
    def w_r(self):

        if self.rho == 0:
            return sp.nan

        return sp.simplify(

            self.p_r

            /

            self.rho

        )

    @property
    def w_t(self):

        if self.rho == 0:
            return sp.nan

        return sp.simplify(

            self.p_t

            /

            self.rho

        )

    # --------------------------------------------------
    # Estado del medio
    # --------------------------------------------------

    @property
    def energia_total(self):

        return sp.simplify(

            self.rho

            +

            self.p_r

            +

            2 * self.p_t

        )

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (

            "LeyConstitutiva("

            f"rho={self.rho}, "

            f"pr={self.p_r}, "

            f"pt={self.p_t})"

        )