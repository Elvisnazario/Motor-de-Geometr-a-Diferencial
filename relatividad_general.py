"""
relatividad_general.py
======================

Implementación de la Relatividad General utilizando el
Motor de Geometría Diferencial (MGD).

Esta clase sirve como referencia para verificar que el
motor reproduce exactamente las ecuaciones clásicas de
Einstein antes de incorporar teorías modificadas.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

from teoria_campo import TeoriaCampo

from einstein import Einstein

from tensor_energia_momento import TensorEnergiaMomento


class RelatividadGeneral(TeoriaCampo):
    """
    Implementación de

        G_{μν} = κ T_{μν}

    donde

        κ = 8πG/c⁴

    permanece simbólico.
    """

    def __init__(
        self,
        variedad,
        riemann,
        metrica,
        tensor_energia,
        kappa,
    ):

        super().__init__(variedad)

        self.riemann = riemann
        self.metrica = metrica
        self.tensor_energia = tensor_energia
        self.kappa = kappa

        self._einstein = Einstein.desde_riemann(
            riemann,
            metrica,
        )

    # --------------------------------------------------
    # Lado geométrico
    # --------------------------------------------------

    def miembro_geometrico(self):

        return self._einstein

    # --------------------------------------------------
    # Lado físico
    # --------------------------------------------------

    def miembro_fisico(self):

        return self.kappa * self.tensor_energia

    # --------------------------------------------------
    # Accesos rápidos
    # --------------------------------------------------

    @property
    def G(self):

        return self._einstein

    @property
    def T(self):

        return self.tensor_energia

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (
            "RelatividadGeneral("
            f"dimension={self.variedad.dimension})"
        )

    def __str__(self):

        return self.__repr__()