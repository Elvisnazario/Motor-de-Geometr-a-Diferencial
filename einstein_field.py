"""
einstein_field.py
=================

Construcción del miembro geométrico de las ecuaciones de Einstein.

    G_{μν}

Este módulo encapsula el lado izquierdo de las ecuaciones
de campo y deja preparado el acoplamiento posterior con el
tensor energía-momento del MGD.

Autor:
    Elvis Omar Nazario Espinoza

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

from einstein import Einstein


class EinsteinField:
    """
    Miembro geométrico de las ecuaciones de campo.

        G_{μν}

    En futuras extensiones del MGD este objeto permitirá
    incorporar términos adicionales sin modificar la clase
    Einstein.
    """

    def __init__(self, einstein):

        if not isinstance(einstein, Einstein):
            raise TypeError(
                "Se esperaba un objeto Einstein."
            )

        self.einstein = einstein
        self.variedad = einstein.variedad

    # --------------------------------------------------
    # Acceso
    # --------------------------------------------------

    def __getitem__(self, clave):
        return self.einstein[clave]

    @property
    def componentes(self):
        return self.einstein.componentes

    @property
    def indices(self):
        return self.einstein.indices

    @property
    def dimension(self):
        return self.variedad.dimension

    # --------------------------------------------------
    # Utilidad
    # --------------------------------------------------

    def es_vacio(self):
        """
        Verifica si G_{μν}=0.
        """
        return self.einstein.es_nulo()

    # --------------------------------------------------
    # Representación
    # --------------------------------------------------

    def __repr__(self):

        return (
            "EinsteinField("
            f"dimension={self.dimension}, "
            f"componentes={len(self.componentes)})"
        )

    def __str__(self):
        return self.__repr__()