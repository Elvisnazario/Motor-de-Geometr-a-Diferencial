"""
nucleo/tensor.py
================

Implementación del objeto Tensor del Motor de Geometría Diferencial (MGD).

El Tensor constituye la estructura algebraica fundamental del motor.
Su única responsabilidad es almacenar componentes de forma dispersa
y garantizar la consistencia matemática de sus índices.

Las operaciones geométricas (contracción, derivadas, subida/bajada
de índices, etc.) pertenecen al paquete operadores/.

Autor:
    Elvis Omar Nazario Espinoza

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

from collections import defaultdict
import sympy as sp

from nucleo.objeto_geometrico import ObjetoGeometrico
from nucleo.excepciones import (
    ErrorDimension,
    ErrorIndicesIncompatibles,
)


class Tensor(ObjetoGeometrico):
    """
    Representación abstracta de un tensor.

    Parámetros
    ----------
    nombre : str

    variedad : Variedad

    componentes : dict

    indices : list
        Lista de tuplas

            ("μ","arriba")

            ("ν","abajo")
    """

    def __init__(
        self,
        nombre,
        variedad,
        componentes,
        indices,
    ):

        super().__init__(
            nombre=nombre,
            variedad=variedad,
            indices=indices,
        )

        self.componentes = defaultdict(
            lambda: 0,
            componentes if componentes is not None else {},
        )

        self._validar_componentes()

    # ==========================================================
    # PROPIEDADES
    # ==========================================================

    @property
    def rango(self):
        """Número de índices del tensor."""
        return len(self.indices)

    @property
    def dimension(self):
        """Dimensión de la variedad."""
        return self.variedad.dimension

    @property
    def forma(self):
        """
        Devuelve únicamente la estructura de covarianza.

        Ejemplo

            ('arriba','abajo','abajo')
        """

        return tuple(
            posicion
            for _, posicion in self.indices
        )

    # ==========================================================
    # VALIDACIÓN
    # ==========================================================

    def _validar_componentes(self):

        for clave in self.componentes.keys():

            if len(clave) != self.rango:

                raise ErrorDimension(
                    f"La componente {clave} posee "
                    f"{len(clave)} índices, "
                    f"pero el tensor '{self.nombre}' "
                    f"es de rango {self.rango}."
                )

            for indice in clave:

                if not isinstance(indice, int):

                    raise ErrorDimension(
                        f"El índice '{indice}' "
                        f"no es un entero."
                    )

                if indice < 0:

                    raise ErrorDimension(
                        f"Índice negativo "
                        f"en {clave}."
                    )

                if indice >= self.dimension:

                    raise ErrorDimension(
                        f"Índice {indice} fuera "
                        f"del rango permitido "
                        f"0..{self.dimension-1}."
                    )

    # ==========================================================
    # ACCESO
    # ==========================================================

    def __getitem__(self, clave):

        return self.componentes[tuple(clave)]

    def __setitem__(self, clave, valor):

        clave = tuple(clave)

        if len(clave) != self.rango:

            raise ErrorDimension(
                "La longitud de la clave "
                "no coincide con el rango "
                "del tensor."
            )

        self.componentes[clave] = valor

    # ==========================================================
    # UTILIDADES
    # ==========================================================

    def copiar(self):
        """
        Devuelve una copia profunda del tensor.
        """

        return Tensor(
            nombre=self.nombre,
            variedad=self.variedad,
            componentes=dict(self.componentes),
            indices=list(self.indices),
        )

            # ==========================================================
    # OPERADORES ALGEBRAICOS
    # ==========================================================

    def __add__(self, otro):

        if not isinstance(otro, Tensor):
            return NotImplemented

        if self.dimension != otro.dimension:
            raise ErrorDimension(
                "No pueden sumarse tensores de distinta dimensión."
            )

        if self.indices != otro.indices:
            raise ErrorIndicesIncompatibles(
                f"No puede sumarse '{self.nombre}' con "
                f"'{otro.nombre}' porque sus índices "
                "no poseen la misma estructura."
            )

        componentes = {}

        claves = (
            set(self.componentes.keys())
            | set(otro.componentes.keys())
        )

        for clave in claves:
            componentes[clave] = (
                self[clave] + otro[clave]
            )

        return Tensor(
            nombre=f"({self.nombre}+{otro.nombre})",
            variedad=self.variedad,
            componentes=componentes,
            indices=list(self.indices),
        )

    def __sub__(self, otro):

        if not isinstance(otro, Tensor):
            return NotImplemented

        if self.dimension != otro.dimension:
            raise ErrorDimension(
                "No pueden restarse tensores de distinta dimensión."
            )

        if self.indices != otro.indices:
            raise ErrorIndicesIncompatibles(
                f"No puede restarse '{self.nombre}' con "
                f"'{otro.nombre}' porque sus índices "
                "no poseen la misma estructura."
            )

        componentes = {}

        claves = (
            set(self.componentes.keys())
            | set(otro.componentes.keys())
        )

        for clave in claves:
            componentes[clave] = (
                self[clave] - otro[clave]
            )

        return Tensor(
            nombre=f"({self.nombre}-{otro.nombre})",
            variedad=self.variedad,
            componentes=componentes,
            indices=list(self.indices),
        )

    def __mul__(self, escalar):

        if isinstance(escalar, Tensor):
            return NotImplemented

        componentes = {}

        for clave, valor in self.componentes.items():
            componentes[clave] = escalar * valor

        return Tensor(
            nombre=f"{escalar}·{self.nombre}",
            variedad=self.variedad,
            componentes=componentes,
            indices=list(self.indices),
        )

    __rmul__ = __mul__

    # ==========================================================
    # COMPARACIÓN
    # ==========================================================

    def __eq__(self, otro):

        if not isinstance(otro, Tensor):
            return False

        return (
            self.indices == otro.indices
            and self.dimension == otro.dimension
            and dict(self.componentes) == dict(otro.componentes)
        )

    # ==========================================================
    # REPRESENTACIÓN
    # ==========================================================

    def __repr__(self):

        return (
            f"Tensor("
            f"nombre='{self.nombre}', "
            f"rango={self.rango}, "
            f"dimension={self.dimension}, "
            f"componentes={len(self.componentes)})"
        )

    def __str__(self):

        return self.__repr__()
