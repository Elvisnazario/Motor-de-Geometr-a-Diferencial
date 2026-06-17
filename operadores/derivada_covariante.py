"""
operadores/derivada_covariante.py
=================================

Operadores para calcular la derivada covariante de un tensor.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from nucleo.tensor import Tensor
from nucleo.excepciones import (
    ErrorDimension,
    ErrorIndicesIncompatibles,
)


def derivada_covariante(tensor, conexion, coordenada):
    """
    Calcula la derivada covariante de un tensor.

    Parámetros
    ----------
    tensor : Tensor
        Tensor sobre el cual se calcula la derivada.

    conexion : Conexion
        Conexión afín (Christoffel).

    coordenada : int
        Índice de la coordenada respecto a la cual se deriva.

    Retorna
    -------
    Tensor

    Notas
    -----
    En esta primera versión (v0.1) se implementa únicamente
    la derivada parcial de las componentes.

    Los términos de Christoffel se incorporarán cuando la
    clase Conexion quede completamente estabilizada.
    """

    if not isinstance(tensor, Tensor):
        raise TypeError(
            "El primer argumento debe ser un Tensor."
        )

    if coordenada < 0:
        raise ErrorDimension(
            "La coordenada no puede ser negativa."
        )

    if coordenada >= tensor.dimension:
        raise ErrorDimension(
            "La coordenada está fuera del rango de la variedad."
        )

    nuevas_componentes = {}

    simbolo = tensor.variedad.coordenadas[coordenada]

    for clave, valor in tensor.componentes.items():

        nuevas_componentes[clave] = sp.diff(
            valor,
            simbolo,
        )

    nuevos_indices = list(tensor.indices)

    nuevos_indices.append(
        (
            f"∂{simbolo}",
            "abajo",
        )
    )

    return Tensor(
        nombre=f"∇_{simbolo}({tensor.nombre})",
        variedad=tensor.variedad,
        componentes=nuevas_componentes,
        indices=nuevos_indices,
    )


def derivada_parcial(tensor, coordenada):
    """
    Alias de derivada_covariante mientras la conexión
    no agregue todavía los términos correctivos.
    """

    return derivada_covariante(
        tensor=tensor,
        conexion=None,
        coordenada=coordenada,
)
