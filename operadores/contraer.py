"""
operadores/contraer.py
======================

Operadores de contracción tensorial del Motor de Geometría Diferencial (MGD).

Este módulo implementa la suma de Einstein de forma completamente
general para cualquier tensor definido sobre una variedad.

Autor:
    Elvis Omar Nazario Espinoza

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

from collections import defaultdict

import sympy as sp

from nucleo.tensor import Tensor
from nucleo.excepciones import (
    ErrorContraccion,
    ErrorDimension,
)


def validar_contraccion(
    tensor_a: Tensor,
    tensor_b: Tensor,
    posicion_a: int,
    posicion_b: int,
):
    """
    Verifica que una contracción sea matemáticamente válida.

    Parámetros
    ----------
    tensor_a : Tensor
        Primer tensor.

    tensor_b : Tensor
        Segundo tensor.

    posicion_a : int
        Posición del índice del primer tensor.

    posicion_b : int
        Posición del índice del segundo tensor.

    Lanza
    -----
    ErrorDimension
        Si las posiciones están fuera del rango permitido.

    ErrorContraccion
        Si los índices no son compatibles para la suma de Einstein.
    """

    if posicion_a < 0 or posicion_a >= tensor_a.rango:
        raise ErrorDimension(
            f"Posición {posicion_a} fuera del rango del tensor '{tensor_a.nombre}'."
        )

    if posicion_b < 0 or posicion_b >= tensor_b.rango:
        raise ErrorDimension(
            f"Posición {posicion_b} fuera del rango del tensor '{tensor_b.nombre}'."
        )

    nombre_a, posicion_indice_a = tensor_a.indices[posicion_a]
    nombre_b, posicion_indice_b = tensor_b.indices[posicion_b]

    if nombre_a != nombre_b:
        raise ErrorContraccion(
            "No es posible contraer índices con nombres diferentes "
            f"('{nombre_a}' y '{nombre_b}')."
        )

    if posicion_indice_a == posicion_indice_b:
        raise ErrorContraccion(
            "La contracción requiere un índice arriba y otro abajo."
        )

    if tensor_a.variedad.dimension != tensor_b.variedad.dimension:
        raise ErrorDimension(
            "Los tensores pertenecen a variedades de distinta dimensión."
        )

    return True
