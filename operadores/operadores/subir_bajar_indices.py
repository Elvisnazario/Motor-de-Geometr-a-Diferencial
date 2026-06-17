"""
operadores/subir_bajar_indices.py
=================================

Operadores para subir y bajar índices utilizando la métrica.

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

from collections import defaultdict

from nucleo.tensor import Tensor
from nucleo.excepciones import (
    ErrorDimension,
)


def subir_indice(tensor, metrica_inversa, posicion):
    """
    Sube un índice covariante usando la métrica inversa.

    Parámetros
    ----------
    tensor : Tensor

    metrica_inversa : Tensor

    posicion : int
        Posición del índice a elevar.

    Nota
    ----
    Esta primera versión implementa únicamente la validación
    estructural. La suma completa de Einstein se incorporará
    cuando el módulo de contracción quede estabilizado.
    """

    if posicion < 0 or posicion >= tensor.rango:
        raise ErrorDimension(
            "La posición del índice está fuera del rango."
        )

    nombre, tipo = tensor.indices[posicion]

    if tipo != "abajo":
        raise ValueError(
            "Solo pueden elevarse índices covariantes."
        )

    nuevos_indices = list(tensor.indices)
    nuevos_indices[posicion] = (nombre, "arriba")

    return Tensor(
        nombre=f"{tensor.nombre}↑",
        variedad=tensor.variedad,
        componentes=dict(tensor.componentes),
        indices=nuevos_indices,
    )


def bajar_indice(tensor, metrica, posicion):
    """
    Baja un índice contravariante usando la métrica.

    Esta versión conserva las componentes y modifica únicamente
    la estructura algebraica del tensor. La contracción completa
    será implementada posteriormente.
    """

    if posicion < 0 or posicion >= tensor.rango:
        raise ErrorDimension(
            "La posición del índice está fuera del rango."
        )

    nombre, tipo = tensor.indices[posicion]

    if tipo != "arriba":
        raise ValueError(
            "Solo pueden bajarse índices contravariantes."
        )

    nuevos_indices = list(tensor.indices)
    nuevos_indices[posicion] = (nombre, "abajo")

    return Tensor(
        nombre=f"{tensor.nombre}↓",
        variedad=tensor.variedad,
        componentes=dict(tensor.componentes),
        indices=nuevos_indices,
  )
