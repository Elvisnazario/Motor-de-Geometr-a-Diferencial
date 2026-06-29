"""
operadores/simetrias.py
=======================

Herramientas para verificar simetrías tensoriales.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from nucleo.tensor import Tensor


def es_simetrico(
    tensor,
    posicion1,
    posicion2,
):
    """
    Comprueba si un tensor es simétrico
    respecto a dos índices.
    """

    if not isinstance(tensor, Tensor):
        raise TypeError(
            "Se esperaba un Tensor."
        )

    for clave, valor in tensor.componentes.items():

        clave_perm = list(clave)

        clave_perm[posicion1], clave_perm[posicion2] = (
            clave_perm[posicion2],
            clave_perm[posicion1],
        )

        clave_perm = tuple(clave_perm)

        otro = tensor.componentes.get(
            clave_perm,
            0,
        )

        if sp.simplify(valor - otro) != 0:
            return False

    return True


def es_antisimetrico(
    tensor,
    posicion1,
    posicion2,
):
    """
    Comprueba si un tensor es antisimétrico.
    """

    if not isinstance(tensor, Tensor):
        raise TypeError(
            "Se esperaba un Tensor."
        )

    for clave, valor in tensor.componentes.items():

        clave_perm = list(clave)

        clave_perm[posicion1], clave_perm[posicion2] = (
            clave_perm[posicion2],
            clave_perm[posicion1],
        )

        clave_perm = tuple(clave_perm)

        otro = tensor.componentes.get(
            clave_perm,
            0,
        )

        if sp.simplify(valor + otro) != 0:
            return False

    return True


def imponer_simetria(
    tensor,
    posicion1,
    posicion2,
):
    """
    Fuerza la simetría de un tensor
    promediando las componentes.

    T' = (T + T^T)/2
    """

    nuevo = tensor.copiar()

    visitados = set()

    for clave, valor in tensor.componentes.items():

        if clave in visitados:
            continue

        clave_perm = list(clave)

        clave_perm[posicion1], clave_perm[posicion2] = (
            clave_perm[posicion2],
            clave_perm[posicion1],
        )

        clave_perm = tuple(clave_perm)

        otro = tensor.componentes.get(
            clave_perm,
            0,
        )

        promedio = sp.simplify(
            (valor + otro) / 2
        )

        nuevo.componentes[clave] = promedio
        nuevo.componentes[clave_perm] = promedio

        visitados.add(clave)
        visitados.add(clave_perm)

    return nuevo


def imponer_antisimetria(
    tensor,
    posicion1,
    posicion2,
):
    """
    Fuerza la antisimetría

    T' = (T - T^T)/2
    """

    nuevo = tensor.copiar()

    visitados = set()

    for clave, valor in tensor.componentes.items():

        if clave in visitados:
            continue

        clave_perm = list(clave)

        clave_perm[posicion1], clave_perm[posicion2] = (
            clave_perm[posicion2],
            clave_perm[posicion1],
        )

        clave_perm = tuple(clave_perm)

        otro = tensor.componentes.get(
            clave_perm,
            0,
        )

        diferencia = sp.simplify(
            (valor - otro) / 2
        )

        nuevo.componentes[clave] = diferencia
        nuevo.componentes[clave_perm] = -diferencia

        visitados.add(clave)
        visitados.add(clave_perm)

    return nuevo