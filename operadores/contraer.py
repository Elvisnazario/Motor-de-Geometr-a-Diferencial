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
def generar_indices_resultantes(
    tensor_a: Tensor,
    tensor_b: Tensor,
    posicion_a: int,
    posicion_b: int,
):
    """
    Genera la lista de índices del tensor resultante después
    de una contracción.

    No realiza ninguna suma de Einstein.
    Únicamente elimina los dos índices contraídos y conserva
    el resto en el orden correcto.

    Ejemplo
    -------
    A^(μ)_(ν)  ×  B^(ν)_(α)

    produce

    C^(μ)_(α)
    """

    indices_resultantes = []

    # Índices del primer tensor
    for indice, estructura in enumerate(tensor_a.indices):
        if indice != posicion_a:
            indices_resultantes.append(estructura)

    # Índices del segundo tensor
    for indice, estructura in enumerate(tensor_b.indices):
        if indice != posicion_b:
            indices_resultantes.append(estructura)

    return indices_resultantes

def contraer(
    tensor_a: Tensor,
    tensor_b: Tensor,
    posicion_a: int,
    posicion_b: int,
):
    """
    Ejecuta la contracción de Einstein entre dos tensores.

    Parámetros
    ----------
    tensor_a : Tensor

    tensor_b : Tensor

    posicion_a : int
        Índice del primer tensor que será contraído.

    posicion_b : int
        Índice del segundo tensor que será contraído.

    Retorna
    -------
    Tensor
        Tensor resultante de la contracción.
    """

    validar_contraccion(
        tensor_a,
        tensor_b,
        posicion_a,
        posicion_b,
    )

    indices_resultantes = generar_indices_resultantes(
        tensor_a,
        tensor_b,
        posicion_a,
        posicion_b,
    )

    componentes_resultantes = defaultdict(lambda: 0)

    dimension = tensor_a.dimension

    for clave_a, valor_a in tensor_a.componentes.items():

        for clave_b, valor_b in tensor_b.componentes.items():

            # Deben coincidir en el índice contraído
            if clave_a[posicion_a] != clave_b[posicion_b]:
                continue

            nueva_clave = []

            # Conservamos todos los índices del tensor A
            for indice, componente in enumerate(clave_a):
                if indice != posicion_a:
                    nueva_clave.append(componente)

            # Conservamos todos los índices del tensor B
            for indice, componente in enumerate(clave_b):
                if indice != posicion_b:
                    nueva_clave.append(componente)

            nueva_clave = tuple(nueva_clave)

            componentes_resultantes[nueva_clave] += (
                valor_a * valor_b
            )

    return Tensor(
        nombre=f"{tensor_a.nombre}_{tensor_b.nombre}",
        componentes=dict(componentes_resultantes),
        indices=indices_resultantes,
        dimension=dimension,
    )
