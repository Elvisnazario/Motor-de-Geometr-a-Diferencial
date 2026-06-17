"""
operadores/producto_tensorial.py
================================

Implementación del producto tensorial general del
Motor de Geometría Diferencial (MGD).

El producto tensorial NO realiza contracciones.

Simplemente construye un nuevo tensor cuyo rango es la suma
de los rangos de los dos tensores de entrada.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

from collections import defaultdict

from nucleo.tensor import Tensor
from nucleo.excepciones import (
    ErrorDimension,
)


def producto_tensorial(
    tensor_a,
    tensor_b,
):
    """
    Calcula el producto tensorial abierto.

    Parámetros
    ----------
    tensor_a : Tensor

    tensor_b : Tensor

    Retorna
    -------
    Tensor
    """

    if not isinstance(tensor_a, Tensor):
        raise TypeError(
            "tensor_a debe ser una instancia de Tensor."
        )

    if not isinstance(tensor_b, Tensor):
        raise TypeError(
            "tensor_b debe ser una instancia de Tensor."
        )

    if tensor_a.dimension != tensor_b.dimension:

        raise ErrorDimension(
            "Los tensores pertenecen a variedades "
            "de distinta dimensión."
        )

    if tensor_a.variedad is not tensor_b.variedad:

        raise ErrorDimension(
            "Los tensores no pertenecen a la misma variedad."
        )

    componentes = defaultdict(lambda: 0)

    for clave_a, valor_a in tensor_a.componentes.items():

        for clave_b, valor_b in tensor_b.componentes.items():

            nueva_clave = tuple(clave_a) + tuple(clave_b)

            componentes[nueva_clave] = valor_a * valor_b

    nuevos_indices = (
        list(tensor_a.indices)
        + list(tensor_b.indices)
    )

    return Tensor(

        nombre=(
            f"{tensor_a.nombre}"
            "⊗"
            f"{tensor_b.nombre}"
        ),

        variedad=tensor_a.variedad,

        componentes=dict(componentes),

        indices=nuevos_indices,
    )
