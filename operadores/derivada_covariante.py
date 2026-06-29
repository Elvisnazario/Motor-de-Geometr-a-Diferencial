"""
operadores/derivada_covariante.py
=================================

Derivada covariante completamente general para tensores de cualquier rango.

Proyecto:
    Motor de Geometría Diferencial (MGD)

Autor:
    Elvis Omar Nazario Espinoza
"""

import sympy as sp

from nucleo.tensor import Tensor
from nucleo.excepciones import ErrorDimension


def derivada_covariante(tensor, conexion, coordenada):
    """
    Calcula la derivada covariante completa

        ∇_λ T

    para un tensor de rango arbitrario.

    Cada índice contravariante aporta un término positivo
    con Christoffel.

    Cada índice covariante aporta un término negativo.

    Parámetros
    ----------
    tensor : Tensor

    conexion : Conexion

    coordenada : int
    """

    if not isinstance(tensor, Tensor):
        raise TypeError(
            "El primer argumento debe ser un Tensor."
        )

    if coordenada < 0 or coordenada >= tensor.dimension:
        raise ErrorDimension(
            "Coordenada fuera del rango."
        )

    dim = tensor.dimension
    x = tensor.variedad.coordenadas[coordenada]

    nuevas_componentes = {}

    for clave, valor in tensor.componentes.items():

        resultado = sp.diff(valor, x)

        # --------------------------------------
        # Correcciones de Christoffel
        # --------------------------------------

        for posicion, (_, tipo) in enumerate(tensor.indices):

            indice_original = clave[posicion]

            if tipo == "arriba":

                for alpha in range(dim):

                    clave_aux = list(clave)
                    clave_aux[posicion] = alpha

                    resultado += (
                        conexion.Gamma[
                            indice_original,
                            coordenada,
                            alpha,
                        ]
                        * tensor[tuple(clave_aux)]
                    )

            else:

                for alpha in range(dim):

                    clave_aux = list(clave)
                    clave_aux[posicion] = alpha

                    resultado -= (
                        conexion.Gamma[
                            alpha,
                            coordenada,
                            indice_original,
                        ]
                        * tensor[tuple(clave_aux)]
                    )

        resultado = sp.simplify(resultado)

        if resultado != 0:
            nuevas_componentes[clave] = resultado

    nuevos_indices = list(tensor.indices)
    nuevos_indices.append(
        ("λ", "abajo")
    )

    return Tensor(
        nombre=f"∇({tensor.nombre})",
        variedad=tensor.variedad,
        componentes=nuevas_componentes,
        indices=nuevos_indices,
    )


def derivada_parcial(tensor, coordenada):
    """
    Derivada parcial de un tensor.

    No utiliza la conexión.
    """

    if coordenada < 0 or coordenada >= tensor.dimension:
        raise ErrorDimension(
            "Coordenada fuera del rango."
        )

    simbolo = tensor.variedad.coordenadas[coordenada]

    componentes = {}

    for clave, valor in tensor.componentes.items():

        nuevo = sp.diff(valor, simbolo)

        if nuevo != 0:
            componentes[clave] = sp.simplify(nuevo)

    nuevos_indices = list(tensor.indices)
    nuevos_indices.append(
        ("λ", "abajo")
    )

    return Tensor(
        nombre=f"∂({tensor.nombre})",
        variedad=tensor.variedad,
        componentes=componentes,
        indices=nuevos_indices,
    )