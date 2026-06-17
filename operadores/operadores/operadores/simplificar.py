"""
operadores/simplificar.py
=========================

Rutinas oficiales de simplificación simbólica del
Motor de Geometría Diferencial (MGD).

Toda simplificación matemática del motor debe realizarse
a través de este módulo para garantizar un comportamiento
uniforme en toda la biblioteca.

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

import sympy as sp


def simplificar(expresion):
    """
    Ejecuta la cascada oficial de simplificación del MGD.

    Orden utilizado:

        expand
        cancel
        factor
        simplify

    Parámetros
    ----------
    expresion

    Retorna
    -------
    Expresión simplificada.
    """

    if expresion == 0:
        return sp.Integer(0)

    return sp.simplify(
        sp.factor(
            sp.cancel(
                sp.expand(expresion)
            )
        )
    )


def simplificar_tensor(tensor):
    """
    Simplifica todas las componentes de un tensor.

    La operación modifica el tensor recibido y además
    lo devuelve para permitir encadenamiento.
    """

    for clave in list(tensor.componentes.keys()):
        tensor.componentes[clave] = simplificar(
            tensor.componentes[clave]
        )

    return tensor


def simplificar_lista(lista):
    """
    Simplifica cada elemento de una lista.
    """

    return [
        simplificar(x)
        for x in lista
    ]


def simplificar_diccionario(diccionario):
    """
    Simplifica todos los valores de un diccionario.
    """

    return {
        clave: simplificar(valor)
        for clave, valor in diccionario.items()
  }
