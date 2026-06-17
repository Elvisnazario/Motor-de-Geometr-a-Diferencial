"""
operadores/sustituciones.py
===========================

Herramientas oficiales para realizar sustituciones
simbólicas sobre expresiones y tensores del
Motor de Geometría Diferencial (MGD).

Proyecto:
    Motor de Geometría Diferencial
"""

from copy import deepcopy

from nucleo.tensor import Tensor


def sustituir(expresion, reglas):
    """
    Aplica sustituciones simbólicas sobre una expresión.

    Parámetros
    ----------
    expresion

    reglas : dict

    Retorna
    -------
    Expresión sustituida.
    """

    return expresion.subs(reglas)


def sustituir_tensor(
    tensor,
    reglas,
):
    """
    Devuelve una copia del tensor con todas sus
    componentes sustituidas.

    El tensor original permanece intacto.
    """

    if not isinstance(tensor, Tensor):
        raise TypeError(
            "Debe proporcionarse un Tensor."
        )

    nuevo = deepcopy(tensor)

    for clave in nuevo.componentes:

        nuevo.componentes[clave] = (
            nuevo.componentes[clave].subs(reglas)
        )

    return nuevo


def sustituir_lista(
    lista,
    reglas,
):
    """
    Sustituye todos los elementos de una lista.
    """

    return [
        elemento.subs(reglas)
        for elemento in lista
    ]


def sustituir_diccionario(
    diccionario,
    reglas,
):
    """
    Sustituye todos los valores de un diccionario.
    """

    nuevo = {}

    for clave, valor in diccionario.items():

        nuevo[clave] = valor.subs(reglas)

    return nuevo
