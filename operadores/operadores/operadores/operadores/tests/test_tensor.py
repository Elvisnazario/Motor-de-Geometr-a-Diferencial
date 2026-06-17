"""
tests/test_tensor.py
====================

Pruebas unitarias del objeto Tensor.

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

import sympy as sp

from nucleo.variedad import Variedad
from nucleo.tensor import Tensor


def construir_tensor_basico():

    t, x, y, z = sp.symbols(
        "t x y z"
    )

    variedad = Variedad(

        coordenadas=(t, x, y, z),

        firma=(1, -1, -1, -1),
    )

    componentes = {

        (0, 0): sp.Integer(1),

        (1, 1): sp.Integer(2),
    }

    indices = [

        ("μ", "arriba"),

        ("ν", "abajo"),
    ]

    return Tensor(

        nombre="T",

        variedad=variedad,

        componentes=componentes,

        indices=indices,
    )


def test_rango():

    tensor = construir_tensor_basico()

    assert tensor.rango == 2


def test_dimension():

    tensor = construir_tensor_basico()

    assert tensor.dimension == 4


def test_componente_existente():

    tensor = construir_tensor_basico()

    assert tensor[(0, 0)] == 1


def test_componente_inexistente():

    tensor = construir_tensor_basico()

    assert tensor[(3, 2)] == 0


def test_copia():

    tensor = construir_tensor_basico()

    copia = tensor.copiar()

    assert copia == tensor

    assert copia is not tensor


def test_suma():

    A = construir_tensor_basico()

    B = construir_tensor_basico()

    C = A + B

    assert C[(0, 0)] == 2

    assert C[(1, 1)] == 4


def test_resta():

    A = construir_tensor_basico()

    B = construir_tensor_basico()

    C = A - B

    assert C[(0, 0)] == 0

    assert C[(1, 1)] == 0


def test_escalar():

    tensor = construir_tensor_basico()

    doble = 2 * tensor

    assert doble[(0, 0)] == 2

    assert doble[(1, 1)] == 4


def test_repr():

    tensor = construir_tensor_basico()

    texto = repr(tensor)

    assert "Tensor" in texto

    assert "T" in texto
