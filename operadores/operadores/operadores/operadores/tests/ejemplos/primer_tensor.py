"""
ejemplos/primer_tensor.py
=========================

Primer ejemplo de utilización del Motor de Geometría
Diferencial (MGD).

Proyecto:
    Motor de Geometría Diferencial (MGD)
"""

import sympy as sp

from nucleo.variedad import Variedad
from nucleo.tensor import Tensor


# --------------------------------------------------
# Coordenadas
# --------------------------------------------------

t, x, y, z = sp.symbols("t x y z")

# --------------------------------------------------
# Variedad
# --------------------------------------------------

M = Variedad(

    coordenadas=(t, x, y, z),

    firma=(1, -1, -1, -1),
)

# --------------------------------------------------
# Componentes del tensor
# --------------------------------------------------

componentes = {

    (0, 0): 1,

    (1, 1): x**2,

    (2, 2): y**2,

    (3, 3): z**2,
}

# --------------------------------------------------
# Índices
# --------------------------------------------------

indices = [

    ("μ", "abajo"),

    ("ν", "abajo"),
]

# --------------------------------------------------
# Construcción
# --------------------------------------------------

g = Tensor(

    nombre="g",

    variedad=M,

    componentes=componentes,

    indices=indices,
)

print()

print("Motor de Geometría Diferencial")

print("--------------------------------")

print()

print(g)

print()

print("Componentes no nulas:")

for clave, valor in g.componentes.items():

    print(clave, "=", valor)

print()

print("Dimensión :", g.dimension)

print("Rango     :", g.rango)

print("Índices   :", g.indices)
