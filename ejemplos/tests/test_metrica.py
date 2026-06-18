import sympy as sp


def comprobar_simetria(tensor):
    """
    Verifica que g_{μν} = g_{νμ}.
    """
    for (i, j), valor in tensor.componentes.items():
        if sp.simplify(valor - tensor.obtener_componente((j, i))) != 0:
            return False
    return True


def determinante_metrica(tensor):
    """
    Calcula el determinante simbólico de la métrica.
    """
    dim = tensor.variedad.dimension

    matriz = sp.Matrix(
        dim,
        dim,
        lambda i, j: tensor.obtener_componente((i, j))
    )

    return sp.factor(matriz.det())
