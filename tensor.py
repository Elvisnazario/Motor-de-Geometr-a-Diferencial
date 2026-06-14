from collections import defaultdict
from nucleo.excepciones import ErrorDimension, ErrorIndicesIncompatibles


class Tensor:
    """
    Representación básica de un tensor en el Motor de Geometría Diferencial (MGD).

    Un tensor se define por:
    - componentes: diccionario disperso de SymPy
    - indices: lista de tuplas (nombre, posicion)
    - dimension: dimensión de la variedad
    """

    def __init__(self, nombre, componentes, indices, dimension):
        self.nombre = nombre
        self.componentes = defaultdict(lambda: 0, componentes)
        self.indices = indices
        self.dimension = dimension

        self._validar()

    # -------------------------
    # VALIDACIÓN
    # -------------------------
    def _validar(self):
        # validar índices
        for idx in self.indices:
            if not isinstance(idx, tuple) or len(idx) != 2:
                raise ErrorDimension(f"Índice inválido en {self.nombre}: {idx}")

            nombre_idx, posicion = idx

            if posicion not in ("arriba", "abajo"):
                raise ErrorDimension(
                    f"Posición inválida en {self.nombre}: {posicion}"
                )

        # validar componentes
        for clave in self.componentes.keys():
            if len(clave) != len(self.indices):
                raise ErrorDimension(
                    f"Componente {clave} no coincide con el rango del tensor {self.nombre}"
                )

    # -------------------------
    # ACCESO
    # -------------------------
    def __getitem__(self, clave):
        return self.componentes[clave]

    def __setitem__(self, clave, valor):
        if len(clave) != len(self.indices):
            raise ErrorDimension("Clave incompatible con el rango del tensor")
        self.componentes[clave] = valor

    # -------------------------
    # OPERACIONES BÁSICAS
    # -------------------------
    def __add__(self, otro):
        if self.indices != otro.indices:
            raise ErrorIndicesIncompatibles("No se pueden sumar tensores incompatibles")

        resultado = Tensor(
            nombre=f"({self.nombre}+{otro.nombre})",
            componentes={},
            indices=self.indices,
            dimension=self.dimension,
        )

        for k in set(self.componentes.keys()).union(otro.componentes.keys()):
            resultado.componentes[k] = self[k] + otro[k]

        return resultado

    def __sub__(self, otro):
        if self.indices != otro.indices:
            raise ErrorIndicesIncompatibles("No se pueden restar tensores incompatibles")

        resultado = Tensor(
            nombre=f"({self.nombre}-{otro.nombre})",
            componentes={},
            indices=self.indices,
            dimension=self.dimension,
        )

        for k in set(self.componentes.keys()).union(otro.componentes.keys()):
            resultado.componentes[k] = self[k] - otro[k]

        return resultado

    # -------------------------
    # MULTIPLICACIÓN POR ESCALAR
    # -------------------------
    def __mul__(self, escalar):
        resultado = Tensor(
            nombre=f"{escalar}*{self.nombre}",
            componentes={},
            indices=self.indices,
            dimension=self.dimension,
        )

        for k in self.componentes:
            resultado.componentes[k] = escalar * self.componentes[k]

        return resultado

    __rmul__ = __mul__

    # -------------------------
    # REPRESENTACIÓN
    # -------------------------
    def __repr__(self):
        return f"Tensor({self.nombre}, rango={len(self.indices)}, dim={self.dimension})"
