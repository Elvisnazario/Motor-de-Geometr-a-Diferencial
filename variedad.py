from nucleo.excepciones import ErrorDimension, ErrorObjetoNoRegistrado


class Variedad:
    """
    Variedad diferenciable M.

    Contiene:
    - coordenadas
    - dimensión
    - firma métrica (+--- o -+++)
    - registro de objetos geométricos activos
    """

    def __init__(self, coordenadas, firma):
        self.coordenadas = coordenadas
        self.dimension = len(coordenadas)
        self.firma = firma

        self._validar_firma()

        # Registro de objetos geométricos activos
        self._objetos = {
            "metrica": None,
            "conexion": None,
            "curvatura": None,
        }

    # -------------------------
    # VALIDACIÓN
    # -------------------------
    def _validar_firma(self):
        if len(self.firma) != self.dimension:
            raise ErrorDimension(
                f"Firma {self.firma} incompatible con dimensión {self.dimension}"
            )

    # -------------------------
    # REGISTRO DE OBJETOS
    # -------------------------
    def registrar_objeto(self, tipo, objeto):
        if tipo not in self._objetos:
            raise ErrorObjetoNoRegistrado(
                f"Tipo de objeto no permitido: {tipo}"
            )

        self._objetos[tipo] = objeto

    def obtener_objeto(self, tipo):
        if tipo not in self._objetos:
            raise ErrorObjetoNoRegistrado(f"Tipo desconocido: {tipo}")

        if self._objetos[tipo] is None:
            raise ErrorObjetoNoRegistrado(
                f"No hay objeto registrado en: {tipo}"
            )

        return self._objetos[tipo]

    # -------------------------
    # UTILIDAD
    # -------------------------
    def activar_metrica(self, metrica):
        self.registrar_objeto("metrica", metrica)

    def activar_conexion(self, conexion):
        self.registrar_objeto("conexion", conexion)

    def activar_curvatura(self, curvatura):
        self.registrar_objeto("curvatura", curvatura)

    # -------------------------
    # REPRESENTACIÓN
    # -------------------------
    def __repr__(self):
        return f"Variedad(dim={self.dimension}, coords={self.coordenadas})"
