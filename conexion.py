import sympy as sp
from nucleo.excepciones import ErrorDimension


class Conexion:
    """
    Conexión afín (tipo Levi-Civita o generalizada).

    NO es un tensor.
    Representa Γ^α_{μν} como objeto geométrico dependiente de una variedad.
    """

    def __init__(self, variedad, metrica):
        self.variedad = variedad
        self.metrica = metrica

        self.coords = variedad.coordenadas
        self.dim = variedad.dimension

        self.g_cov = metrica.g_cov
        self.g_inv = metrica.inversa()

        self._calcular_derivadas_metrica()
        self._calcular_christoffel()

    # -------------------------
    # DERIVADAS DE LA MÉTRICA
    # -------------------------
    def _calcular_derivadas_metrica(self):
        self.dg = sp.MutableDenseNDimArray.zeros(
            self.dim, self.dim, self.dim
        )

        for lam in range(self.dim):
            for mu in range(self.dim):
                for nu in range(self.dim):
                    self.dg[lam, mu, nu] = sp.diff(
                        self.g_cov[mu, nu],
                        self.coords[lam]
                    )

    # -------------------------
    # CÁLCULO DE CHRISTOFFEL
    # -------------------------
    def _calcular_christoffel(self):
        self.Gamma = sp.MutableDenseNDimArray.zeros(
            self.dim, self.dim, self.dim
        )

        for alpha in range(self.dim):
            for mu in range(self.dim):
                for nu in range(self.dim):

                    suma = 0

                    for sigma in range(self.dim):
                        suma += (
                            0.5
                            * self.g_inv[alpha, sigma]
                            * (
                                self.dg[mu, nu, sigma]
                                + self.dg[nu, mu, sigma]
                                - self.dg[sigma, mu, nu]
                            )
                        )

                    self.Gamma[alpha, mu, nu] = sp.simplify(suma)

    # -------------------------
    # ACCESO
    # -------------------------
    def __getitem__(self, clave):
        alpha, mu, nu = clave
        return self.Gamma[alpha, mu, nu]

    # -------------------------
    # SIMETRÍA (TORSIÓN)
    # -------------------------
    def es_sin_torsion(self):
        for a in range(self.dim):
            for m in range(self.dim):
                for n in range(self.dim):
                    if sp.simplify(
                        self.Gamma[a, m, n] - self.Gamma[a, n, m]
                    ) != 0:
                        return False
        return True

    # -------------------------
    # REPRESENTACIÓN
    # -------------------------
    def __repr__(self):
        return f"Conexion(dim={self.dim}, tipo=LeviCivita?)"
