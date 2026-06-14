import sympy as sp


class Riemann:
    """
    Tensor de curvatura de Riemann:

        R^ρ_{σμν} =
        ∂_μ Γ^ρ_{νσ}
        - ∂_ν Γ^ρ_{μσ}
        + Γ^ρ_{μλ} Γ^λ_{νσ}
        - Γ^ρ_{νλ} Γ^λ_{μσ}

    Construido a partir de una conexión afín (Christoffel).
    """

    def __init__(self, variedad, conexion):
        self.variedad = variedad
        self.conexion = conexion

        self.coords = variedad.coordenadas
        self.dim = variedad.dimension

        self.Gamma = conexion.Gamma

        self._calcular_riemann()

    # -------------------------
    # CÁLCULO DEL TENSOR DE RIEMANN
    # -------------------------
    def _calcular_riemann(self):
        self.R = sp.MutableDenseNDimArray.zeros(
            self.dim, self.dim, self.dim, self.dim
        )

        for rho in range(self.dim):
            for sigma in range(self.dim):
                for mu in range(self.dim):
                    for nu in range(self.dim):

                        # términos derivados
                        term1 = sp.diff(
                            self.Gamma[rho, nu, sigma],
                            self.coords[mu]
                        )

                        term2 = sp.diff(
                            self.Gamma[rho, mu, sigma],
                            self.coords[nu]
                        )

                        # términos no lineales
                        suma1 = 0
                        suma2 = 0

                        for lam in range(self.dim):
                            suma1 += (
                                self.Gamma[rho, mu, lam]
                                * self.Gamma[lam, nu, sigma]
                            )

                            suma2 += (
                                self.Gamma[rho, nu, lam]
                                * self.Gamma[lam, mu, sigma]
                            )

                        self.R[rho, sigma, mu, nu] = sp.simplify(
                            term1 - term2 + suma1 - suma2
                        )

    # -------------------------
    # ACCESO
    # -------------------------
    def __getitem__(self, clave):
        rho, sigma, mu, nu = clave
        return self.R[rho, sigma, mu, nu]

    # -------------------------
    # PROPIEDAD IMPORTANTE
    # -------------------------
    def es_plano(self):
        """Detecta si la variedad es localmente plana (R = 0)."""
        return self.R.is_zero_matrix

    # -------------------------
    # REPRESENTACIÓN
    # -------------------------
    def __repr__(self):
        return f"Riemann(dim={self.dim})"
