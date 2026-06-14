import sympy as sp
from nucleo.excepciones import ErrorDimension


class Metrica:
    """
    Métrica disforme:

    g_{μν} = A(ψ, ω) η_{μν} + B(ψ, ω) ψ_μ ψ_ν
    """

    def __init__(self, variedad, A, B, psi, eta):
        self.variedad = variedad
        self.A = A
        self.B = B
        self.psi = psi
        self.eta = eta

        self.dimension = variedad.dimension

        self._construir_metrica_cov()

    # -------------------------
    # NORMA ω = ψ^μ ψ_μ
    # -------------------------
    def _omega(self, psi_sup, psi_cov):
        return sum(psi_sup[i] * psi_cov[i] for i in range(self.dimension))

    # -------------------------
    # SUBIR ÍNDICE (Minkowski)
    # -------------------------
    def _subir_indice(self, vector):
        return [
            sum(self.eta[i, j] * vector[j] for j in range(self.dimension))
            for i in range(self.dimension)
        ]

    # -------------------------
    # CONSTRUCCIÓN g_{μν}
    # -------------------------
    def _construir_metrica_cov(self):
        psi_cov = self.psi
        psi_sup = self._subir_indice(psi_cov)

        omega = self._omega(psi_sup, psi_cov)

        self.g_cov = sp.Matrix.zeros(self.dimension, self.dimension)

        for mu in range(self.dimension):
            for nu in range(self.dimension):
                self.g_cov[mu, nu] = (
                    self.A * self.eta[mu, nu]
                    + self.B * psi_cov[mu] * psi_cov[nu]
                )

        self._psi_sup = psi_sup
        self._omega_val = omega

    # -------------------------
    # INVERSA SHERMAN-MORRISON
    # -------------------------
    def inversa(self):
        psi_sup = self._psi_sup
        omega = self._omega_val

        Delta = self.A + self.B * omega

        g_inv = sp.Matrix.zeros(self.dimension, self.dimension)

        for mu in range(self.dimension):
            for nu in range(self.dimension):
                g_inv[mu, nu] = (
                    (1 / self.A) * self.eta[mu, nu]
                    - (self.B / (self.A * Delta)) * psi_sup[mu] * psi_sup[nu]
                )

        return g_inv

    # -------------------------
    # SUBIR ÍNDICES DE UN VECTOR
    # -------------------------
    def subir_vector(self, vector):
        g_inv = self.inversa()

        return [
            sum(g_inv[i, j] * vector[j] for j in range(self.dimension))
            for i in range(self.dimension)
        ]

    # -------------------------
    # VALIDACIÓN
    # -------------------------
    def validar(self):
        if self.g_cov.shape != (self.dimension, self.dimension):
            raise ErrorDimension("Métrica mal formada")

    # -------------------------
    # REPRESENTACIÓN
    # -------------------------
    def __repr__(self):
        return f"Metrica(dim={self.dimension})"
