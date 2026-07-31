import numpy as np
import matplotlib.pyplot as plt

# =====================================================================
# 1. CONSTANTES FÍSICAS
# =====================================================================
G = 4.30091e-6          # (km/s)^2 * kpc / M_sun
c = 299792.458          # km/s

# =====================================================================
# 2. PARÁMETROS FÍSICOS FIJOS (Galaxia genérica)
# =====================================================================
M_bar = 5e10            # Masa bariónica (M_sun)
r_0   = 3.0             # Radio de escala bariónico (kpc)
rho_0 = 1e7             # Densidad del vacío central (M_sun / kpc^3)
r_c   = 12.0            # Radio de escala del núcleo elástico (kpc)

r = np.linspace(0.1, 40.0, 1000) # Radio de 0.1 a 40 kpc

# =====================================================================
# 3. COMPONENTES DEL MODELO
# =====================================================================
# Masa bariónica regulada
M_b_reg = M_bar * (r**2) / (r**2 + r_0**2)

# Masa del vacío integrada: rho(r) = rho_0 / (1 + (r/r_c)^2)
x = r / r_c
M_vac = 4.0 * np.pi * rho_0 * (r_c**3) * (x - np.arctan(x))

# Masa efectiva total
M_eff = M_b_reg + M_vac

# Componente métrica A(r)
A = 1.0 - (2.0 * G * M_eff) / ((c**2) * r)

# Derivada numérica limpia dA/dr (Diferencias finitas centrales)
dr = r[1] - r[0]
dA_dr = np.gradient(A, dr)

# Velocidad circular derivada estrictamente de la métrica:
# v^2 = (c^2 * r * A') / (2 * A)
v2 = ((c**2) * r * dA_dr) / (2.0 * A)

# Limpieza de valores no físicos (si v^2 < 0)
v_metric = np.sqrt(np.maximum(0.0, v2))

# Velocidad newtoniana estándar de control (v^2 = G * M_eff / r)
v_newton = np.sqrt(G * M_eff / r)

# =====================================================================
# 4. GRAFICACIÓN DE DIAGNÓSTICO
# =====================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

# Graph 1: Masa Efectiva
axes[0, 0].plot(r, M_eff / 1e10, 'b-', label=r'$M_{eff}(r)$')
axes[0, 0].plot(r, M_b_reg / 1e10, 'k--', label=r'$M_{bar}(r)$')
axes[0, 0].set_ylabel(r'Masa ($10^{10} M_{\odot}$)')
axes[0, 0].set_title('1. Perfil de Masa Integrada')
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].legend()

# Graph 2: Métrica A(r)
axes[0, 1].plot(r, 1 - A, 'r-')
axes[0, 1].set_ylabel(r'$1 - A(r) = \frac{2GM}{c^2r}$')
axes[0, 1].set_title('2. Potencial Métrico (Desviación de 1)')
axes[0, 1].grid(True, alpha=0.3)

# Graph 3: Derivada dA/dr
axes[1, 0].plot(r, dA_dr, 'g-')
axes[1, 0].set_xlabel('Radio r (kpc)')
axes[1, 0].set_ylabel(r'$A\'(r)$')
axes[1, 0].set_title(r'3. Gradiente Métrico $A\'(r)$')
axes[1, 0].grid(True, alpha=0.3)

# Graph 4: Curva de Velocidad
axes[1, 1].plot(r, v_metric, 'r-', linewidth=2, label=r'Métrica: $v^2 = \frac{c^2 r A\'}{2 A}$')
axes[1, 1].plot(r, v_newton, 'b--', linewidth=1.5, label=r'Newton: $v^2 = \frac{G M}{r}$')
axes[1, 1].set_xlabel('Radio r (kpc)')
axes[1, 1].set_ylabel('Velocidad v (km/s)')
axes[1, 1].set_title('4. Curva de Rotación Resultante')
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].legend()

plt.tight_layout()
plt.show()
