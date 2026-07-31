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

r = np.linspace(0.1, 40.0, 1000)

# =====================================================================
# 3. COMPONENTES DEL MODELO CONTINUO
# =====================================================================
M_b_reg = M_bar * (r**2) / (r**2 + r_0**2)
x = r / r_c
M_vac = 4.0 * np.pi * rho_0 * (r_c**3) * (x - np.arctan(x))
M_eff = M_b_reg + M_vac

# A(r) derivado del potencial de distribución continua
A = 1.0 - (2.0 * G * M_eff) / ((c**2) * r)

# Gradiente métrico consistente con distribución de masa extendida: dA/dr = 2 G M_eff / (c^2 r^2)
dA_dr = (2.0 * G * M_eff) / ((c**2) * (r**2))

# Velocidad circular de la métrica: v^2 = c^2 r (dA/dr) / (2 A)
v2 = ((c**2) * r * dA_dr) / (2.0 * A)
v_metric = np.sqrt(np.maximum(0.0, v2))
v_newton = np.sqrt(G * M_eff / r)

# =====================================================================
# 4. GRAFICACIÓN DE DIAGNÓSTICO
# =====================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

# Graph 1: Masa Efectiva
axes[0, 0].plot(r, M_eff / 1e10, 'b-', label='M_eff(r)')
axes[0, 0].plot(r, M_b_reg / 1e10, 'k--', label='M_bar(r)')
axes[0, 0].set_ylabel('Masa (10^10 M_sun)')
axes[0, 0].set_title('1. Perfil de Masa Integrada')
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].legend()

# Graph 2: Métrica A(r)
axes[0, 1].plot(r, 1 - A, 'r-')
axes[0, 1].set_ylabel('1 - A(r)')
axes[0, 1].set_title('2. Potencial Métrico')
axes[0, 1].grid(True, alpha=0.3)

# Graph 3: Derivada dA/dr
axes[1, 0].plot(r, dA_dr, 'g-')
axes[1, 0].set_xlabel('Radio r (kpc)')
axes[1, 0].set_ylabel('dA/dr')
axes[1, 0].set_title('3. Gradiente Métrico Continuo dA/dr')
axes[1, 0].grid(True, alpha=0.3)

# Graph 4: Curva de Velocidad
axes[1, 1].plot(r, v_metric, 'r-', linewidth=2, label='Métrica Continua')
axes[1, 1].plot(r, v_newton, 'b--', linewidth=1.5, label='Newtoniano')
axes[1, 1].set_xlabel('Radio r (kpc)')
axes[1, 1].set_ylabel('Velocidad v (km/s)')
axes[1, 1].set_title('4. Curva de Rotación Resultante')
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].legend()

plt.tight_layout()
plt.show()
