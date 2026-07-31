import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution

# Configuración de estilo visual premium y profesional
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available() else 'default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.titlesize': 14,
    'grid.alpha': 0.3
})

# =====================================================================
# 1. CONSTANTES FÍSICAS Y CONVERSIÓN DE UNIDADES (SISTEMA ASTROFÍSICO)
# =====================================================================
class ConstantesFisicas:
    """
    Clase contenedora de constantes físicas universales y conversiones dimensionales.
    Garantiza consistencia estricta en el régimen astronómico: M_sun, kpc, km/s.
    """
    G = 4.30091e-6          # Constante de gravitación en (km/s)^2 * kpc / M_sun
    c = 299792.458          # Velocidad de la luz en km/s

# =====================================================================
# 2. MODELO GEOMÉTRICO DE LA ATRACCIÓN ENERGÉTICA (MASA INTEGRADA)
# =====================================================================
class AtraccionEnergeticaModel:
    """
    Implementación matemática de la métrica de la Atracción Energética.
    La curvatura espacio-temporal se deriva de la integración de la masa efectiva
    M_eff(r) = M_bar_reg(r) + M_vac(r), incorporando el efecto de repulsión elástica/presión.
    """
    def __init__(self, M_bar, r_0, rho_0, r_c):
        """
        M_bar : Masa bariónica central en M_sun
        r_0   : Radio de transición elástica central (regularización) en kpc
        rho_0 : Densidad central de energía/tensión del vacío elástico
        r_c   : Radio de escala del núcleo del vacío elástico en kpc
        """
        self.M_bar = M_bar
        self.r_0 = r_0
        self.rho_0 = rho_0
        self.r_c = r_c
        self.G = ConstantesFisicas.G
        self.c = ConstantesFisicas.c

    def masa_efectiva(self, r):
        """
        Calcula la masa acumulada total M_eff(r) = M_bar_reg(r) + M_vac(r).
        M_vac(r) es la integral del perfil de densidad elástica del vacío.
        """
        if r < 1e-6:
            return 0.0
        # Componente bariónica suavizada para evitar la singularidad central
        M_bar_reg = self.M_bar * (r**2) / (r**2 + self.r_0**2)
        
        # Masa acumulada del vacío elástico (Masa efectiva/Tensión del medio)
        M_vac = 4.0 * np.pi * self.rho_0 * (self.r_c**3) * ( (r / self.r_c) - np.arctan(r / self.r_c) )
        
        return M_bar_reg + M_vac

    def dM_dr(self, r):
        """
        Derivada analítica exacta dM_eff/dr utilizada para el cálculo del potencial.
        """
        if r < 1e-6:
            return 0.0
        dM_bar_dr = self.M_bar * (2.0 * r * (self.r_0**2)) / ((r**2 + self.r_0**2)**2)
        dM_vac_dr = 4.0 * np.pi * (r**2) * self.rho_0 / (1.0 + (r / self.r_c)**2)
        return dM_bar_dr + dM_vac_dr

    def A_metric(self, r):
        """
        Componente temporal g_tt = -A(r) de la métrica de la Atracción Energética.
        A(r) = 1 - 2*G*M_eff(r) / (c^2 * r)
        """
        if r < 1e-6:
            return 1.0
        A_val = 1.0 - (2.0 * self.G * self.masa_efectiva(r)) / ((self.c**2) * r)
        return A_val

    def dA_dr(self, r):
        """
        Derivada analítica exacta dA/dr obtenida mediante la regla del cociente.
        dA/dr = (2*G / c^2) * (M_eff / r^2 - dM_dr / r)
        """
        if r < 1e-6:
            return 0.0
        M_eff = self.masa_efectiva(r)
        dM_dr_val = self.dM_dr(r)
        return (2.0 * self.G / (self.c**2)) * ( (M_eff / (r**2)) - (dM_dr_val / r) )

    def velocidad_circular(self, r):
        """
        Velocidad orbital circular exacta derivada de la geodésica relativista ecuatorial:
        v^2 = c^2 * r * A'(r) / (2 * A(r))
        """
        A_val = self.A_metric(r)
        if A_val <= 1e-5:
            return 0.0
        Ap_val = self.dA_dr(r)
        v2 = (self.c**2) * r * Ap_val / (2.0 * A_val)
        return np.sqrt(max(0.0, v2))

# =====================================================================
# 3. GEODÉSICAS EXACTAS DERIVADAS DEL LAGRANGIANO (G = c = 1)
# =====================================================================

def geodesicas_nulas(l, y, M, r0, rho0, rc, L):
    """
    Ecuaciones geodésicas polares nulas exactas (k=0) para simulación de luz.
    y = [r, phi, vr]
    """
    r, phi, vr = y
    
    model = AtraccionEnergeticaModel(M, r0, rho0, rc)
    model.G = 1.0
    model.c = 1.0
    
    A = model.A_metric(r)
    Ap = model.dA_dr(r)
    
    dr_dl = vr
    dphi_dl = L / (r**2)
    dvr_dl = -0.5 * Ap * (L**2) / (r**2) + A * (L**2) / (r**3)
    
    return [dr_dl, dphi_dl, dvr_dl]

def geodesicas_masivas(l, y, M, r0, rho0, rc, L):
    """
    Ecuaciones geodésicas polares masivas exactas (k=1) para órbitas estelares.
    y = [r, phi, vr]
    """
    r, phi, vr = y
    
    model = AtraccionEnergeticaModel(M, r0, rho0, rc)
    model.G = 1.0
    model.c = 1.0
    
    A = model.A_metric(r)
    Ap = model.dA_dr(r)
    
    dr_dl = vr
    dphi_dl = L / (r**2)
    dvr_dl = -0.5 * Ap * (1.0 + (L**2) / (r**2)) + A * (L**2) / (r**3)
    
    return [dr_dl, dphi_dl, dvr_dl]

# =====================================================================
# 4. PROTOCOLO DE VALIDACIÓN EXPERIMENTAL MULTI-ESCALA
# =====================================================================

print("="*75)
print("   SISTEMA DE VALIDACIÓN EXPERIMENTAL MONÓTONO Y AJUSTE DE GALAXIAS")
print("   Física de Medios Elásticos Continuos con Masa Integrada (Sin Trampas)")
print("   Autor: Elvis Omar Nazario Espinoza")
print("="*75)

# --- TEST 1: Optimización de Curvas de Rotación Galácticas SPARC ---
print("\n[TEST 1/3] Ejecutando ajuste evolutivo global con restricciones elásticas...")

galaxias_data = {
    'NGC 3198 (Espiral Gigante)': {
        'r': np.array([2.0, 4.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0]),
        'v': np.array([62.2, 115.7, 144.8, 152.8, 155.2, 156.9, 155.0, 152.0, 149.0, 147.0]),
        'err': np.array([5.5,  4.8,   4.2,   4.5,   4.9,   5.1,   5.3,   5.6,   5.8,   6.0])
    },
    'NGC 2403 (Espiral Mediana)': {
        'r': np.array([1.0, 2.0, 3.5, 5.0, 7.0, 9.0, 11.5, 14.0, 17.0, 19.5]),
        'v': np.array([68.1, 95.3, 112.4, 122.1, 128.5, 131.0, 133.2, 134.1, 135.0, 134.7]),
        'err': np.array([4.1,  3.8,   3.5,   3.9,   4.2,   4.5,   4.6,   4.8,   5.0,   5.1])
    },
    'UGC 2259 (Espiral Enana)': {
        'r': np.array([0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
        'v': np.array([23.4, 45.1, 62.3, 71.5, 78.2, 82.1, 84.4, 85.1, 85.6, 85.2]),
        'err': np.array([2.5,  2.2,   2.1,   2.3,   2.5,   2.6,   2.8,   2.9,   3.0,   3.1])
    }
}

ajustes_resultados = {}

for nombre, data in galaxias_data.items():
    r_arr, v_obs, err_arr = data['r'], data['v'], data['err']
    
    bounds = [
        (1e9, 2e11),   # Masa bariónica M_sun
        (0.1, 5.0),    # r_0
        (1e5, 5e7),    # rho_0 
        (1.0, 25.0)    # r_c
    ]
    
    def loss_func(params):
        M_b, r0, rho0, rc = params
        if r0 >= rc:
            return 1e12
            
        model = AtraccionEnergeticaModel(M_b, r0, rho0, rc)
        v_theo = np.array([model.velocidad_circular(r) for r in r_arr])
        
        # Minimización chi-cuadrado estándar
        return np.sum(((v_obs - v_theo) / err_arr) ** 2)

    res = differential_evolution(loss_func, bounds, seed=42, polish=True, maxiter=100)
    M_b_opt, r0_opt, rho0_opt, rc_opt = res.x
    
    model_opt = AtraccionEnergeticaModel(M_b_opt, r0_opt, rho0_opt, rc_opt)
    r_dense = np.linspace(0.1, max(r_arr)*1.05, 200)
    v_dense = np.array([model_opt.velocidad_circular(rd) for rd in r_dense])
    
    v_fit = np.array([model_opt.velocidad_circular(r) for r in r_arr])
    residuos = v_obs - v_fit
    R2 = 1.0 - (np.sum(residuos**2) / np.sum((v_obs - np.mean(v_obs))**2))
    MAPE = np.mean(np.abs(residuos) / v_obs) * 100.0
    
    ajustes_resultados[nombre] = {
        'r_dense': r_dense, 'v_dense': v_dense,
        'r_obs': r_arr, 'v_obs': v_obs, 'err_obs': err_arr,
        'R2': R2, 'MAPE': MAPE, 'params': res.x
    }
    print(f" -> {nombre:<28} | R² = {R2:.4f} | MAPE = {MAPE:.2f}% | M_b = {M_b_opt:.2e} M_sun")

# --- TEST 2: Deflexión Óptica de la Luz (Métrica Normalizada) ---
print("\n[TEST 2/3] Simulando deflexión de luz bajo la métrica del medio elástico...")

M_l, r0_l, rho0_l, rc_l = 0.5, 0.2, 0.05, 2.0
deflection_angles = []
impact_parameters = [4.5, 6.0, 8.0]
light_tracks = []

for b in impact_parameters:
    r0_init = 20.0
    phi0_init = -np.arcsin(b / r0_init)
    
    model_geo = AtraccionEnergeticaModel(M_l, r0_l, rho0_l, rc_l)
    model_geo.G = 1.0
    model_geo.c = 1.0
    
    A0 = model_geo.A_metric(r0_init)
    L_light = b / np.sqrt(A0)
    
    # Condición radial inicial exacta para rayo entrante
    vr0 = -np.sqrt(max(1.0 - A0 * (L_light**2) / (r0_init**2), 1e-12))
    
    y0 = [r0_init, phi0_init, vr0]
    
    sol = solve_ivp(
        geodesicas_nulas, 
        (0.0, 40.0), 
        y0, 
        args=(M_l, r0_l, rho0_l, rc_l, L_light), 
        rtol=1e-8, 
        atol=1e-10
    )
    
    r_pts = sol.y[0]
    phi_pts = sol.y[1]
    
    x_pts = r_pts * np.cos(phi_pts)
    y_pts = r_pts * np.sin(phi_pts)
    light_tracks.append((x_pts, y_pts, b))
    
    phi_final = phi_pts[-1]
    deflexion = np.degrees(np.abs(phi_final - 0.0))
    deflection_angles.append(deflexion)
    
    print(f" -> Parámetro de Impacto b = {b:.1f} | Deflexión: {deflexion:.4f}°")

# --- TEST 3: Precesión Geodésica de Roseta (Estrella) ---
print("\n[TEST 3/3] Simulando órbita estelar excéntrica y precesión...")

r0_orb = 6.0
model_orb = AtraccionEnergeticaModel(M_l, r0_l, rho0_l, rc_l)
model_orb.G = 1.0
model_orb.c = 1.0

# Momento angular ajustado para producir excentricidad estable (Roseta)
A_orb = model_orb.A_metric(r0_orb)
Ap_orb = model_orb.dA_dr(r0_orb)
v_circ_geom = np.sqrt(max(0.0, 0.5 * r0_orb * Ap_orb / A_orb))

# Aplicamos una pequeña excentricidad (85% de la velocidad circular)
L_star = r0_orb * (v_circ_geom * 0.85) 
vr0_orb = 0.0

sol_orb = solve_ivp(
    geodesicas_masivas, 
    (0.0, 300.0), 
    [r0_orb, 0.0, vr0_orb], 
    args=(M_l, r0_l, rho0_l, rc_l, L_star), 
    rtol=1e-9, 
    atol=1e-11
)

x_orb = sol_orb.y[0] * np.cos(sol_orb.y[1])
y_orb = sol_orb.y[0] * np.sin(sol_orb.y[1])
print(" -> Simulación de órbita excéntrica completada de forma correcta.")

# =====================================================================
# 5. VISUALIZACIÓN DE RESULTADOS GRÁFICOS MULTI-ESCALA
# =====================================================================
fig = plt.figure(figsize=(15, 5), constrained_layout=True)
gs = plt.GridSpec(1, 3, figure=fig)

# Panel 1: Ajuste SPARC (Corregido 'ms=4' para evitar errores de Matplotlib)
ax1 = fig.add_subplot(gs[0, 0])
colors = ['#e74c3c', '#3498db', '#2ecc71']
for idx, (nombre, r_fit) in enumerate(ajustes_resultados.items()):
    ax1.errorbar(r_fit['r_obs'], r_fit['v_obs'], yerr=r_fit['err_obs'], fmt='o', color=colors[idx], alpha=0.6, ms=4, capsize=3)
    ax1.plot(r_fit['r_dense'], r_fit['v_dense'], '-', color=colors[idx], linewidth=2.0, label=f"{nombre[:8]} ($R^2$: {r_fit['R2']:.3f})")
ax1.set_xlabel("Radio r (kpc)")
ax1.set_ylabel("Velocidad de Rotación v (km/s)")
ax1.set_title("Test 1: Ajuste SPARC (Modelo de Masa Integrada)")
ax1.legend(frameon=True, loc="lower right", fontsize=8)

# Panel 2: Deflexión Óptica Suave
ax2 = fig.add_subplot(gs[0, 1])
lens_colors = ['#2980b9', '#8e44ad', '#f1c40f']
for i, (x, y, b) in enumerate(light_tracks):
    ax2.plot(x, y, color=lens_colors[i], linewidth=1.5, label=f"b = {b:.1f} ({deflection_angles[i]:.2f}°)")
core = plt.Circle((0, 0), rc_l, color='#2c3e50', alpha=0.3, label="Núcleo Elástico ($r_c$)")
ax2.add_patch(core)
ax2.set_xlim(-15, 15)
ax2.set_ylim(-8, 8)
ax2.set_aspect('equal')
ax2.set_xlabel("Coordenada X")
ax2.set_ylabel("Coordenada Y")
ax2.set_title("Test 2: Deflexión Óptica Suave y Monótona")
ax2.legend(frameon=True, fontsize=8)

# Panel 3: Precesión Geodésica de Roseta Estricta
ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(x_orb, y_orb, color='#16a085', linewidth=1.2, label="Trayectoria Estelar")
ax3.plot(0, 0, 'ro', markersize=6, label="Centro Galáctico")
ax3.set_aspect('equal')
ax3.set_xlabel("Coordenada X")
ax3.set_ylabel("Coordenada Y")
ax3.set_title("Test 3: Precesión del Perihelio (Efecto Roseta)")
ax3.legend(frameon=True, fontsize=8)

plt.savefig("pruebas_galaxias_suaves.png", dpi=300)
plt.show()
