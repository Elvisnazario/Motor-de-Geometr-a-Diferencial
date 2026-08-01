import numpy as np
import matplotlib
matplotlib.use('Agg')  # Previene congelamientos de terminal en VS Code/Linux
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad
from scipy.optimize import differential_evolution

plt.style.use('default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'grid.alpha': 0.3
})

# =====================================================================
# 1. CONSTANTES FÍSICAS (Unidades Astronómicas: kpc, km/s, M_sun)
# =====================================================================
class ConstantesFisicas:
    G = 4.30091e-6          # (km/s)^2 * kpc / M_sun
    c = 299792.458          # km/s

# =====================================================================
# 2. MODELO CONSISTENTE BASADO EN EL POTENCIAL INTEGRADO Phi(r)
# =====================================================================
class AtraccionEnergeticaModel:
    def __init__(self, M_bar, r_0, rho_0, r_c):
        self.M_bar = M_bar
        self.r_0 = r_0
        self.rho_0 = rho_0
        self.r_c = r_c
        self.G = ConstantesFisicas.G
        self.c = ConstantesFisicas.c

    def masa_efectiva(self, r):
        r = np.maximum(r, 1e-6)
        M_bar_reg = self.M_bar * (r**2) / (r**2 + self.r_0**2)
        x = r / self.r_c
        M_vac = 4.0 * np.pi * self.rho_0 * (self.r_c**3) * (x - np.arctan(x))
        return M_bar_reg + M_vac

    def dPhi_dr(self, r):
        """ Fuerza gravitacional por unidad de masa: dPhi/dr = G*M_eff(r)/r^2 (en (km/s)^2/kpc) """
        r = np.maximum(r, 1e-6)
        return (self.G * self.masa_efectiva(r)) / (r**2)

    def potencial_Phi(self, r):
        """ Potencial integrado Phi(r) = \int_0^r (G*M_eff(s)/s^2) ds """
        if np.isscalar(r):
            val, _ = quad(lambda s: self.dPhi_dr(s), 1e-6, max(r, 1e-6))
            return val
        else:
            return np.array([quad(lambda s: self.dPhi_dr(s), 1e-6, max(ri, 1e-6))[0] for ri in r])

    def Phi_adimensional(self, r):
        """ Potencial adimensional Phi/c^2 """
        return self.potencial_Phi(r) / (self.c**2)

    def dPhi_dr_adimensional(self, r):
        """ Gradiente adimensional (1/c^2) * dPhi/dr (1/kpc) """
        return self.dPhi_dr(r) / (self.c**2)

    def A_metric(self, r):
        """ Componente temporal de la métrica g_tt = -(1 + 2*Phi/c^2) """
        return 1.0 + 2.0 * self.Phi_adimensional(r)

    def dA_dr(self, r):
        """ Derivada exacta dA/dr = (2/c^2) * dPhi/dr """
        return 2.0 * self.dPhi_dr_adimensional(r)

    def velocidad_circular(self, r):
        """ Velocidad circular exacta v^2 = r * dPhi/dr / A(r) en km/s """
        r = np.maximum(r, 1e-6)
        A_val = self.A_metric(r)
        v2 = (r * self.dPhi_dr(r)) / A_val
        return np.sqrt(np.maximum(0.0, v2))

# =====================================================================
# 3. ECUACIONES GEODÉSICAS EN UNIDADES ADIMENSIONALES (c = 1)
# =====================================================================
def geodesicas_nulas(l, y, model, b):
    """
    Fotones viajando a c=1.
    y = [r, phi, dr/dl]
    b = parámetro de impacto en kpc
    """
    r, phi, vr = y
    A = model.A_metric(r)
    Ap = model.dA_dr(r) # en 1/kpc
    
    # Momento angular adimensional L = b (para fotón desde el infinito)
    L = b
    
    dr_dl = vr
    dphi_dl = L / (r**2)
    # Ecuación geodesica para fotones en límite débil (c=1)
    dvr_dl = -0.5 * Ap + (A * L**2) / (r**3) - 0.5 * Ap * (L**2) / (r**2)
    return [dr_dl, dphi_dl, dvr_dl]

def geodesicas_masivas(l, y, model, L_star):
    r, phi, vr = y
    A = model.A_metric(r)
    Ap = model.dA_dr(r)
    
    dr_dl = vr
    dphi_dl = L_star / (r**2)
    # Ecuación para partículas masivas lentas (v << c)
    dvr_dl = -0.5 * Ap + (A * L_star**2) / (r**3)
    return [dr_dl, dphi_dl, dvr_dl]

# =====================================================================
# 4. EJECUCIÓN Y VALIDACIÓN
# =====================================================================
print("="*75)
print("   MOTOR DE GEOMETRÍA DIFERENCIAL: MODELO NORMALIZADO (c = 1)")
print("="*75)

# --- TEST 1: SPARC ---
print("\n[TEST 1/3] Ajustando datos galácticos de SPARC...")

galaxias_data = {
    'NGC 3198': {
        'r': np.array([2.0, 4.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0]),
        'v': np.array([62.2, 115.7, 144.8, 152.8, 155.2, 156.9, 155.0, 152.0, 149.0, 147.0]),
        'err': np.array([5.5,  4.8,   4.2,   4.5,   4.9,   5.1,   5.3,   5.6,   5.8,   6.0])
    },
    'NGC 2403': {
        'r': np.array([1.0, 2.0, 3.5, 5.0, 7.0, 9.0, 11.5, 14.0, 17.0, 19.5]),
        'v': np.array([68.1, 95.3, 112.4, 122.1, 128.5, 131.0, 133.2, 134.1, 135.0, 134.7]),
        'err': np.array([4.1,  3.8,   3.5,   3.9,   4.2,   4.5,   4.6,   4.8,   5.0,   5.1])
    },
    'UGC 2259': {
        'r': np.array([0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]),
        'v': np.array([23.4, 45.1, 62.3, 71.5, 78.2, 82.1, 84.4, 85.1, 85.6, 85.2]),
        'err': np.array([2.5,  2.2,   2.1,   2.3,   2.5,   2.6,   2.8,   2.9,   3.0,   3.1])
    }
}

ajustes_resultados = {}

for nombre, data in galaxias_data.items():
    r_arr, v_obs, err_arr = data['r'], data['v'], data['err']
    
    bounds = [
        (1e9, 2e11),    # M_bar (M_sun)
        (0.5, 10.0),    # r_0 (kpc)
        (1e5, 5e7),     # rho_0 (M_sun / kpc^3)
        (2.0, 30.0)     # r_c (kpc)
    ]
    
    def loss_func(params):
        M_b, r0, rho0, rc = params
        model = AtraccionEnergeticaModel(M_b, r0, rho0, rc)
        v_theo = np.array([model.velocidad_circular(r) for r in r_arr])
        return np.sum(((v_obs - v_theo) / err_arr) ** 2)

    res = differential_evolution(loss_func, bounds, seed=42, polish=True, maxiter=150)
    M_b_opt, r0_opt, rho0_opt, rc_opt = res.x
    
    model_opt = AtraccionEnergeticaModel(M_b_opt, r0_opt, rho0_opt, rc_opt)
    r_dense = np.linspace(0.1, max(r_arr)*1.05, 100)
    v_dense = np.array([model_opt.velocidad_circular(rd) for rd in r_dense])
    
    v_fit = np.array([model_opt.velocidad_circular(r) for r in r_arr])
    residuos = v_obs - v_fit
    R2 = 1.0 - (np.sum(residuos**2) / np.sum((v_obs - np.mean(v_obs))**2))
    MAPE = np.mean(np.abs(residuos) / v_obs) * 100.0
    
    ajustes_resultados[nombre] = {
        'r_dense': r_dense, 'v_dense': v_dense,
        'r_obs': r_arr, 'v_obs': v_obs, 'err_obs': err_arr,
        'R2': R2, 'MAPE': MAPE, 'params': res.x, 'model': model_opt
    }
    print(f" -> {nombre:<10} | R² = {R2:.4f} | MAPE = {MAPE:.2f}% | M_b = {M_b_opt:.2e} | r_c = {rc_opt:.2f} kpc")

# --- TEST 2: Deflexión Óptica (Lente Gravitacional Real) ---
print("\n[TEST 2/3] Calculando trayectoria de fotones (Normalización c=1)...")
model_lens = ajustes_resultados['NGC 3198']['model']
light_tracks = []
deflection_angles = []
impact_parameters = [10.0, 15.0, 20.0]

for b in impact_parameters:
    r0_init = 40.0
    phi0_init = -np.arcsin(b / r0_init)
    vr0 = -np.cos(phi0_init)  # Componente radial para fotón plano viajando a c=1
    
    y0 = [r0_init, phi0_init, vr0]
    sol = solve_ivp(
        geodesicas_nulas, 
        (0.0, 80.0), 
        y0, 
        args=(model_lens, b), 
        rtol=1e-8, 
        atol=1e-10
    )
    
    r_pts, phi_pts = sol.y[0], sol.y[1]
    x_pts, y_pts = r_pts * np.cos(phi_pts), r_pts * np.sin(phi_pts)
    light_tracks.append((x_pts, y_pts, b))
    
    # Deflexión analítica de lente en campo débil: alpha = 4 * G * M_eff / (c^2 * b)
    M_eff_b = model_lens.masa_efectiva(b)
    alpha_rad = (4.0 * ConstantesFisicas.G * M_eff_b) / ((ConstantesFisicas.c**2) * b)
    alpha_arcsec = alpha_rad * (180.0 / np.pi) * 3600.0
    
    print(f" -> Parámetro b = {b:.1f} kpc | M_eff({b}kpc) = {M_eff_b:.2e} M_sun | Ángulo de deflexión: {alpha_arcsec:.2f} arcsec ({alpha_rad:.2e} rad)")

# --- TEST 3: Precesión Geodésica ---
print("\n[TEST 3/3] Simulando órbita estelar...")
r0_orb = 8.0
v_circ_target = model_lens.velocidad_circular(r0_orb) / ConstantesFisicas.c # En unidades de c
L_star = r0_orb * (v_circ_target * 0.95)

sol_orb = solve_ivp(
    geodesicas_masivas, 
    (0.0, 2000.0), 
    [r0_orb, 0.0, 0.0], 
    args=(model_lens, L_star), 
    rtol=1e-8, 
    atol=1e-10
)

x_orb = sol_orb.y[0] * np.cos(sol_orb.y[1])
y_orb = sol_orb.y[0] * np.sin(sol_orb.y[1])
print(" -> Simulación de órbita estelar completada.")

# =====================================================================
# 5. VISUALIZACIÓN Y GUARDADO
# =====================================================================
fig = plt.figure(figsize=(15, 5), constrained_layout=True)
gs = plt.GridSpec(1, 3, figure=fig)

# Panel 1: SPARC
ax1 = fig.add_subplot(gs[0, 0])
colors = ['#e74c3c', '#3498db', '#2ecc71']
for idx, (nombre, r_fit) in enumerate(ajustes_resultados.items()):
    ax1.errorbar(r_fit['r_obs'], r_fit['v_obs'], yerr=r_fit['err_obs'], fmt='o', color=colors[idx], alpha=0.6, capsize=3)
    ax1.plot(r_fit['r_dense'], r_fit['v_dense'], '-', color=colors[idx], linewidth=2.0, label=f"{nombre} ($R^2$: {r_fit['R2']:.3f})")
ax1.set_xlabel("Radio r (kpc)")
ax1.set_ylabel("Velocidad v (km/s)")
ax1.set_title("Test 1: SPARC (Potencial Integrado)")
ax1.legend(frameon=True, loc="lower right", fontsize=8)

# Panel 2: Deflexión Óptica
ax2 = fig.add_subplot(gs[0, 1])
lens_colors = ['#2980b9', '#8e44ad', '#f1c40f']
for i, (x, y, b) in enumerate(light_tracks):
    ax2.plot(x, y, color=lens_colors[i], linewidth=1.5, label=f"b = {b:.1f} kpc")
core = plt.Circle((0, 0), model_lens.r_c, color='#2c3e50', alpha=0.2, label=f"Núcleo ($r_c = {model_lens.r_c:.1f}$ kpc)")
ax2.add_patch(core)
ax2.set_xlim(-40, 40)
ax2.set_ylim(-20, 20)
ax2.set_aspect('equal')
ax2.set_xlabel("X (kpc)")
ax2.set_ylabel("Y (kpc)")
ax2.set_title("Test 2: Deflexión Óptica Corregida")
ax2.legend(frameon=True, fontsize=8)

# Panel 3: Precesión
ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(x_orb, y_orb, color='#16a085', linewidth=1.0, label="Órbita Estelar")
ax3.plot(0, 0, 'ro', markersize=5, label="Centro Galáctico")
ax3.set_aspect('equal')
ax3.set_xlabel("X (kpc)")
ax3.set_ylabel("Y (kpc)")
ax3.set_title("Test 3: Precesión Geodésica")
ax3.legend(frameon=True, fontsize=8)

plt.savefig("pruebas_galaxias_suaves.png", dpi=300)
print("\n[ÉXITO] Unidades normalizadas. Imagen guardada como 'pruebas_galaxias_suaves.png'.")
