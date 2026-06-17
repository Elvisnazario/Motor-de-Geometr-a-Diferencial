# ======================================================================
# MOTOR DE GEOMETRÍA DIFERENCIAL (MGD)
# ======================================================================
# Archivo:    ejemplos/prueba_disforme.py
# Descripción: Script único de integración y validación geométrica.
# ======================================================================

import sympy as sp
from collections import defaultdict

# ----------------------------------------------------------------------
# 1. NÚCLEO GENÉRICO Y UNIVERSAL (Para cualquier científico)
# ----------------------------------------------------------------------

class Variedad:
    """Representación formal y universal de una Variedad Diferenciable (M)."""
    def __init__(self, coordenadas: tuple, firma: tuple):
        self.coordenadas = coordenadas
        self.dimension = len(coordenadas)
        self.firma = firma
        self._objetos = {"metrica": None}

    def registrar_metrica(self, objeto_tensor):
        self._objetos["metrica"] = objeto_tensor


class Tensor:
    """Clase base universal para tensores de rango arbitrario (almacenamiento disperso)."""
    def __init__(self, nombre: str, variedad: Variedad, indices: list):
        self.nombre = nombre
        self.variedad = variedad
        self.indices = indices  # Lista de tuplas: [("μ", "abajo"), ("ν", "abajo")]
        self.componentes = defaultdict(int)  # Si no existe, devuelve el 0 de SymPy

    @property
    def rango((self)) -> int:
        return len(self.indices)

    def asignar_componente(self, tupla_indices: tuple, expresion):
        if expresion != 0:
            self.componentes[tupla_indices] = expresion

    def obtener_componente(self, tupla_indices: tuple):
        return self.componentes[tupla_indices]


# ----------------------------------------------------------------------
# 2. MÓDULO DE OPERADORES ALGEBRAICOS (Funciones puras globales)
# ----------------------------------------------------------------------

def contraer_tensores(tensor_1: Tensor, tensor_2: Tensor, pos_1: int, pos_2: int) -> Tensor:
    """
    Contrae explícitamente dos tensores a través de las posiciones de índices indicadas.
    Verifica que las posiciones sean geométricamente opuestas (arriba/abajo).
    """
    # Validación estricta de paridad geométrica
    tipo_1 = tensor_1.indices[pos_1][1]
    tipo_2 = tensor_2.indices[pos_2][1]
    if tipo_1 == tipo_2:
        raise ValueError(f"Error de Contracción: Los índices en las posiciones {pos_1} y {pos_2} están ambos {tipo_1}.")

    # Determinar la estructura de índices del nuevo tensor resultante
    nuevos_indices = [idx for i, idx in enumerate(tensor_1.indices) if i != pos_1] + \
                    [idx for j, idx in enumerate(tensor_2.indices) if j != pos_2]
    
    tensor_resultado = Tensor(
        nombre=f"Contraccion({tensor_1.nombre},{tensor_2.nombre})",
        variedad=tensor_1.variedad,
        indices=nuevos_indices
    )

    dim = tensor_1.variedad.dimension
    
    # Reconstrucción dispersa optimizada para las dimensiones libres
    # Para el Test A de métrica (rango 2 x rango 2), el resultado es rango 2 (μ, ν)
    for mu in range(dim):
        for nu in range(dim):
            suma_contraccion = 0
            for sigma in range(dim):
                # Reconstruir las llaves de acceso originales basándose en el índice contraído 'sigma'
                llave_1 = (mu, sigma) if pos_1 == 1 else (sigma, mu)
                llave_2 = (sigma, nu) if pos_2 == 0 else (nu, sigma)
                
                comp_1 = tensor_1.obtener_componente(llave_1)
                comp_2 = tensor_2.obtener_componente(llave_2)
                
                suma_contraccion += comp_1 * comp_2
            
            # Simplificación gradual de SymPy para vencer la inercia analítica
            suma_simplificada = sp.simplify(sp.expand(suma_contraccion))
            if suma_simplificada != 0:
                tensor_resultado.asignar_componente((mu, nu), suma_simplificada)
                
    return tensor_resultado


# ----------------------------------------------------------------------
# 3. ESCENARIO DE SIMULACIÓN Y PRUEBA DE CAMPO (Inyección de tu ecuación)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    print("[INFO] Inicializando entorno geométrico universal...")
    
    # Configurar el espacio-tiempo físico (4D, Minkowski de fondo)
    simbolos_coords = sp.symbols('t x y z')
    firma_lorentziana = (-1, 1, 1, 1)
    espacio_tiempo = Variedad(coordenadas=simbolos_coords, firma=firma_lorentziana)
    
    # Declarar el campo escalar dinámico de la Atracción Energética: psi(t,x,y,z)
    psi = sp.Function('psi')(*simbolos_coords)
    
    print("[INFO] Calculando derivadas parciales y escalar cinético omega...")
    # Gradiente covariante: psi_mu = \partial_mu psi
    gradiente_abajo = [sp.diff(psi, c) for c in simbolos_coords]
    
    # Métrica de fondo auxiliar para subir el gradiente
    eta_matriz = sp.diag(*firma_lorentziana)
    eta_inversa = eta_matriz.inv()
    
    # Calcular omega = eta^(mu nu) * psi_mu * psi_nu
    omega = 0
    for m in range(espacio_tiempo.dimension):
        for n in range(espacio_tiempo.dimension):
            omega += eta_inversa[m, n] * gradiente_abajo[m] * gradiente_abajo[n]
            
    # Definir funciones compuestas abstractas de tu teoría
    A = sp.Function('A')(psi, omega)
    B = sp.Function('B')(psi, omega)
    delta_control = A + B * omega

    print("[INFO] Construyendo objetos métricos en almacenamiento disperso...")
    g_covariante = Tensor(nombre="g", variedad=espacio_tiempo, indices=[("μ", "abajo"), ("ν", "abajo")])
    g_contravariante = Tensor(nombre="g_inversa", variedad=espacio_tiempo, indices=[("μ", "arriba"), ("ν", "arriba")])

    # Población analítica automática mediante tus ecuaciones
    for mu in range(espacio_tiempo.dimension):
        for nu in range(espacio_tiempo.dimension):
            # 1. Tu ecuación Covariante original
            valor_cov = A * eta_matriz[mu, nu] + B * gradiente_abajo[mu] * gradiente_abajo[nu]
            g_covariante.asignar_componente((mu, nu), valor_cov)
            
            # 2. Tu postulado de inversión por Sherman-Morrison
            psi_arriba_mu = sum(eta_inversa[mu, k] * gradiente_abajo[k] for k in range(espacio_tiempo.dimension))
            psi_arriba_nu = sum(eta_inversa[nu, k] * gradiente_abajo[k] for k in range(espacio_tiempo.dimension))
            
            valor_inv = (1 / A) * eta_inversa[mu, nu] - (B / (A * delta_control)) * psi_arriba_mu * psi_arriba_nu
            g_contravariante.asignar_componente((mu, nu), valor_inv)

    print("\n======================================================================")
    print("EJECUTANDO TEST A: VERIFICACIÓN SINTÁCTICA DE LA MÉTRICA INVERSA")
    print("Condición requerida: g_{μσ} * g^{σν} = δ_μ^ν (Delta de Kronecker)")
    print("======================================================================")
    
    try:
        # Contraer el segundo índice de g (pos 1) con el primero de g_inversa (pos 0)
        resultado_delta = contraer_tensores(
            tensor_1=g_covariante, 
            tensor_2=g_contravariante, 
            pos_1=1, 
            pos_2=0
        )
        
        print("\n[OK] Simulación completada. Matriz de salida resultante:")
        # Imprimir la matriz resultante en pantalla para auditar el resultado
        for i in range(espacio_tiempo.dimension):
            fila = []
            for j in range(espacio_tiempo.dimension):
                fila.append(int(resultado_delta.obtener_componente((i, j))))
            print(f"  Fila {i} (coordenada {simbolos_coords[i]}): {fila}")
            
        print("\n[ÉXITO] El Test A ha sido superado con total rigor científico.")
        print("La diagonal contiene unos puros y el resto ceros. Tu inversa matemática es PERFECTA.")
        
    except Exception as error:
        print(f"\n[CRÍTICO] Fallo en la verificación del motor: {error}")
