# Arquitectura del Motor de Geometría Diferencial (MGD)

## 🏛️ Visión general

El MGD es un sistema modular para el cálculo simbólico en geometría diferencial.

Su diseño separa completamente:

- Estructuras matemáticas
- Operaciones algebraicas
- Modelos físicos
- Construcción de curvatura

---

## 🧱 Capas del sistema

### 1. Núcleo geométrico (nucleo/)

Contiene las entidades base:

- Tensor
- ObjetoGeometrico
- Variedad

Responsabilidad:
Representar objetos matemáticos puros sin lógica algebraica pesada.

---

### 2. Operadores (operadores/)

Contiene las operaciones matemáticas:

- Contracción de índices
- Subida y bajada de índices
- Derivadas parciales
- Derivadas covariantes

Responsabilidad:
Ejecutar álgebra tensorial sin mutar objetos originales.

---

### 3. Modelos (modelos/)

Define teorías físicas específicas:

- Métrica disforme
- Campo escalar ψ
- Funciones A(ψ, ω), B(ψ, ω)

Responsabilidad:
Inyectar física al sistema.

---

### 4. Geometría diferencial (geometria/)

Construcciones avanzadas:

- Símbolos de Christoffel
- Tensor de Riemann
- Tensor de Ricci
- Tensor de Einstein

Responsabilidad:
Construcción de curvatura y geometría.

---

## ⚙️ Principios del sistema

- Inmutabilidad de objetos geométricos
- Separación estricta entre datos y operaciones
- Evaluación simbólica con SymPy
- Errores explícitos (no silenciosos)
- Diseño modular extensible

---

## 🔄 Flujo lógico

1. Definir variedad
2. Definir métrica
3. Construir tensores base
4. Aplicar operadores
5. Generar curvatura

---

## 🧠 Filosofía

El sistema representa la geometría como flujo computacional:

g → Γ → Riemann → Ricci

Sin ambigüedad matemática ni computacional.
