# Motor de Geometría Diferencial (MGD)

## 🧠 Descripción general

El **Motor de Geometría Diferencial (MGD)** es una biblioteca simbólica en desarrollo diseñada para la representación, manipulación y cálculo de estructuras fundamentales de la geometría diferencial y la relatividad general.

El objetivo del proyecto es construir un sistema computacional riguroso, modular y extensible capaz de trabajar con:

- Tensores de cualquier rango
- Métricas pseudo-riemannianas generales
- Conexiones afines (Símbolos de Christoffel)
- Tensores de curvatura (Riemann, Ricci y Einstein)
- Geometrías deformadas y modelos gravitacionales modificados

---

## 🏛️ Fundamentación matemática

El motor está basado en la teoría estándar de geometría diferencial sobre variedades:

- Variedad diferenciable \( \mathcal{M} \)
- Espacio tangente \( T_p\mathcal{M} \)
- Métrica \( g_{\mu\nu} \)
- Conexión de Levi-Civita (cuando aplica)
- Convención de suma de Einstein

Se contempla la extensión futura hacia geometrías con:

- Torsión
- No metricidad
- Métricas disformes
- Teorías gravitacionales alternativas

---

## ⚙️ Arquitectura del sistema

El MGD se estructura en tres capas principales:

### 1. Núcleo geométrico

Contiene las entidades matemáticas fundamentales:

- Tensor
- ObjetoGeometrico
- Variedad

Responsable de la representación estructural de los objetos geométricos.

---

### 2. Capa de operadores

Encargada de la manipulación algebraica:

- Contracción de índices
- Subida y bajada de índices
- Derivadas (ordinarias y covariantes)
- Simplificación simbólica

---

### 3. Capa de geometría diferencial

Construcción de objetos derivados:

- Símbolos de Christoffel
- Tensor de Riemann
- Tensor de Ricci
- Tensor de Einstein

---

## 🧩 Filosofía del diseño

El proyecto sigue los siguientes principios:

- Separación estricta entre datos y operaciones
- Representación simbólica mediante SymPy
- Inmutabilidad de los objetos geométricos
- Errores explícitos en lugar de resultados silenciosos
- Modularidad extrema del sistema

---

## ⚠️ Nota importante sobre conexiones

Los Símbolos de Christoffel:

- ❌ NO son tensores
- ❌ NO transforman como tensores
- ✔ Son conexiones afines

Por lo tanto, se implementan como objetos geométricos independientes dentro del sistema.

---

## 📦 Estado del proyecto

**Versión:** v0.1.0 (fase de diseño)

Actualmente el proyecto se encuentra en etapa de:

- Definición de arquitectura
- Especificación matemática
- Diseño de interfaces
- Preparación del núcleo tensorial

---

## 👤 Autor

**Investigador Principal:**  
Elvis Omar Nazario Espinoza

---

## 🛠️ Herramientas utilizadas

- Python
- SymPy
- GitHub
- Asistencia de herramientas de inteligencia artificial (ChatGPT)

---

## 📌 Agradecimientos

Durante el diseño de la arquitectura y la documentación se utilizaron herramientas de inteligencia artificial como apoyo conceptual y técnico.

Todas las decisiones finales de diseño e implementación pertenecen al autor del proyecto.
