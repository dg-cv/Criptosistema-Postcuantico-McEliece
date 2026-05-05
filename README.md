# Códigos de Goppa y su aplicación al criptosistema post-cuántico McEliece

---

## Descripción

Este repositorio contiene el código desarrollado para el Trabajo de Fin de Grado en Ingeniería Matemática:

**“Códigos de Goppa y su aplicación al criptosistema post-cuántico McEliece”**

El objetivo del trabajo ha sido estudiar en profundidad los códigos de Goppa y su aplicación al criptosistema de McEliece, el cual ha resurgido en los últimos años debido a su potencial resistencia frente a ataques mediante computación cuántica.

Como parte central del proyecto, se ha realizado una implementación completa en Python tanto de:
- los códigos de Goppa  
- como del criptosistema de McEliece  

El código desarrollado supera las 3000 líneas, por lo que este repositorio permite organizarlo y facilitar su acceso.

---

## Contenido técnico

El proyecto incluye la implementación de:

- Operaciones en cuerpos finitos  
- Aritmética de polinomios cuyos coeficientes son elementos de cuerpos finitos 
- Operaciones matriciales sobre cuerpos finitos  
- Generación y decodificación de códigos de Goppa  
- Implementación del criptosistema de McEliece  

---

## Estructura del repositorio

El repositorio se compone de 13 archivos `.py`:

### Cuerpos finitos y polinomios
- `z_pz.py` : Operaciones en el cuerpo finito ℤ/pℤ  
- `z_pz_pol.py` : Polinomios en ℤ/pℤ[x]  
- `f_q.py` : Operaciones en Fq (q = pⁿ)  
- `f_q_pol.py` : Polinomios en Fq[x]  
- `f_q_m.py` : Operaciones en Fqᵐ  
- `f_q_m_pol.py` : Polinomios en Fqᵐ[x]  

### Álgebra lineal
- `mat_f_q.py` : Matrices sobre Fq  
- `mat_f_q_m.py` : Matrices sobre Fqᵐ  

### Códigos de Goppa y criptosistema
- `goppa.py` : Generación y decodificación de códigos de Goppa  
- `mceliece.py` : Implementación del criptosistema de McEliece  

### Ejemplos de uso
- `ejemplo_decodificacion_goppa_tfg.py`  
- `ejemplo_mceliece_tfg.py`  
- `ejemplo_mceliece_parametros_originales.py`  

---

## Requisitos

- Python 3.11.5  
- Librerías utilizadas:
  - `random`
  - `sympy`

---

## Ejecución

Las implementaciones de cuerpos finitos, polinomios y matrices son necesarias para ejecutar:

- `goppa.py` : generación y decodificación de códigos de Goppa  
- `mceliece.py` : implementación del criptosistema  

Para probar el funcionamiento, se incluyen los siguientes scripts:

- `ejemplo_decodificacion_goppa_tfg.py`  
- `ejemplo_mceliece_tfg.py`  
- `ejemplo_mceliece_parametros_originales.py`  

---

## Tiempos de ejecución

- Ejemplos estándar: menos de 1 minuto  
- `ejemplo_mceliece_parametros_originales.py`: aproximadamente 50 minutos  

---

## Resultados

El proyecto demuestra la viabilidad de implementar:

- Códigos de Goppa arbitrarios  
- El criptosistema de McEliece completo  

desde cero y con rigor matemático.

---

## Posibles mejoras

- Optimización del rendimiento computacional  
- Implementaciones más eficientes  
- Paralelización de ciertas operaciones  

---

## Autor

Diana Gómez Moreno  
Trabajo de Fin de Grado en Ingeniería Matemática Trabajo de Fin de Grado en Ingeniería Matemática de la Universidad Complutense de Madrid, 'Códigos de Goppa y su aplicación al criptosistema post-cuántico McEliece'.
