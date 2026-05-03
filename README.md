# TFG-Codigos-de-Goppa-y-su-aplicacion-al-criptosistema-post-cuantico-McEliece

# Objetivo del TFG y repositorio
Este repositorio contiene el código desarrollado para el Trabajo de Fin de Grado (TFG) en Ingeniería Matemática 'Códigos de Goppa y su aplicación al criptosistema post-cuántico McEliece'.

El objetivo de este TFG ha sido estudiar los códigos de Goppa y su aplicación al criptosistema de McEliece, que ha resurgido recientemente debido a su potencial, aunque aún no demostrada, resistencia frente a ataques realizados mediante ordenadores cuánticos. Como parte del trabajo, han sido implementados tanto los códigos de Goppa como el criptosistema de McEliece en el lenguaje de programación Python. Puesto que el código resultante supera las 3000 líneas, se ha creado este repositorio para organizarlo y facilitar su acceso.

## Estructura del repositorio
El repositorio se compone de un total de 13 archivos .py:
- z_pz.py : operaciones en el cuerpo finito de p elementos, con p número primo
- z_pz_pol.py : operaciones en el anillo de polinomios Z/pZ[x]
- f_q.py : operaciones en el cuerpo finito de q elementos Fq, con q=p^n, n>=1
- f_q_pol.py : operaciones en el anillo de polinomios Fq[x]
- f_q_m.py : operaciones en el cuerpo de q^m elementos Fq^m, con q=p^n, m>=1
- f_q_m_pol.py : operaciones en el anillo de polinomios Fq^m[x]
- mat_f_q_m.py : operaciones de matrices con elementos en el cuerpo finito Fq^m
- mat_f_q.py : operaciones de matrices con elementos en el cuerpo finito Fq
- goppa.py : funciones necesarias para la generación y decodificación de códigos de Goppa, según lo expuesto en la memoria del TFG. 
- mceliece.py : funciones necesarias para la implementación del criptosistema de McEliece, según lo expuesto en la memoria del TFG. 
- ejemplo_decodificacion_goppa_tfg.py : código del ejemplo de decodificación de un código de goppa expuesto en la memoria del TFG, capítulo 4
- ejemplo_mceliece_tfg.py : código del ejemplo de encriptado y desencriptado de McEliece expuesto en la memoria del TFG, capítulo 5
- ejemplo_mceliece_parametros_originales.py : ejemplo de encriptado y desencriptado de McEliece con los parámetros que este mismo sugirió en su artículo de 1978, donde expone dicho criptosistema.

## Requisitos
- Python 3.11.5
- Bibliotecas utilizadas: 'random', 'sympy'

## Ejecución
La implementación de las operaciones de cuerpos finitos, anillos de polinomios y matrices con elementos en cuerpos finitos es necesaria para poder ejecutar los archivos
goppa.py y mceliece.py, donde están definidas las funciones pertinentes a la generación y decodificación de códigos de Goppa arbitrarios según lo expuesto en la memoria del TFG, así como las funciones necesarias para definir correctamente el criptosistema de McEliece, respectivamente.

Se pueden probar las funciones anteriores en los archivos 'ejemplo_decodificacion_goppa_tfg.py', 'ejemplo_mceliece_tfg.py' y 'ejemplo_mceliece_parametros_originales.py'. 

Cabe mencionar que el tiempo de ejecución de 'ejemplo_decodificacion_goppa_tfg.py', 'ejemplo_mceliece_tfg.py' no llega al minuto, mientras que 'ejemplo_mceliece_parametros_originales.py' tarda aproximadamente 50 minutos en ejecutarse. 

## Autor
Trabajo realizado por [Diana Gómez Moreno] para el TFG en Ingeniería Matemática 'Códigos de Goppa y su aplicación al criptosistema post-cuántico McEliece'.
