# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Las matrices están definidas por filas

Contiene las siguientes funciones: 
    mat_h(g,a,p,f,h)
    mult_mat(a,b,p,f,h)    
"""
import f_q_m, f_q_m_pol

def mat_h(g,a,p,f,h):
    '''
    Parameters
    ----------
    g : LIST
        Lista de listas de listas que representa un polinomio de grado t con 
        coeficientes en el cuerpo Fq^m
    a : LIST
        lista de listas de listas, donde cada elemento de la lista principal
        es un elemento del cuerpo Fq^m. 
        concretamente, son cada uno de los 'l' ai's necesarios para definir 
        el código de Goppa C
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ
    h : LIST
        lista de listas que representa un polinomio de grado m >= 1 irreducible, 
        mónico, con coeficientes en el cuerpo Fq
    
    Returns
    -------
    Lista de listas de listas de listas, que representa una matriz de orden txl, 
    definida por filas, cuyos elementos pertenecen al cuerpo Fq^m.
    Es una de las expresiones de la matriz de control de paridad del código 
    de Goppa C definido con los elementos de 'a' y polinomio g(x) del anillo Fq^m[x]
    '''
    t = len(g) - 1
    l = len(a)
    
    H = [0]*t
    
    g_a_inv = [0]*l
    
    for i in range(l):
        g_ai = f_q_m_pol.eval_pol(g, a[i], p, f, h)
        g_a_inv[i] = f_q_m.inv_mult(g_ai, p, f, h)
        
    # Vamos definiendo una a una las t filas de la matriz H
    for i in range(t):
        H[i] = [0]*l
        
        # Los elementos de la fila i son los siguientes: 
        for j in range(l):
            pot_aj = f_q_m.pot(a[j], i, p, f, h)
            H[i][j] = f_q_m.mult(g_a_inv[j], pot_aj, p, f, h)            
            
    return (H)


def mult_mat(a,b,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas de listas que representa una matriz de dimensión 
        fa x ca definida por filas, cuyos elementos pertenecen al cuerpo Fq^m
    b : LIST
        Lista de listas de listas de listas que representa una matriz de dimensión 
        ca x cb definida por filas cuyos elementos pertenecen al cuerpo Fq^m
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ
    h : LIST
        lista de listas que representa un polinomio de grado m >= 1 irreducible, 
        mónico con coeficientes en el cuerpo Fq

    Returns
    -------
    sol, lista de listas de listas de listas, que representa una matriz de 
    dimensión fa x cb cuyos elementos pertenecen al cuerpo Fq^m, 
    resultado del producto de matrices a*b
    '''
    fa = len(a) # número de filas de a
    ca = len(a[0]) # número de columnas de a
    cb = len(b[0]) # número de columnas de b
    sol = [0] * fa # definimos el número de filas de sol
    
    for i in range(fa):
        sol[i] = [f_q_m.cero(p, f, h)] * cb # inicializamos cada fila para que contenga en cada entrada el elemento cero del cuerpo Fq
    
    for i in range(fa):
        for j in range(cb):
            for k in range(ca):
                prod = f_q_m.mult(a[i][k], b[k][j], p, f, h)
                sol[i][j] = f_q_m.suma(sol[i][j], prod, p, f, h)
    
    return(sol)
