# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Contiene las siguientes funciones: 
    gen(g,a,p,f,h)
    enc(v,Gpu,t,p,f)
    dec(v,P,S,a,g,p,f,h)    

Funciones auxiliares: 
    incluir_error(c,t,p,f)
    permut(n,p,f)
"""
import mat_f_q, goppa, f_q, random

def gen(g,a,p,f,h):
    '''
    Parameters
    ----------
    g : LIST
        Lista de listas de listas que representa un polinomio de grado t con 
        coeficientes en el cuerpo Fq^m
    a : LIST
        lista de listas de listas, donde cada elemento de la lista principal
        es un elemento de Fq^m.
        concretamente, son cada uno de los 'l' ai's necesarios para definir 
        el código de Goppa C
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
    sk y pk, donde sk es la clave secreta del criptosistema de McEliece definido
    con los parámetros anteriores, y pk es la clave pública del mismo
    '''
    G = goppa.mat_g_can(g, a, p, f, h)
    P = permut(len(G),p,f)
    k = len(G[0]) # dimensión del código
    S = mat_f_q.mat_cuad_al_inv(k, p, f)[0]
    t = len(g) - 1
    
    Gpk_0 = mat_f_q.mult_mat(G,S,p,f)
    Gpk = mat_f_q.mult_mat(P,Gpk_0,p,f)

    sk = [G,P,S,g,a]
    pk = [Gpk,t]
    return(sk,pk)
    

def enc(v,Gpu,t,p,f):
    '''
    Parameters
    ----------
    v : LIST
        lista de listas que representa un vector columna cuyos elementos pertenecen
        al cuerpo Fq.
        es el mensaje que se quiere transmitir
    Gpu : LIST
        lista de listas de listas que representa una matriz cuyos elementos pertenecen
        al cuerpo Fq.
        es la clave pública del criptosistema de McEliece
    t : INT
        grado del polinomio g que define el código de Goppa en el que se basa
        el encriptado
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ

    Returns
    -------
    v3, lista de listas que representa un vector columna cuyos elementos pertenecen
    al cuerpo Fq. Es el mensaje que se quería transmitir, encriptado según McEliece
    '''
    v_mat = mat_f_q.vect_a_mat(v)
    v1 = mat_f_q.mult_mat(Gpu,v_mat,p,f) # v1 = (G clave pública) * v
    
    v2 = mat_f_q.mat_a_vect(v1)
    v3 = incluir_error(v2,t,p,f)[0] # v encriptado = m + e, con e error aleatorio 
    return(v3)


def dec(v,P,S,a,g,p,f,h):
    '''
    Parameters
    ----------
    v : LIST
        lista de listas que representa un vector columna cuyos elementos pertenecen
        al cuerpo Fq.
        es un mensaje encriptado
    P : LIST
        lista de listas de listas que representa una matriz de dimensión lxl, con l
        longitud del código en el que se basa el encriptado, cuyos elementos pertenecen
        al cuerpo Fq.
        es una permutación
    S : LIST
        lista de listas de listas que representa una matriz de dimensión kxk, con k
        dimensión del código en el que se basa el encriptado, cuyos elementos pertenecen
        al cuerpo Fq.
        es una matriz invertible 
    a : LIST
        lista de listas de listas, donde cada elemento de la lista principal
        es un elemento de Fq^m.
        concretamente, son cada uno de los 'l' ai's necesarios para definir 
        el código de Goppa C
    g : LIST
        Lista de listas de listas que representa un polinomio de grado t con 
        coeficientes en el cuerpo Fq^m
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
    corr, lista de listas que representa un vector columna cuyos elementos pertenecen
    al cuerpo Fq. Es el mensaje que se quería transmitir, desencriptado según McEliece
    '''
    l = len(P)
    k = len(S)
    
    P_inv = mat_f_q.mat_inv(P,p,f) # la inversa de la permutación P
    S_inv = mat_f_q.mat_inv(S,p,f) # la inversa de la matriz S
    
    v0 = mat_f_q.vect_a_mat(v)
    v1_0 = mat_f_q.mult_mat(P_inv, v0, p, f) # P^{-1}*v
    v1 = mat_f_q.mat_a_vect(v1_0)
    
    # Eliminamos el error (corregimos la palabra recibida para que pertenezca al código)
    t = len(g) - 1 # grado del polinomio g 
    v2_0,e = goppa.decode(v1,a,g,t,p,f,h)
    
    # Debemos quedarnos con S*m, y por la estructura de la matriz G, 
    # sabemos qué coordenadas forman el mensaje original
    Sm_0 = v2_0[l-k:] 
    for i in range(len(Sm_0)):
        Sm_0[i] = f_q.inv_adit(Sm_0[i], p, f)
    
    Sm = mat_f_q.vect_a_mat(Sm_0)
    corr_0 = mat_f_q.mult_mat(S_inv, Sm, p, f)
    corr = mat_f_q.mat_a_vect(corr_0)
    return(corr)


# Funciones auxiliares:
def incluir_error(c,t,p,f):
    '''
    Parameters
    ----------
    c : LIST
        lista de listas que representa un vector columna cuyos elementos pertenecen
        al cuerpo Fq
    t : INT
        número entero >= 1
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ

    Returns
    -------
    r, peso, error; peso es un int, r y error son listas de listas que representan 
    vectores columna cuyos elementos pertenecen al cuerpo Fq, tales que: el peso de Hamming 
    de error es 'peso', y r = c + e
    '''
    l = len(c)
    r = goppa.copy_vec(c) # aquí se guardará c tras añadirle un error
    error = [f_q.cero(p,f)] * l # aquí almacenaremos los términos de error
    max_err = (t) // 2 # máximo número de errores que podemos decodificar con 
                       # nuestro código
                       
    peso = 0
    
    for i in range(l):
        rand = random.randint(0, 4)
        
        if (rand == 1 and peso < max_err):
            error[i] = f_q.rand_elem(p, f)
            r[i] = f_q.suma(c[i],error[i],p,f)
            peso += 1
        
        if peso == max_err:
            break
    
    return(r,peso,error) 


def permut(n,p,f):
    '''
    Parameters
    ----------
    n : INT
        dimensión de la matriz resultante
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ

    Returns
    -------
    P, lista de listas de listas que representa una matriz de orden nxn, 
    definida por filas, cuyos elementos pertenecen al cuerpo Fq.
    '''
    P = [0] * n
    
    for i in range(n):
        P[i] = []
        j = 0
        while j < n:
            if i==j:
                P[i].append([1])
            else:
                P[i].append([])
            j+=1
            
    for w in range(50*n):
        f1 = random.randint(0,n-1)
        f2 = random.randint(0,n-1)
    
        while f1 == f2:
            f2 = random.randint(0,n-1)
    
        pf1 = P[f1]
        pf2 = P[f2]
    
        P[f1] = pf2
        P[f2] = pf1

    return(P)
