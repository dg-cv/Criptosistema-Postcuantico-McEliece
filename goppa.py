# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Contiene las siguientes funciones: 
    gen_ai(g,l,p,f,h)
    mat_h(g,a,p,f,h)
    mat_h_stan(g,a,p,f,h)
    mat_g_can(g,a,p,f,h)
    sindrome(r,a,g,p,f,h)
    encode(m,G,p,f)
    decode(r,a,g,k,p,f,h)

Funciones auxiliares: 
    gcd_truncado(a,b,k,p,f,h) 
    copy_vec(v)
"""
import f_q, mat_f_q, f_q_m, f_q_m_pol

def gen_ai(g,l,p,f,h):
    '''
    Parameters
    ----------
    g : LIST
        Lista de listas de listas que representa un polinomio de grado t con 
        coeficientes en el cuerpo Fq^m
    l : INT
        Número de ai's que queremos generar, con cada ai un elemento del cuerpo Fq^m
        tal que g(ai)/=0, y con todos los ai's distintos entre sí. 
        l es la longitud de las palabras del código
    p : INT
        número primo >= 2
    f : LIST
        lista que representa un polinomio de grado n >= 1 irreducible, mónico, 
        con coeficientes en el cuerpo Z/pZ
    h : LIST
        lista de listas que representa un polinomio de grado m >= 1 irreducible, 
        mónico con coeficientes en el cuerpo Fq

    Returns
    -------
    ai, una lista de l elementos, donde cada elemento es uno de los ai's mencionados 
    anteriormente. Así, es una lista de listas de listas
    '''  
    ai = [0]*l
    
    j = 0
    while j < l:
        # Generamos un posible candidato para un nuevo ai:
        ai[j] = f_q_m.rand_elem(p, f, h)
        dif = True
        
        for k in range(j):
            if ai[k] == ai[j]:
                 dif = False 

        if not dif:
            continue

        # Para añadirlo a la lista, debe no ser raíz de g y ser diferente a los ai ya generados
        no_raiz = (f_q_m_pol.eval_pol(g, ai[j], p, f, h) != f_q_m.cero(p, f, h))
        
        if no_raiz:
            j = j + 1

    return(ai)


def mat_h(g,a,p,f,h):
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
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes 
        en el cuerpo Z/pZ
    h : LIST
        lista de listas que representa un polinomio de grado m >= 1 irreducible, 
        mónico, con coeficientes en el cuerpo Fq

    Returns
    -------
    Lista de listas de listas que representa una matriz de dimensión mtxl, definida por filas,
    cuyos elementos pertenecen al cuerpo Fq.
    Es una de las expresiones de la matriz de control de paridad del código 
    de Goppa C definido con los elementos de 'a' y polinomio g
    '''
    return(mat_f_q.mat_h(g,a,p,f,h))


def mat_h_stan(g,a,p,f,h):
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
    Lista de listas de listas que representa una matriz definida por filas,
    cuyos elementos pertenecen al cuerpo Fq.
    Es la forma canónica de la matriz de control de paridad del código de Goppa C 
    definido con los elementos de 'a' y polinomio g
    '''
    H = mat_h(g,a,p,f,h)

    m = len(h) - 1 # grado de h
    t = len(g) - 1 # grado de g
    l = len(a) # número de ai's
    
    fi = 0 # contador de la fila de H en la que estamos
    continuar = True # valdrá False cuando tengamos H en su forma escalonada reducida
    
    while continuar: 
        if fi >= m*t: # el caso en que tiene rango máximo de entre los posibles
            continuar = False
            fi = fi - 1
            
        elif H[fi][fi] == f_q.cero(p, f): # el elemento (fi,fi) de la diagonal es el cero del cuerpo f_q    
                                        
            # Vemos con un bucle si algún elemento de esa columna es dinstinto de cero, pues permitiremos intercambio de filas
            # para hacer la reducción gaussiana
            
            change = False # nos indicará si hay que hacer intercambio de filas
            r = fi + 1 
            while (r<m*t and (not change)):
                if H[r][fi] != f_q.cero(p,f):
                    change = True
                else:
                    r = r + 1
                
            if change: # Podemos continuar con la reducción gaussiana intercambiando las filas fi y r
                # Intercambiamos las filas:
                fi_orig = H[fi]
                fr_orig = H[r]
                
                H[fi] = fr_orig
                H[r] = fi_orig
                
                # Realizamos la siguiente resta para continuar la reducción gaussiana por la misma fila fi, que
                # es la antigua fila r, en la siguiente vuelta del bucle
                fi = fi - 1

            else: # No podemos continuar con la reducción gaussiana; veamos si hemos acabado o si tocaría intercambiar columnas
                  # cosa que no haremos, parando aquí el proceso
                for i in range(fi,m*t): # empezando en la propia fila fi
                    for j in range(fi+1,l): # las siguientes columnas
                        if H[i][j] != f_q.cero(p, f):
                            raise ValueError('No es posible calcular H en forma canónica sin intercambiar columnas')
                        
                
                continuar = False # la matriz está en forma triangular y salimos de este primer bucle

        else: 
            inv1 = f_q.inv_mult(H[fi][fi], p, f)
            
            for i in range(fi,l): # multiplicamos la fila fi por el inverso del elemento H[fi][fi] para que su primer elemento sea = f_q.uno()
                H[fi][i] = f_q.mult(H[fi][i], inv1, p, f)
            
                
            
            for i in range(fi+1,m*t): # eliminamos todos los elementos en la columna fi mediante resta de filas
                
                if H[i][fi] != f_q.cero(p, f):
                    
                    for j in range(l-1,fi-1,-1):
                        mult = f_q.mult(H[fi][j], H[i][fi], p, f)
                        inv_adit = f_q.inv_adit(mult, p, f)
                        H[i][j] = f_q.suma(H[i][j], inv_adit, p, f)
                
        
        # Seguimos con la reducción gaussiana
        fi = fi + 1

    for ix in range(len(H)): # Tenemos solo la matriz en forma escalonada, queremos su forma escalonada reducida
        i = len(H)-1-ix  
        for j in range(i): # todas las filas a las que les deberemos restar la fila i
            for k in range(l-1,i-1,-1):
                mult = f_q.mult(H[i][k],H[j][i], p, f)
                inv_adit = f_q.inv_adit(mult, p, f)
                H[j][k] = f_q.suma(H[j][k], inv_adit, p, f) 
           
    return H


def mat_g_can(g,a,p,f,h):
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
    Lista de listas de listas que representa una matriz, definida por filas,
    cuyos elementos pertenecen al cuerpo Fq.
    Es la forma canónica de la matriz generadora del código de Goppa C 
    definido con los elementos de 'a' y polinomio g
    '''
    H = mat_h_stan(g,a,p,f,h)
    
    l = len(a) # número de ai's = número de filas de G
    l_k = len(H) # número de filas de H = número de columnas de G
    k = l - l_k # número de columnas de G
    
    G = [0] * l # dotamos a G de su número de filas

    for f in range(l):
        G[f] = [f_q.cero(p, f)] * (k) # inicializamos cada fila de modo que contenga en cada entrada el elemento cero del cuerpo Fq
    
    # Las primeras l_k filas son las de la matriz A
    for f in range(l_k):
        for c in range(k):
            G[f][c] = H[f][l_k + c]
    
    # Las últimas k son -Id(k):
    i = 0 # índice del elemento de la matriz -Id no nulo
    for f in range(l_k,l):
        G[f][i] = f_q.inv_adit(f_q.uno(p, f), p, f)
        i = i + 1
        
    return(G)


def sindrome(r,a,g,p,f,h):
    '''
    Parameters
    ----------
    r : LIST
        Vector columna de tamaño l de elementos del cuerpo Fq que representa una palabra
        recibida. Es, por tanto, una lista de listas. 
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
    sind, una lista de listas de listas que representa el polinomio con coeficientes en
    el cuerpo Fq^m conocido como síndrome
    '''
    sind = f_q_m_pol.cero(p, f, h)
    
    for i in range(len(a)):
        ri = [r[i]] # ri como elemento de Fq^m, para poder hacer las operaciones pertinentes
        pol_ai = [f_q_m.inv_adit(a[i],p,f,h) , f_q_m.uno(p, f, h)] # x - ai
        inv_i = f_q_m_pol.gcd_ext(pol_ai, g, p, f, h)[1] # inverso módulo g(x) de x - ai
        mult_i = f_q_m_pol.mult([ri], inv_i, p, f, h) # ri / (x - ai); hay que tratar ri 
                                                      # como polinomio de grado 0 de Fq^m[x]
        sumando = f_q_m_pol.div_pol(mult_i, g, p, f, h)[1] # (ri / (x - ai)) mod (g(x))
        sind = f_q_m_pol.suma(sind, sumando, p, f, h)

    return(sind)


def encode(m,G,p,f):
    '''
    Parameters
    ----------
    m : LIST
        Vector columna de tamaño k de elementos del cuerpo Fq que representa un mensaje
        que se quiere transmitir
    G : LIST
        Lista de listas de listas que representa una matriz, definida por filas,
        cuyos elementos pertenecen al cuerpo Fq. Es una matriz generadora de un 
        código de Goppa C
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ

    Returns
    -------
    r, vector columna de tamaño l de elementos del cuerpo Fq que representa la palabra
    m codificada para pertenecer al código de Goppa C de matriz generadora G
    '''
    m_mat = mat_f_q.vect_a_mat(m)
    m1 = mat_f_q.mult_mat(G,m_mat,p,f)
    r = mat_f_q.mat_a_vect(m1)
    return(r)


def decode(r,a,g,k,p,f,h):
    '''
    Parameters
    ----------
    r : LIST
        Vector columna de tamaño l de elementos del cuerpo Fq que representa una palabra
        recibida. Esta puede pertenecer al código o no. 
    a : LIST
        lista de listas de listas, donde cada elemento de la lista principal
        es un elemento de Fq^m.
        concretamente, son cada uno de los 'l' ai's necesarios para definir 
        el código de Goppa C
    g : LIST
        Lista de listas de listas que representa un polinomio de grado k con 
        coeficientes en el cuerpo Fq^m
    k : INT
        Grado del polinomio g
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
    c y e, listas de listas que representan vectores columna cuyos elementos pertenecen
    al cuerpo Fq, con c la palabra corregida del código C y e el error
    '''
    l = len(r)
    e = [f_q.cero(p, f)] * l # donde guardaremos el vector error
    sind = sindrome(r,a,g,p,f,h) # síndrome de r
    
    if sind != f_q_m_pol.cero(p, f, h): # la palabra recibida no pertenece al código
        c = [f_q.cero(p, f)] * l; # donde guardaremos la palabra que sí pertenece al código
        
        n,o = gcd_truncado(g,sind,k,p,f,h) # los polinomios necesarios para 
                                           # el algoritmo de decodificación
        od = f_q_m_pol.deriv(o, p, f, h) # o'(x)
        
        for i in range(l): # obtenemos los términos de error
            
            if f_q_m_pol.eval_pol(o, a[i], p, f, h) == f_q_m.cero(p, f, h): # en el término i hay un error
                na = f_q_m_pol.eval_pol(n, a[i], p, f, h) # n(ai)
                oda = f_q_m_pol.eval_pol(od, a[i], p, f, h) # o'(ai)
                oda_inv = f_q_m.inv_mult(oda, p, f, h) # inverso multiplicativo de o'(ai)
                
                ei = f_q_m.mult(na, oda_inv, p, f, h)
                e[i] = ei[0]

        for i in range(l): # c = r - e
            c[i] = f_q.suma(r[i],f_q.inv_adit(e[i], p, f), p, f)
        
        return(c,e)
    
    
    else: # la palabra recibida pertenece al código
        c = copy_vec(r)
        return(c,e)


#  Funciones auxiliares:
def gcd_truncado(a,b,k,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        lista de listas de listas que representa un polinomio del anillo Fq^m[x]
        mónico irreducible de grado k
    b : LIST
        lista de listas de listas que representa un polinomio del anillo Fq^m[x]
        de grado < k, congruente mod(a(x)) con sumatorio en i de ri/(x-alpha_i)
    k : INT
        Grado del polinomio a
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
    g1 y t1, polinomios con coeficientes en el cuerpo Fq^m que permitirán, mediante 
    el Algoritmo de Decodificación, determinar el vector de error y así recuperar el mensaje original
    '''
    s0 = f_q_m_pol.uno(p, f, h)
    t0 = f_q_m_pol.cero(p, f, h)
    g0 = f_q_m_pol.cop_pol(a, p, f, h)
    
    s1 = f_q_m_pol.cero(p, f, h)
    t1 = f_q_m_pol.uno(p, f, h)
    g1 = f_q_m_pol.cop_pol(b, p, f, h)
    
    while (len(g1) - 1) >= (k/2):
        s1_c = f_q_m_pol.cop_pol(s1, p, f, h)
        t1_c = f_q_m_pol.cop_pol(t1, p, f, h)
        g1_c = f_q_m_pol.cop_pol(g1, p, f, h)
        
        # calculamos los nuevos g1, s1 y t1
        q,g1 = f_q_m_pol.div_pol(g0, g1, p, f, h)
        s1 = f_q_m_pol.suma(s0,  f_q_m_pol.inv_adit(f_q_m_pol.mult(q, s1, p, f, h), p, f, h), p, f, h)
        t1 = f_q_m_pol.suma(t0,  f_q_m_pol.inv_adit(f_q_m_pol.mult(q, t1, p, f, h), p, f, h), p, f, h)
        
        # asignamos los nuevos valores a los g0, s0 y t0
        g0 = g1_c
        t0 = t1_c
        s0 = s1_c
    
    # nos aseguramos de que sean coprimos
    gcd = f_q_m_pol.gcd(t1, g1, p, f, h)    
    if gcd != f_q_m_pol.uno(p, f, h):
        t1 = f_q_m_pol.div_pol(t1, gcd, p, f, h)
        g1 = f_q_m_pol.div_pol(g1, gcd, p, f, h)
        
    return(g1,t1)


def copy_vec(v):
    '''
    Parameters
    ----------
    v : LIST
        lista de listas que representa un vector columna cuyos elementos pertenecen
        al cuerpo Fq

    Returns
    -------
    cop, una lista de listas que representa una copia del vector v, de manera que
    cop y v no estén relacionadas
    '''
    long = len(v)
    cop = [0]*long
    
    for i in range(long):
        cop[i] = v[i] + []
    
    return(cop)
