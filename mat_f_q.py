# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Contiene las siguientes funciones: 
    mat_h(g,a,p,f,h)
    mat_cuad_al(k,p,f)
    mult_mat(a,b,p,f)
    mat_inv(M,p,f)
    vect_a_mat(v)
    mat_a_vect(m)
"""
import f_q, mat_f_q_m, random

def mat_h(g,a,p,f,h):
    '''
    Parameters
    ----------
    g : LIST
        lista de listas de listas que representa un polinomio de grado t con 
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
    Lista de listas de listas que representa una matriz de orden mtxl, definida por filas,
    cuyos elementos pertenecen al cuerpo Fq.
    Es una de las expresiones de la matriz de control de paridad del código de Goppa C 
    definido con los elementos de 'a' y polinomio g(x) del anillo Fq^m[x]
    '''
    m = len(h) - 1  # grado de h
    t = len(g) - 1  # grado de g
    l = len(a) # número de ai's
    
    # H será de dimensión mt x l
    H = [0]*(m*t) # definimos su número de filas
    
    for i in range(m*t):
        H[i] = [f_q.cero(p, f)] * l # inicializamos cada fila de modo que contenga 
                                    # en cada entrada el elemento cero del cuerpo Fq
    Hm = mat_f_q_m.mat_h(g,a,p,f,h) #Hm tiene dimensión t x l

    for i in range(t): # cada fila de la matriz Hm
        for j in range(l): # cada elemento de la fila i de la matriz Hm
            for k in range(len(Hm[i][j])): # cada coeficiente del elemento [i][j] en la matriz Hm
                H[i*m+k][j] = f_q.suma(Hm[i][j][k], f_q.cero(p, f), p, f) # para que no estén las dos matrices 
                                                                          # H y Hm relacionadas
    return(H)


def mat_cuad_al_inv(k,p,f):
    '''
    Parameters
    ----------
    k : INT
        dimensión de la matriz cuadrada resultante
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ

    Returns
    -------
    M y M_inv, listas de listas de listas, que representan matrices de dimensión
    kxk de elementos aleatorios de Fq tales que M_inv es la inversa de M
    '''
    M = [0]*k # número de filas
    
    for i in range(k):
        M[i] = [0]*k
        for j in range(k):
            M[i][j] = f_q.rand_elem(p, f)
    
    inv = False
    while not inv:
        try:
            M_inv = mat_inv(M,p,f)
            inv = True
        except ValueError:
            fila = random.randint(0,k-1)
            col = random.randint(0,k-1)
            M[fila][col] = f_q.rand_elem(p, f)
            
    return(M,M_inv)


def mult_mat(a,b,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa una matriz de dimensión fa x ca 
        definida por filas, cuyos elementos pertenecen al cuerpo Fq
    b : LIST
        Lista de listas de listas que representa una matriz de dimensión ca x cb 
        definida por filas cuyos elementos pertenecen al cuerpo Fq
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ
        
    Returns
    -------
    sol, lista de listas de listas, que representa una matriz de dimensión 
    fa x cb cuyos elementos pertenecen al cuerpo Fq, resultado del producto 
    de matrices a*b
    '''
    fa = len(a) # número de filas de a
    ca = len(a[0]) # número de columnas de a
    cb = len(b[0]) # número de columnas de b
    sol = [0] * fa # definimos el número de filas de sol
    
    for i in range(fa):
        sol[i] = [f_q.cero(p, f) for j in range(cb)] # inicializamos cada fila 
                                                     # para que contenga en cada entrada el elemento 
                                                     # cero del cuerpo Fq
    for i in range(fa):
        for j in range(cb):
            for k in range(ca):
                prod = f_q.mult(a[i][k], b[k][j], p, f)
                sol[i][j] = f_q.suma(sol[i][j], prod, p, f)
    
    return(sol)


def mat_inv(M,p,f):
    '''
    Parameters
    ----------
    M : LIST
        Lista de listas de listas que representa una matriz cuadrada, 
        definida por filas, cuyos elementos pertenecen al cuerpo Fq
    p : INT
        número primo >= 2
    f : LIST
        polinomio de grado n >= 1 irreducible, mónico, con coeficientes en
        el cuerpo Z/pZ

    Returns
    -------
    M_inv, lista de listas de listas que representa la matriz inversa de m, 
    definida por filas
    '''
    n = len(M) # M es una matriz de dimensiones nxn
    M_i = [0]*n # adjuntamos a M la matriz identidad nxn
                # luego, M_i tendrá dimensiones nx(2n)
    
    unidad = n 
    for i in range(n):
        M_i[i] = []
        for j in range(2*n):
            if j<n:
                M_i[i].append(M[i][j]+[]) 
            else:
                if j==unidad:
                    M_i[i].append(f_q.uno(p, f))
                else:
                    M_i[i].append(f_q.cero(p, f))
        unidad+=1 
    # Aplicamos reducción gaussiana a M_i:
    l = 2*n
    fi = 0 # contado = fila de M en la que estamos
    continuar = True # valdrá False cuando terminemos
    while continuar: 
        if fi >= n: # tiene rango máximo y hemos acabado
            continuar = False
            fi = fi - 1
            
        elif M_i[fi][fi] == f_q.cero(p, f): # el elemento (fi,fi) de la 
                                            # diagonal es el cero del cuerpo f_q                             
            # Vemos con un bucle si algún elemento de esa columna es dinstinto de cero
            # y precisamos hacer un intercambio de filas
            change = False # nos indicará si hay que hacer intercambio de filas
            r = fi + 1 
            while (r<n and (not change)):
                if M_i[r][fi] != f_q.cero(p,f):
                    change = True
                else:
                    r = r + 1
                
            if change: # Continuamos con la reducción gaussiana 
                       # intercambiando las filas fi y r
                fi_orig = M_i[fi]
                fr_orig = M_i[r]
                
                M_i[fi] = fr_orig
                M_i[r] = fi_orig
                # Realizamos la siguiente resta para continuar 
                # la reducción gaussiana por la misma fila fi, que
                # es la antigua fila r, en la siguiente vuelta del bucle
                fi = fi - 1

            else: # No hay candidado a pivote: la matriz no tiene inversa
                raise ValueError('No es invertible')
                continuar = False # la matriz está en forma triangular y 
                                  # salimos de este primer bucle

        else: 
            inv1 = f_q.inv_mult(M_i[fi][fi], p, f)
            for i in range(fi,l): # multiplicamos la fila fi por el inverso 
                                  # del elemento M_i[fi][fi] para que 
                                  # su primer elemento sea = f_q.uno()
                M_i[fi][i] = f_q.mult(M_i[fi][i], inv1, p, f)
               
            for i in range(fi+1,n): # eliminamos todos los elementos en 
                                    # la columna fi mediante resta de filas
                if M_i[i][fi] != f_q.cero(p, f):
                    for j in range(l-1,fi-1,-1):
                        mult = f_q.mult(M_i[fi][j], M_i[i][fi], p, f)
                        inv_adit = f_q.inv_adit(mult, p, f)
                        M_i[i][j] = f_q.suma(M_i[i][j], inv_adit, p, f)
                
        # Seguimos con la reducción gaussiana:
        fi = fi + 1

    for ix in range(n): # Tenemos solo la matriz en forma escalonada
                        # y queremos su forma escalonada reducida      
        i = n-1-ix   
        for j in range(i): # todas las filas a las que les deberemos restar la fila i
            for k in range(l-1,i-1,-1):
                mult = f_q.mult(M_i[i][k],M_i[j][i], p, f)
                inv_adit = f_q.inv_adit(mult, p, f)
                M_i[j][k] = f_q.suma(M_i[j][k], inv_adit, p, f) 
          
    M_inv = [M_i[i][n:] for i in range(n)]    
    return M_inv


def vect_a_mat(v):
    '''
    Parameters
    ----------
    v : LIST
        lista de listas que representa un vector columna cuyos elementos 
        pertenecen al cuerpo Fq

    Returns
    -------
    mat, una lista de listas de listas que representa una matriz con una única columna,
    cuyos elementos son los del vector v
    '''
    l = len(v)
    mat = [0] * l # donde guardaremos el vector en forma de matriz
    
    for i in range(l):
        mat_i = v[i]+[] # copiamos el elemento para que m y v no estén relacionadas
        mat[i] = [mat_i]
    
    return(mat)


def mat_a_vect(m):
    '''
    Parameters
    ----------
    m : LIST
        lista de listas de listas que representa una matriz columna cuyos elementos
        pertenecen al cuerpo Fq
    Returns
    -------
    v, una lista de listas que representa un vector columna cuyos elementos son los
    de la matriz m
    '''
    l = len(m)
    v = [0]*l
    for i in range(l):
        v_i = m[i][0]+[] # copiamos el elemento para que m y v no estén relacionadas
        v[i] = v_i
 
    return(v)
