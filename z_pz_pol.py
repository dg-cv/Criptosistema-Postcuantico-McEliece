# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Anillo de polinomios con coeficientes en el cuerpo Z/pZ: Z/pZ[x]
p número primo >= 2

Contiene las siguientes funciones: 
    cero(p)
    uno(p)
    suma(a,b,p)
    inv_adit(a,b,p)
    mult(a,b,p)
    gcd(a,b,p)
    gcd_ext(a,b,p)
    pot_mod(a,k,p,g) 

Funciones auxiliares: 
    vd_len(a,p)
    div_pol(f,g,p)
    gcd_ext_no_mon(a,b,p)
"""
import z_pz

def cero(p):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2

    Returns
    -------
    El elemento neutro de la suma del anillo Z/pZ[x]
    '''
    return([])


def uno(p):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2

    Returns
    -------
    El elemento unidad del anillo Z/pZ[x]
    '''

    return([z_pz.uno(p)])


def suma(a,b,p):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    b : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    p : INT
        número primo >= 2

    Returns
    -------
    Lista que representa el polinomio del anillo (Z/pZ)[x] resultado de sumar
    a y b
    '''
    la = len(a)
    lb = len(b)
    
    if la>=lb:
        
        sol = a + []
        
        for i in range(lb):
            sol[i] = z_pz.suma(a[i],b[i],p)
        
        sol = vd_len(sol,p)
        
    else:
        sol = suma(b,a,p)
    
    return(sol)


def inv_adit(a,p):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    p : INT
        número primo >= 2

    Returns
    -------
    Lista que representa el inverso aditivo del polinomio a del anillo (Z/pZ)[x]
    '''
    la = len(a)
    sol = [0]*la
    
    for i in range(la):
        sol[i] = z_pz.inv_adit(a[i],p)
    
    return(sol)


def mult(a,b,p):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    b : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    p : INT
        número primo >= 2

    Returns
    -------
    Lista que representa el resultado de multiplicar los polinomios a y b
    del anillo (Z/pZ)[x]
    '''   
    if a==cero(p) or b==cero(p):
        return(cero(p))
    
    else:
        la = len(a)
        lb = len(b)
        l = la + lb - 1
        sol = l*[z_pz.cero(p)]
    
        for i in range(la):
            for j in range(lb):
                sol[i+j] = z_pz.suma(sol[i+j],z_pz.mult(a[i],b[j],p),p)
        
        return(sol)


def gcd(a,b,p):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    b : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    p : INT
        número primo >= 2
        
    Returns
    -------
    g lista que representa el polinomio mónico del anillo (Z/pZ)[x] tal que 
    gcd(a,b) = g
    '''
    return(gcd_ext(a,b,p)[0])


def gcd_ext(a,b,p):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    b : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x], pero tal que
        deg(b)<=deg(a)
    p : INT
        número primo >= 2
        
    Returns
    -------
    g,x,y listas que representan los polinomio del anillo (Z/pZ)[x] tal que 
    gcd(a,b) = g, con a*x + b*y = g, con g mónico
    '''
    g,x,y = gcd_ext_no_mon(a,b,p)
    
    if g[-1] != z_pz.uno(p):
        coef = [z_pz.inv_mult(g[-1],p)]
        g = mult(g,coef,p)
        x = mult(x,coef,p)
        y = mult(y,coef,p)
        
    return(g,x,y)


def pot_mod(a,k,p,g):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    k : INT
        número entero, exponente de la potencia
    p : INT
        número primo >= 2
    g : LIST
        polinomio mónico irreducible de grado arbitrario del anillo (Z/pZ)[x]

    Returns
    -------
    Lista que representa el polinomio con coeficientes en (Z/pZ) a^k (mod g)
    '''   
    if (k==0):
        return(uno(p))
    
    elif (k>0):
        if k%2==0: #k es par
            x = pot_mod(a,k//2,p,g)
            aux = mult(x,x,p)
            sol = div_pol(aux,g,p)[1]
        
        else: #k es impar
            x = pot_mod(a,k-1,p,g)
            aux = mult(a,x,p)
            sol = div_pol(aux,g,p)[1]
        
        return(sol)
        
    else: #solo podremos calcular esta potencia si existe el inverso multiplicativo de a
          #modulo g
        aux = gcd_ext(a,g,p)[1]
        inv_mult = div_pol(aux,g,p)[1]
        
        return(pot_mod(inv_mult,-k,p,g))


# Funciones auxiliares:
def vd_len(a,p):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    p : INT
        número primo >= 2

    Returns
    -------
    Lista que representa un polinomio de grado arbitrario del anillo (Z/pZ)[x]
    de longitud únicamente hasta su mayor coeficiente no nulo   
    '''
    if len(a) == 0:
        return (a)
    
    else:
        while (len(a)>0 and a[-1]==z_pz.cero(p)):
            a.pop()

        return(a)


def div_pol(f,g,p):
    '''
    Parameters
    ----------
    f : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    g : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    p : INT
        número primo >= 2
        
    Returns
    -------
    q y r, listas que representan los polinomios del anillo (Z/pZ)[x] tal que
    f = q*g + r
    '''
    if (f == cero(p) or len(f)<len(g)):
        q = cero(p)
        r = f + []
        
        return(q,r)
    
    else:
        lq = len(f) - len(g) + 1 #longitud del coeficiente
        
        q = [z_pz.cero(p)]*(lq) 
        
        cg_i = z_pz.inv_mult(g[-1],p) # inverso del coeficiente principal de g
        
        while len(f)>=len(g): #mientras podamos dividir f entre g
            cf = f[-1] # coeficiente principal de f
            
            paso = [z_pz.cero(p)]*lq
            paso[len(f)-len(g)] = z_pz.mult(cf,cg_i,p)
            q = suma(q,paso,p)
            
            sumar = inv_adit(mult(g,paso,p),p)
            f = suma(f,sumar,p)
                   
        r = f + []

        return(q,r)
        
        
def gcd_ext_no_mon(a,b,p):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    b : LIST
        polinomio de grado arbitrario del anillo (Z/pZ)[x]
    p : INT
        número primo >= 2
        
    Returns
    -------
    g,x,y listas que representan a los polinomio del anillo (Z/pZ)[x] tal que 
    gcd(a,b) = g, con a*x + b*y = g, con g no necesariamente mónico
    '''
    if b == cero(p):
        return(a+[],uno(p),cero(p))
    
    else: 
        q,r = div_pol(a,b,p)
        
        g,a2,b2 = gcd_ext_no_mon(b,r,p)
        
        return(g, b2, suma(a2,inv_adit(mult(q,b2,p),p),p))

