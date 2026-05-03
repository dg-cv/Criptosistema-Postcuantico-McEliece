# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Anillo de polinomios con coeficientes en el cuerpo de q^m elementos, Fq^m, con q = p^n 
p número primo
n número entero >= 1
m número entero >= 1

Se define Fq^m como Fq[s]/<h>, h irreducible mónico tal que deg(h) = m 

Contiene las siguientes funciones: 
    cero(p,f,h)
    uno(p,f,h)
    suma(a,b,p,f,h)
    inv_adit(a,p,f,h)
    mult(a,b,p,f,h)
    gcd(a,b,p,f,h)
    gcd_ext(a,b,p,f,h)
    pot_mod(a,k,p,f,g,h)
    eval_pol(a,u,p,f,h)
    deriv(a,p,f,h)

Funciones auxiliares: 
    vd_len(a,p,f,h)
    div_pol(a,b,p,f,h)
    gcd_ext_no_mon(a,b,p,f,h)
    cop_pol(a,p,f,h)
"""
import f_q_m, f_q_pol

def cero(p,f,h):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    El elemento neutro de la suma del anillo Fq^m[x]
    '''
    return([])


def uno(p,f,h):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    El elemento unidad del anillo Fq^m[x]
    '''
    return( [ f_q_m.uno(p,f,h) ] )


def suma(a,b,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    b : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    Lista de listas de listas que representa al polinomio del anillo Fq^m[x]
    resultado de sumar a y b
    '''
    la = len(a)
    lb = len(b)
    
    if la>=lb:
        sol = [0]*la
        for i in range(la):
            if i<lb:
                sol[i] = f_q_m.suma(a[i],b[i],p,f,h)
            else: 
                sol[i] = f_q_m.suma(a[i],f_q_m.cero(p,f,h),p,f,h)
                
        sol = vd_len(sol,p,f,h)
        
    else: 
        sol = suma(b,a,p,f,h)
    
    return(sol)


def inv_adit(a,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    Lista de listas de listas que representa el inverso aditivo del elemento a 
    del anillo Fq^m[x]
    '''
    la = len(a)
    sol = [0]*la
    
    for i in range(la):
        sol[i] = f_q_m.inv_adit(a[i],p,f,h)
    
    return(sol)


def mult(a,b,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    b : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    Lista de listas de listas que representa el elemento del anillo Fq^m[x] 
    resultado de multiplicar a y b
    '''
    if a==cero(p,f,h) or b==cero(p,f,h):
        return(cero(p,f,h))
    
    else:
        la = len(a)
        lb = len(b)
        l = la + lb - 1
        sol = l*[f_q_m.cero(p,f,h)]
    
        for i in range(la):
            for j in range(lb):
                sol[i+j] = f_q_m.suma(sol[i+j],f_q_m.mult(a[i],b[j],p,f,h),p,f,h)
        
        return(sol)


def gcd(a,b,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    b : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  
        
    Returns
    -------
    g lista de listas de listas que representa el polinomio del anillo Fq^m[x] 
    tal que gcd(a,b) = g, con g mónico
    '''
    return(gcd_ext(a,b,p,f,h)[0])


def gcd_ext(a,b,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    b : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  
        
    Returns
    -------
    g,x,y listas de listas de listas que representan los polinomio del anillo 
    Fq^m[x] tal que gcd(a,b) = g, con g mónico, y a*x + b*y = g
    '''
    g,x,y = gcd_ext_no_mon(a,b,p,f,h)
    
    if g[-1] != f_q_m.uno(p,f,h):
        coef = [ f_q_m.inv_mult(g[-1],p,f,h) ] #polinomio de Fq^m[x]
        g = mult(g,coef,p,f,h)
        x = mult(x,coef,p,f,h)
        y = mult(y,coef,p,f,h)
        
    return(g,x,y)


def pot_mod(a,k,p,f,h,g):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    b : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  
    g : LIST
        lista de listas de listas que representa un polinomio mónico irreducible de 
        grado arbitrario, perteneciente al anillo Fq^m[x]

    Returns
    -------
    Lista de listas de listas que representa el polinomio con coeficientes en el 
    cuerpo Fp^m a^k (mod g)
    '''
    if (k==0):
        return(uno(p,f,h))
    
    elif (k>0):
        if k%2==0: #k es par
            x = pot_mod(a,k//2,p,f,h,g)
            aux = mult(x,x,p,f,h)
            sol = div_pol(aux,g,p,f,h)[1]
        
        else: #k es impar
            x = pot_mod(a,k-1,p,f,h,g)
            aux = mult(a,x,p,f,h)
            sol = div_pol(aux,g,p,f,h)[1]
        
        return(sol)
        
    else: #solo podremos calcular esta potencia si existe el inverso multiplicativo de a
          #modulo g
        aux = gcd_ext(a, g, p, f, h)[1]
        inv_mul = div_pol(aux, g, p, f, h)[1]
        return(pot_mod(inv_mul,-k,p,f,h,g))


def eval_pol(a,u,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    u : LIST
        lista de listas que representa un elemento del cuerpo Fq^m
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    Lista de listas de tamaño < m que representa el elemento del cuerpo Fq^m que 
    se obtiene al evaluar el polinomio a en u
    '''
    sol = f_q_m.cero(p, f, h)
    la = len(a)
    
    for i in range(la):
        potc = f_q_m.pot(u, i, p, f, h)
        multip = f_q_m.mult(a[i], potc, p, f, h)
        sol = f_q_m.suma(sol, multip, p, f, h)
    
    return(sol)


def deriv(a,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    Lista de listas de listas que representa el elemento del anillo Fq^m[x] que 
    se obtiene al derivar el polinomio a
    '''
    if len(a) == 0:
       return cero(p,f,h)
    
    ad = [f_q_m.cero(p, f, h)]*(len(a)-1)
    
    for i in range(1,len(a)):
        ai = a[i]
        for j in range(i):
            ad[i-1] = f_q_m.suma(ad[i-1],ai,p,f,h)
            
    ad = vd_len(ad,p,f,h)
    return(ad)


# Funciones auxiliares:
def vd_len(a,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    Lista de listas de listas que representa el mismo polinomio a, únicamente hasta su 
    coeficiente de mayor grado no nulo 
    '''
    if len(a) == 0:
        return(cero(p,f,h))
    
    else:
        while(len(a)>0 and a[-1]==f_q_m.cero(p,f,h)):
            a.pop()
            
        return(a)    
    

def div_pol(a,b,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    b : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    q y r, listas de listas de listas que representan los polinomios del anillo 
    Fq^m[x] tales que a = q*b + r
    '''
    if (a == cero(p,f,h) or len(a)<len(b)):
        q = cero(p,f,h)
        r = cop_pol(a,p,f,h) # es necesario copiar a, pues es el resto
        return(q,r)
    
    else:
        lq = len(a) - len(b) + 1 #longitud del cociente
        q = [f_q_m.cero(p,f,h)]*lq 
        cb_i = f_q_m.inv_mult(b[-1],p,f,h) # inverso del coeficiente principal de b
        
        while len(a)>=len(b): #mientras podamos dividir a entre b
            ca = a[-1] # coeficiente principal de a
            
            paso = [f_q_m.cero(p,f,h)]*lq
            paso[len(a)-len(b)] = f_q_m.mult(ca,cb_i,p,f,h)
            q = suma(q,paso,p,f,h)
            
            sumar = inv_adit(mult(b,paso,p,f,h),p,f,h)
            a = suma(a,sumar,p,f,h)
                   
        r = a # no me hace falta copiar a, porque a ya la hemos cambiado en el bucle, 
              # de modo que la variable a ya no está asociada a la lista de listas de listas original
        return(q,r)
        
        
def gcd_ext_no_mon(a,b,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    b : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  
    
    Returns
    -------
    g,x,y listas de listas de listas que representan los polinomio del anillo 
    Fq^m[x] tales que gcd(a,b) = g, con a*x + b*y = g, g no necesariamente mónico
    '''
    if b == cero(p,f,h):
        return(cop_pol(a,p,f,h),uno(p,f,h),cero(p,f,h))
    
    else: 
        q,r = div_pol(a,b,p,f,h)
        g,a2,b2 = gcd_ext_no_mon(b,r,p,f,h)
        return(g, b2, suma(a2,inv_adit(mult(q,b2,p,f,h),p,f,h),p,f,h))


def cop_pol(a,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas de listas que representa un polinomio de grado arbitrario del
        anillo Fq^m[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : TYPE
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[s]  

    Returns
    -------
    Lista de listas de listas que es una copia del polinomio a, 
    independiente del original.
    '''
    la = len(a)
    sol = [0]*la
    
    for i in range(la):
        sol[i] = f_q_pol.cop_pol(a[i],p,f)
    return(sol)
