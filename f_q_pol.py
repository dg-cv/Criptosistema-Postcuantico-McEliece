# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Anillo de polinomios con coeficientes en el cuerpo de q elementos, Fq, con q = p^n 
p número primo
n número entero >= 1

Se define Fq como (Z/pZ)[t]/<f> con deg(f) = n, f irreducible mónico

Contiene las siguientes funciones: 
    cero(p,f)
    uno(p,f)
    suma(a,b,p,f)
    inv_adit(a,p,f)
    mult(a,b,p,f)
    gcd(a,b,p,f)
    gcd_ext(a,b,p,f)
    pot_mod(a,k,p,f,g)

Funciones auxiliares: 
    vd_len(a,p,f)
    div_pol(a,b,p,f)
    gcd_ext_no_mon(a,b,p,f)
    cop_pol(a,p,f)
"""
import f_q

def cero(p,f):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]

    Returns
    -------
    El elemento neutro de la suma del anillo Fq[x]
    '''
    return([])


def uno(p,f):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
        
    Returns
    -------
    El elemento unidad del anillo Fq[x]
    '''
    return([f_q.uno(p,f)])


def suma(a,b,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    b : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]

    Returns
    -------
    Lista de listas que representa el polinomio del anillo Fq[x]
    resultado de sumar a y b
    '''
    la = len(a)
    lb = len(b)
    
    if la>=lb:
        sol = [0]*la
        for i in range(la):
            if i<lb:
                sol[i] = f_q.suma(a[i],b[i],p,f)
            else: 
                sol[i] = a[i] + []
                
        sol = vd_len(sol,p,f)
        
    else: 
        sol = suma(b,a,p,f)
    
    return(sol)


def inv_adit(a,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]

    Returns
    -------
    Lista de listas que representa el inverso aditivo del elemento a del
    anillo de polinomios Fq[x]
    '''
    la = len(a)
    sol = [0]*la
    
    for i in range(la):
        sol[i] = f_q.inv_adit(a[i],p,f)
    
    return(sol)


def mult(a,b,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    b : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]

    Returns
    -------
    Lista de listas que representa el elemento del anillo de polinomios Fq[x]
    resultado de multiplicar a y b
    '''
    if a==cero(p,f) or b==cero(p,f):
        return(cero(p,f))
    
    else:
        la = len(a)
        lb = len(b)
        l = la + lb - 1
        sol = l*[f_q.cero(p,f)]
    
        for i in range(la):
            for j in range(lb):
                sol[i+j] = f_q.suma(sol[i+j],f_q.mult(a[i],b[j],p,f),p,f)
        
        return(sol)


def gcd(a,b,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    b : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
        
    Returns
    -------
    g lista de listas que representa el polinomio del anillo Fq[x] tal que 
    gcd(a,b) = g, con g mónico
    '''
    return(gcd_ext(a,b,p,f)[0])


def gcd_ext(a,b,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x] 
    b : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
        
    Returns
    -------
    g,x,y listas de listas que representan los polinomio del anillo Fq[x] 
    tal que gcd(a,b) = g, con g mónico, tal que a*x + b*y = g
    '''
    g,x,y = gcd_ext_no_mon(a,b,p,f)
    
    if g[-1] != f_q.uno(p,f):
        coef = [f_q.inv_mult(g[-1],p,f)]
        g = mult(g,coef,p,f)
        x = mult(x,coef,p,f)
        y = mult(y,coef,p,f)
        
    return(g,x,y)


def pot_mod(a,k,p,f,g):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    k : INT
        número entero. Es el exponente de la potencia
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    g : LIST
        lista de listas que representa un polinomio mónico irreducible de grado
        arbitrario del anillo Fq[x]

    Returns
    -------
    Lista de listas que representa el polinomio con coeficientes en el 
    cuerpo Fq a^k (mod g)
    '''
    if (k==0):
        return(uno(p,f))
    
    elif (k>0):
        if k%2==0: #k es par
            x = pot_mod(a,k//2,p,f,g)
            aux = mult(x,x,p,f)
            sol = div_pol(aux,g,p,f)[1]
        
        else: #k es impar
            x = pot_mod(a,k-1,p,f,g)
            aux = mult(a,x,p,f)
            sol = div_pol(aux,g,p,f)[1]
        
        return(sol)
        
    else: # solo podremos calcular esta potencia si existe el inverso multiplicativo de a
          # modulo g
        aux = gcd_ext(a, g, p, f)[1]
        inv_mul = div_pol(aux, g, p, f)[1]
        return(pot_mod(inv_mul,-k,p,f,g))


# Funciones auxiliares:
def vd_len(a,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]

    Returns
    -------
    Lista de listas que representa el mismo polinomio a, únicamente hasta su 
    coeficiente de mayor grado no nulo 
    '''
    if len(a) == 0:
        return(cero(p,f))
    else:
        while(len(a)>0 and a[-1]==f_q.cero(p,f)):
            a.pop()
            
        return(a)    
    

def div_pol(a,b,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    b : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]

    Returns
    -------
    q y r, listas de listas que representan los polinomios del anillo Fq[x] 
    tales que a = q*b + r
    '''
    if (a == cero(p,f) or len(a)<len(b)):
        q = cero(p,f)
        r = cop_pol(a,p,f) # es necesario copiar a, pues es el resto
        return(q,r)
    
    else:
        lq = len(a) - len(b) + 1 #longitud del cociente
        
        q = [f_q.cero(p,f)]*(lq) 
        
        cb_i = f_q.inv_mult(b[-1],p,f) # inverso del coeficiente principal de b
        
        while len(a)>=len(b): #mientras podamos dividir a entre b
            ca = a[-1] # coeficiente principal de a
            
            paso = [f_q.cero(p,f)]*lq
            paso[len(a)-len(b)] = f_q.mult(ca,cb_i,p,f)
            q = suma(q,paso,p,f)
            
            sumar = inv_adit(mult(b,paso,p,f),p,f)
            a = suma(a,sumar,p,f)
                   
        r = a # no me hace falta copiar a, porque a ya la hemos cambiado en el bucle, 
              # de modo que la variable a ya no está asociada a la lista de listas original
        return(q,r)        


def gcd_ext_no_mon(a,b,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    b : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
        
    Returns
    -------
    g,x,y listas de listas que representan a los polinomio del anillo Fq[x] 
    tal que gcd(a,b) = g, con a*x + b*y = g, g no necesariamente mónico
    '''
    if b == cero(p,f):
        return(cop_pol(a,p,f),uno(p,f),cero(p,f))
    
    else: 
        q,r = div_pol(a,b,p,f)
        g,a2,b2 = gcd_ext_no_mon(b,r,p,f)
        return(g, b2, suma(a2,inv_adit(mult(q,b2,p,f),p,f),p,f))


def cop_pol(a,p,f):
    '''
    Parameters
    ----------
    a : LIST
        Lista de listas que representa un polinomio de grado arbitrario del
        anillo Fq[x]
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]

    Returns
    -------
    Lista de listas que es una copia del polinomio a, independiente del 
    original.
    '''
    la = len(a)
    sol = [0]*la
    
    for i in range(la):
        sol[i] = a[i] + []

    return(sol)
