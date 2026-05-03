# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Cuerpo de q elementos, Fq, con q = p^n 
p número primo >= 2
n número entero >= 1

Se define Fq como (Z/pZ)[x]/<f> con deg(f) = n, f irreducible mónico

Contiene las siguientes funciones: 
    cero(p,f)
    uno(p,f)
    suma(a,b,p,f)
    inv_adit(a,p,f)
    mult(a,b,p,f)
    inv_mult(a,p,f)
    pot(a,k,p,f)
    rand_elem(p,f)
"""
import z_pz_pol, z_pz

def cero(p,f):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[x]

    Returns
    -------
    El elemento neutro de la suma del cuerpo Fq
    '''
    return(z_pz_pol.cero(p))


def uno(p,f):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[x]
        
    Returns
    -------
    El elemento unidad del cuerpo Fp
    '''
    return(z_pz_pol.uno(p))


def suma(a,b,p,f):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado < n perteneciente al anillo (Z/pZ)[x], el cual 
        constituye un elemento del cuerpo Fq
    b : LIST
        polinomio de grado < n perteneciente al anillo (Z/pZ)[x], el cual 
        constituye un elemento del cuerpo Fq
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[x]

    Returns
    -------
    Lista de tamaño < n que representa el elemento de Fq resultado de sumar
    a y b.
    '''
    sol = z_pz_pol.suma(a,b,p)
    return(sol)


def inv_adit(a,p,f):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado < n perteneciente al anillo (Z/pZ)[x], el cual 
        constituye un elemento del cuerpo Fq
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[x]
        
    Returns
    -------
    Lista de tamaño < n que representa el inverso aditivo del elemento a 
    en el cuerpo Fq
    '''
    sol = z_pz_pol.inv_adit(a,p)
    return(sol)


def mult(a,b,p,f):
   '''
   Parameters
   ----------
   a : LIST
       polinomio de grado < n perteneciente al anillo (Z/pZ)[x], el cual 
       constituye un elemento del cuerpo Fq
   b : LIST
       polinomio de grado < n perteneciente al anillo (Z/pZ)[x], el cual 
       constituye un elemento del cuerpo Fq
   p : INT
       número primo >= 2
   f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[x]
        
   Returns
   -------
   Lista de tamaño < n que representa el elemento del cuerpo Fq resultado 
   de multiplicar a y b.
   '''
   aux = z_pz_pol.mult(a,b,p)
   sol = z_pz_pol.div_pol(aux, f, p)[1]
   return (sol)


def inv_mult(a,p,f):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado < n perteneciente al anillo (Z/pZ)[x], el cual 
        constituye un elemento del cuerpo Fq
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[x]

    Returns
    -------
    sol, lista de tamaño < n que representa el inverso multiplicativo del 
    elemento a en el cuerpo Fq
    '''
    sol = z_pz_pol.gcd_ext(f,a,p)[2]
    sol = z_pz_pol.div_pol(sol, f, p)[1] #debe pertenecer al cuerpo
    return(sol)


def pot(a,k,p,f):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado < n perteneciente al anillo (Z/pZ)[x], el cual 
        constituye un elemento del cuerpo Fq
    k : INT
        número entero
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[x]

    Returns
    -------
    Lista que representa el elemento del cuerpo Fq resultado de calcular a^k 
    '''
    sol = z_pz_pol.pot_mod(a,k,p,f)
    return(sol)


def rand_elem(p,f):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[x]

    Returns
    -------
    Lista que representa un elemento aleatorio del cuerpo Fq
    '''
    n = len(f) - 1 #grado de f
    sol = [z_pz.cero(p)]*n
    
    for i in range(n):
        sol[i] = z_pz.rand_elem(p)
            
    sol = z_pz_pol.vd_len(sol, p)        
    return(sol)
