# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Cuerpo de q^m elementos, Fq^m, con q = p^n 
p número primo >= 2
n número entero >= 1
m número entero >= 1

Se define Fq^m como Fq[x]/<h>, h irreducible mónico tal que deg(h) = m 

Contiene las siguientes funciones: 
    cero(p,f,h)
    uno(p,f,h)
    suma(a,b,p,f,h)
    inv_adit(a,p,f,h)
    mult(a,b,p,f,h)
    inv_mult(a,p,f,h)
    pot(a,k,p,f,h)
    rand_elem(p,f,h)   
"""
import f_q_pol, f_q

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
        de grado m >= 1, perteneciente al anillo Fq[x]  

    Returns
    -------
    El elemento neutro de la suma del cuerpo Fq^m
    '''
    return(f_q_pol.cero(p,f))


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
        de grado m >= 1, perteneciente al anillo Fq[x]  

    Returns
    -------
    El elemento unidad del cuerpo Fq^m
    '''
    return(f_q_pol.uno(p,f))


def suma(a,b,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado < m perteneciente al anillo Fq[x], el cual 
        constituye un elemento del cuerpo Fq^m
    b : LIST
        polinomio de grado < m perteneciente al anillo Fq[x], el cual 
        constituye un elemento del cuerpo Fq^m
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[x]  
        
    Returns
    -------
    Lista de listas de tamaño < m que representa el elemento de Fq^m 
    resultado de sumar a y b.
    '''
    sol = f_q_pol.suma(a,b,p,f)
    return(sol)


def inv_adit(a,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado < m perteneciente al anillo Fq[x], el cual 
        constituye un elemento del cuerpo Fq^m
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[x]  
    
    Returns
    -------
    Lista de listas de tamaño < m que representa el inverso aditivo del 
    elemento a del cuerpo Fq^m
    '''
    sol = f_q_pol.inv_adit(a,p,f)
    return(sol)


def mult(a,b,p,f,h):
   '''
   Parameters
   ----------
   a : LIST
       polinomio de grado < m perteneciente al anillo Fq[x], el cual 
       constituye un elemento del cuerpo Fq^m
   b : LIST
       polinomio de grado < m perteneciente al anillo Fq[x], el cual 
       constituye un elemento del cuerpo Fq^m
   p : INT
       número primo >= 2
   f : LIST
       polinomio mónico, irreducible, de grado n >= 1, perteneciente 
       al anillo (Z/pZ)[t]
   h : LIST
       lista de listas que representa un polinomio mónico, irreducible, 
       de grado m >= 1, perteneciente al anillo Fq[x]  

   Returns
   -------
   Lista de listas de tamaño < m que representa el elemento del cuerpo Fq^m 
   resultado de multiplicar a y b.
   '''
   aux = f_q_pol.mult(a,b,p,f)
   sol = f_q_pol.div_pol(aux, h, p, f)[1]
   return (sol)


def inv_mult(a,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado < m perteneciente al anillo Fq[x], el cual 
        constituye un elemento del cuerpo Fq^m
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[x]  

    Returns
    -------
    sol, lista de tamaño < m que representa el inverso multiplicativo del 
    elemento a del cuerpo Fq^m
    '''
    sol = f_q_pol.gcd_ext(h,a,p,f)[2]
    sol = f_q_pol.div_pol(sol, h, p, f)[1] #debe pertenecer al cuerpo
    return(sol)


def pot(a,k,p,f,h):
    '''
    Parameters
    ----------
    a : LIST
        polinomio de grado < m perteneciente al anillo Fq[x], el cual 
        constituye un elemento del cuerpo Fq^m
    k : INT
        número entero
    p : INT
        número primo >= 2
    f : LIST
        polinomio mónico, irreducible, de grado n >= 1, perteneciente 
        al anillo (Z/pZ)[t]
    h : LIST
        lista de listas que representa un polinomio mónico, irreducible, 
        de grado m >= 1, perteneciente al anillo Fq[x]  

    Returns
    -------
    lista de listas que representa al elemento del cuerpo Fq^m resultado de
    calcular a^k 
    '''
    sol = f_q_pol.pot_mod(a,k,p,f,h)
    return(sol)


def rand_elem(p,f,h):
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
        de grado m >= 1, perteneciente al anillo Fq[x]  

    Returns
    -------
    Una lista de listas que representa un elemento aleatorio del cuerpo Fq^m
    '''
    n = len(h) - 1 #grado de h
    sol = [f_q.cero(p,f)]*n
    
    for i in range(n):
        sol[i] = f_q.rand_elem(p,f)
            
    sol = f_q_pol.vd_len(sol, p, f)        

    return(sol)
