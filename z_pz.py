# -*- coding: utf-8 -*-
"""
Diana Gómez Moreno

Cuerpo de p elementos, con p>=2 número primo, Z/pZ

Contiene las siguientes funciones: 
    cero(p)
    uno(p)
    suma(a,b,p)
    inv_adit(a,p)
    mult(a,b,p)
    inv_mult(a,p)
    pot(a,k,p)
    rand_elem(p)
    
Funciones auxiliares: 
    gcd(a,b)
"""
import random, sympy

def cero(p):
    '''    
    Parameters
    ----------
    p : INT
        número primo >= 2

    Returns
    -------
    El elemento neutro de la suma del cuerpo Z/pZ
    '''
    return(0)


def uno(p):
    '''    
    Parameters
    ----------
    p : INT
        número primo >= 2

    Returns
    -------
    El elemento unidad del cuerpo Z/pZ
    '''
    return(1)


def suma(a,b,p):
    '''
    Parameters
    ----------
    a : INT
        elemento del cuerpo Z/Zp
    b : INT
        elemento del cuerpo Z/Zp
    p : INT
        número primo >= 2

    Returns
    -------
    El elemento del cuerpo Z/pZ resultado de sumar a y b
    '''
    return((a+b) % p)


def inv_adit(a,p):
    '''
    Parameters
    ----------
    a : INT
        elemento del cuerpo Z/Zp
    p : INT
        número primo >= 2

    Returns
    -------
    El inverso aditivo de a en el cuerpo Z/pZ
    '''
    return((-a) % p)


def mult(a,b,p):
    '''
    Parameters
    ----------
    a : INT
        elemento del cuerpo Z/Zp
    b : INT
        elemento del cuerpo Z/Zp
    p : INT
        número primo >= 2

    Returns
    -------
    El elemento del cuerpo Z/pZ resultado de multiplicar a y b
    '''
    return((a*b) % p)


def inv_mult(a,p):
    '''
    Parameters
    ----------
    a : INT
        elemento del cuerpo Z/Zp
    p : INT
        número primo >= 2

    Returns
    -------
    El inverso multiplicativo de a en el cuerpo Z/pZ
    '''
    inv = gcd(a,p)[1]
    return(inv % p)


def pot(a,k,p):
    '''
    Parameters
    ----------
    a : INT
        elemento del cuerpo Z/Zp
    k : INT
        número entero
    p : INT
        número primo >= 2

    Returns
    -------
    El elemento del cuerpo Z/pZ resultado de calcular a^k
    '''    
    if k==0:
        return(uno(p))
    
    elif k > 0:
        
        if k%2==0: #k es par
            x = pot(a,k//2,p)
            sol = mult(x,x,p)
            return(sol)
        
        else: #k es impar
            sol = mult(a,pot(a,k-1,p),p)
            return(sol)
    else: 
        return(pot(inv_mult(a,p),(-1)*k,p))


def rand_elem(p):
    '''
    Parameters
    ----------
    p : INT
        número primo >= 2

    Returns
    -------
    Un elemento aleatorio del cuerpo Z/pZ 
    '''
    a = random.randint(0,p-1)   
    return(a)


# Funciones auxiliares:
def gcd(a,b):
    '''
    Parameters
    ----------
    a : INT
        número entero positivo
    b : TYPE
        número entero positivo

    Returns
    -------
    Los números enteros g,x,y tal que gcd(a,b)=g, con 
    a*x + y*b = g
    '''
    x, y, g = sympy.gcdex(a,b)              
    return(g,x,y)






























