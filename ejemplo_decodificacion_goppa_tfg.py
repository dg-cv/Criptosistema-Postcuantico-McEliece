import goppa, f_q, random

random.seed(14)

# Función auxiliar:
def incluir_error(c, t, p, f):
    l = len(c)
    r = goppa.copy_vec(c)
    error = [f_q.cero(p, f)] * l  # vector error
    max_err = (t) // 2  # máximo número de errores
    peso = 0
    for i in range(l):
        rand = random.randint(0, 1)
        if (rand == 1 and peso < max_err):
            error[i] = f_q.rand_elem(p, f)
            r[i] = f_q.suma(c[i], error[i], p, f)
            peso += 1
        if peso == max_err:
            break
    return (r, peso, error)

# Ejemplo:
p = 911
q = 911**3
m = 2
f = [209, 840, 191, 1]
h = [[38, 906, 76], [832, 683, 311], [1]]
g = [[[2, 1]], [[2, 0, 1], [1]], [[], [1]], [[1]]]
t = len(g) - 1  # grado de g
a = [[[288, 512, 260], [104, 124, 837]], [[91, 285, 287], [121, 28, 165]],
     [[758, 421, 116], [702, 538, 600]], 
     [[95, 435, 491], [701, 172, 547]],
     [[395, 466, 316], [496, 470, 447]], 
     [[427, 609, 846], [97, 262, 485]],
     [[403, 241, 460], [886, 630, 500]]]

G = goppa.mat_g_can(g, a, p, f, h)
print(f'La matriz G es {G}\n')

c = [G[j][0] for j in range(len(G))]
print(f'Palabra transmitida {c}\n')

r, peso, error = incluir_error(c, t, p, f)
print(f'Error real {error}\n')
print(f'Palabra con error {r}\n')

sind = goppa.sindrome(r, a, g, p, f, h)
print(f'Síndrome de r {sind}\n')

n, o = goppa.gcd_truncado(g, sind, t, p, f, h)
print(f'sigma es {o}\n')
print(f'eta es {n}\n')

corr, err_decode = goppa.decode(r, a, g, t, p, f, h) # palabra corregida 
                                                     # y error
print(f'Error determinado con decode {err_decode}\n')
print(f'Palabra decodificada {corr}')

