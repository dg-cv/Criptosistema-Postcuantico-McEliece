import mceliece, random

random.seed(25)

# Ejemplo:
p = 2
f = [1,1]
q = 2
h = [[1], [1], [1], [1], [1]]
m = 4

g = [ [[],[1]], [[1]], [[1]], [[1]] ] # irreducible en F2^4[x], t = 3
a = [[], [[1]], [[1], [], [1], [1]], [[], [], [1], [1]],
    [[1], [], [], [1]], [[1], [1], [1]], [[1], [], [1]],
    [[1], [1]], [[1], [1], [], [1]], [[], [], [], [1]],
    [[], [1], [1], [1]], [[], [1], [], [1]],
    [[], [1], [1]], [[], [], [1]]]

# Gen:
sk,pk = mceliece.gen(g, a, p, f, h)

# Enc: 
r = [[1],[1]]
ct = mceliece.enc(r, pk[0], pk[1], p, f)

# Dec:
r_recuperado = mceliece.dec(ct,sk[1],sk[2],sk[4],sk[3],p,f,h)
print('\n')
print(f'El mensaje enviado era {r}')
print(f'El mensaje recuperado es {r_recuperado}')