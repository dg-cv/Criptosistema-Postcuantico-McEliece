import mceliece, goppa, f_q, random, time

random.seed(14)

inicio = time.time()

# Ejemplo:
p = 2
f = [1,1]
q = 2
h = [[1], [], [], [1], [], [], [], [], [], [], [1]]
m = 10


g = [ [[1]] ] # irreducible en F2^10[x], t = 50

for i in range(1, 51):
    
    if i==35:
        g.append([[1], [1], [1], [], [], [1]])
    
    elif i==50:
        g.append([[1]])
    
    else:
        g.append([])

#a = goppa.gen_ai(g,600,p,f,h)
#a = goppa.gen_ai(g,700,p,f,h)
a = goppa.gen_ai(g,1024,p,f,h)

# Gen:
sk,pk = mceliece.gen(g, a, p, f, h)


# Enc: 
dim = len(sk[0][0]) # obtenemos la dimensión del código de Goppa resultante
print(f'La dimensión del código es {dim}')
r = [0]*dim
for i in range(dim):
    r[i] = f_q.rand_elem(p, f)
ct = mceliece.enc(r, pk[0], pk[1], p, f)

# Dec:
r_recuperado = mceliece.dec(ct,sk[1],sk[2],sk[4],sk[3],p,f,h)
print(f'El mensaje enviado era {r}')
print(f'El mensaje recuperado es {r_recuperado}\n')
print(f'¿Es r=r_recuperado? : {r==r_recuperado}\n')

# Tiempo de ejecución
fin = time.time()
print(f"Tiempo de ejecución: {fin - inicio:.4f} segundos")