from matplotlib import pyplot as plt
from scipy import integrate
import numpy as np

a = 1
b = 5
eps = 0.005

def func(x):
    return 1/(x+2)

v, err = integrate.quad(func, a, b)
print("Точний розв'язок інтеграла: ", v, "\n")

def visualize_function():
    x = np.linspace(1, 5)
    plt.plot(x, 1/(2+x))
    plt.show()

visualize_function()

def Simpson_Method_residual_terms():
    print("Використовуючи оцінку залишкових членів")
    h = 1
    n = 4
    nodes = []
    for i in range(n + 1):
        nodes.append(a + i * h)
    print("Рівновіддалені вузли з кроком", h, ":", nodes)
    p1 = 0
    for i in range(1, int(n/2 + 1)):
        p1 += func(nodes[2 * i - 1])
    p2 = 0
    for i in range(1, int(n/2)):
        p2 += func(nodes[2 * i])
    result = (h / 3) * (func(nodes[0]) + 4 * p1 + 2 * p2 + func(nodes[len(nodes) - 1]))
    print("I = ", result)

Simpson_Method_residual_terms()

def Simpson_Method_Runge():
    print("\nВикористовуючи правило Рунге")
    n = 2
    h1 = (b-a)/n
    nodes1 = []
    for i in range(n+1):
        nodes1.append(a + i * h1)
    print("Рівновіддалені вузли з кроком", h1, ":", nodes1)
    result = (h1 / 3) * (func(nodes1[0]) + 4 * func(nodes1[1]) + func(nodes1[len(nodes1) - 1]))
    print("I1 = ", result)

    iter = 1
    result_prev = result
    n_prev = n
    posteriori_eval = False
    while not posteriori_eval:
        n_new = 2*n_prev
        h = (b-a)/n_new
        nodes = []
        for i in range(n_new + 1):
            nodes.append(a + i * h)
        print("Рівновіддалені вузли з кроком", h, ":", nodes)
        p1 = 0
        for i in range(1, n_prev + 1):
            p1 += func(nodes[2 * i - 1])
        p2 = 0
        for i in range(1, n_prev):
            p2 += func(nodes[2 * i])
        result_next = (h / 3) * (func(nodes[0]) + 4 * p1 + 2 * p2 + func(nodes[len(nodes) - 1]))
        iter += 1
        print(f"I{iter} = ", result_next)
        posteriori_eval = abs(result_prev - result_next)/(2**4 - 1) <= eps
        if posteriori_eval:
            break
        result_next = result_prev
        n_new = n_prev
    print("Відповідь: ", result_next, "з кроком", int(h))

Simpson_Method_Runge()