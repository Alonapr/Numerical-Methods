# find the largest root x^3 - 4x^2 - 4x + 13 = 0

import numpy as np
from matplotlib import pyplot as plt
import math

def visualize_function():
    x = [-6, 0, 8]
    y = [0, 0, 0]
    plt.plot(x, y)
    x = np.linspace(-4, 6, 48)
    plt.plot(x, x**3 - 4*x*x - 4*x + 13)
    plt.show()
visualize_function()

# the interval with the largest root [4, 5]

class Resolve:
    def __init__(self):
        self.eps = 10**(-4)
        self.a1 = 4
        self.b1 = 5

    def func(self, x):
        return x**3 - 4*x*x - 4*x + 13

    def deriv(self, x):
        return 3*x*x - 8*x - 4

    def deriv2(self, x):
        return 6*x - 8

    def Relaxation_iterative_method(self, f, a, b, eps):
        print("\n" + "Relaxation iterative method:" + f"Interval: [{a}, {b}] ")

        m1 = min(abs(self.deriv(a)), abs(self.deriv(b)))
        M1 = max(abs(self.deriv(a)), abs(self.deriv(b)))
        tau = 2 / (M1 + m1)

        print("Verification of the sufficient convergence condition")
        print('===========================')
        if tau <= 0 or tau >= 2/M1:
            print("Narrow the interval.")
        else:
            print("The iterative process converges.")

            iteration_n = 1
            x0 = self.a1
            x_next = x0 - tau * f(x0)
            x_prev = x0
            print(f'Iteration {iteration_n}: x_{iteration_n} = {x_prev:.4f}')
            print('===========================')
            iteration_n += 1
            while abs(x_next - x_prev) > eps:
                print(f'Iteration {iteration_n}: x_{iteration_n} = {x_next:.4f}')
                print('===========================')
                x_prev = x_next
                x_next = x_prev - tau * f(x_prev)
                iteration_n += 1
            print(f'RESULT: {x_next}' + '\n')

            print(f'A POSTERIORY estimate of the number of iterations: {iteration_n - 1}')

            q = (M1 - m1) / (M1 + m1)
            n = int(math.log(((b - a) / eps), math.e) / math.log((1 / q), math.e)) + 1
            print(f'A PRIORI estimate of the number of iterations: {n}')

            return x_next



    def Modified_Newton_method(self, f, a, b, eps):
        print('\n' + '***************************')
        print("Modified Newton's method:" + f"Interval [{a}, {b}] ")

        print("Verification of the sufficient convergence condition")
        print('===========================')
        if self.deriv(a)*self.deriv(b) > 0 and self.deriv2(a)*self.deriv2(b) > 0:
            print("The iterative process converges.")

            iteration_n = 1
            x0 = self.a1
            x_next = x0 - f(x0) / self.deriv(x0)
            x_prev = x0
            print(f'Iteration {iteration_n}: x_{iteration_n} = {x_prev:.4f}')
            print('===========================')
            iteration_n += 1
            while abs(x_next - x_prev) > eps:
                print(f'Iteration {iteration_n}: x_{iteration_n} = {x_next:.4f}')
                print('===========================')
                x_prev = x_next
                x_next = x_prev - f(x_prev) / self.deriv(x0)
                iteration_n += 1
            print(f'RESULT: {x_next}' + '\n')

            print(f'A POSTERIORY estimate of the number of iterations: {iteration_n - 1}')
            return x_next
        else:
            print("Narrow the interval.")

    def result(self):
        res1 = self.Relaxation_iterative_method(self.func, self.a1, self.b1, self.eps)
        print("Largest root with the 1st method: ", res1)
        res2 = self.Modified_Newton_method(self.func, self.a1, self.b1, self.eps)
        print("Largest root with the 2nd method: ", res2)

r = Resolve()
r.result()