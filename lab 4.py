import numpy as np
from matplotlib import pyplot as plt

class Resolve:
    def __init__(self):
        self.a = -7.5
        self.b = 7.5
        self.n = 15
        self.h = (self.b - self.a)/self.n
        self.x = -0.5

    def func(self, x):
        return 1 / (x**2 + 1)

    def x_values(self):
        x = []
        x_prev = self.a
        for k in range(0, self.n + 1):
            x_next = x_prev + self.h
            x.insert(k, x_prev)
            x_prev = x_next
        return x

    def y_values(self):
        y = []
        for k in range(0, self.n + 1):
            y.insert(k, self.func(self.x_values()[k]))
        return y

    def Lagrange_Polinom(self, x_values, y_values, x):
        L_x = 0
        for k in range(0, self.n + 1):
            basic_pol = 1
            for j in range(0, self.n + 1):
                if j != k:
                    basic_pol *= (x - x_values[j])/(x_values[k] - x_values[j])
            L_x += basic_pol * y_values[k]
        return L_x

    def table(self, x_values, y_values):
        """
        :return: елементи таблиці зі значеннями функції у вузлах і розділеними різницями
        """
        k = [[0] * (self.n + 1) for i in range(self.n + 1)]
        for i in range(0, self.n + 1):
            k[i][0] = y_values[i]
        for i in range(1, self.n + 1):
            for j in range(i, self.n + 1):
                k[j][i] = (k[j][i - 1] - k[j - 1][i - 1]) / (x_values[j] - x_values[j - i])
        return k

    def split_dif(self, table):
        """
        :return: діагональні елементи таблиці - розділені різниці, необхідні для застосування
        інтерполяційної формули Ньютона
        """
        el = []
        for i in range(0, self.n + 1):
            el.append(table[i][i])
        return el

    def Newton_Polinom(self, split_dif, x_values, x):
        P_x = split_dif[0]
        for i in range(1, (self.n + 1)):
            basic_pol = split_dif[i]
            for j in range(i):
                basic_pol *= (x - x_values[j])
            P_x += basic_pol
        return P_x

    def visualize_function(self):
        x = np.linspace(-7.5, 7.5)
        plt.plot(x, 1 / (x**2 + 1))
        plt.plot(self.x_values(), self.y_values(), color="red")
        plt.show()

    def result(self):
        print("Крок h = ", self.h, "\n")
        print("Інтерполяційний поліном Лагранжа")
        print("x_k: ", self.x_values())
        print("y_k: ", self.y_values())
        Lagrange = self.Lagrange_Polinom(self.x_values(), self.y_values(), self.x)
        print("Наближене значення в точці -0.5: ", Lagrange, "\n")

        print("Інтерполяційний поліном Ньютона")
        #print(self.table(self.x_values(), self.y_values()))
        split_dif = self.split_dif(self.table(self.x_values(), self.y_values()))
        print("Необхідні розділені різниці: ", split_dif)
        Newton = self.Newton_Polinom(split_dif, self.x_values(), self.x)
        print("Наближене значення в точці -0.5: ", Newton)

r = Resolve()
r.result()
r.visualize_function()