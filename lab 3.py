import numpy as np

A = np.array([[7, 2, 3, 1], [2, 5, 1, 2], [3, 1, 9, 4], [1, 2, 4, 8]])

def Lamda_max(A):
    eps = 0.05
    m = 0
    iteration_n = 1
    x0 = np.array([4, 0, 0, 1])
    x_prev = x0
    x_next = A @ x0
    lam_next = x_next[m] / x_prev[m]
    print(f'Ітерація {iteration_n}: lamda_{iteration_n} = {lam_next:.4f}')
    print('===========================')
    iteration_n += 1
    x_prev = x_next
    x_next = A @ x_prev
    lam_prev = lam_next
    lam_next = x_next[m] / x_prev[m]
    print(f'Ітерація {iteration_n}: lamda_{iteration_n} = {lam_next:.4f}')
    iteration_n += 1
    while abs(lam_next - lam_prev) > eps:
        print('===========================')
        print(f'Ітерація {iteration_n}: lamda_{iteration_n} = {lam_next:.4f}')
        x_prev = x_next
        x_next = A @ x_prev
        lam_prev = lam_next
        lam_next = x_next[m] / x_prev[m]
        iteration_n += 1
    lam_max_A = lam_prev
    return lam_max_A

def Exponentiation_method(A):
    print("Знаходження максимального власного значення матриці A:", "\n", A)
    lam_max_A = Lamda_max(A)
    print("Максимальне власне значення: ", lam_max_A)
    B = np.array([[lam_max_A, 0, 0, 0], [0, lam_max_A, 0, 0], [0, 0, lam_max_A, 0], [0, 0, 0, lam_max_A]]) - A
    print("\n")
    print("Знаходження максимального власного значення матриці B:", "\n", B)
    lam_max_B = Lamda_max(B)
    lam_min_A = lam_max_A - lam_max_B
    print("Мінімальне власне значення: ", lam_min_A)

Exponentiation_method(A)