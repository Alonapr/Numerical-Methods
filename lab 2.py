import numpy as np

A_gauss = np.array([[1, 5, 3, -4, 20], [3, 1, -2, 0, 9], [5, -7, 0, 10, -9], [0, 3, -5, 0, 1]])
A_g = np.array([[1, 5, 3, -4], [3, 1, -2, 0], [5, -7, 0, 10], [0, 3, -5, 0]])

def calc_solution(A):
    x = np.zeros(A[0].shape[0] - 1)
    for i in range(len(A)):
        x[-i - 1] = A[-i - 1, -1] - (x * A[-i - 1, :-1]).sum()
    return x

def Gauss_method(A):
    A_cur = A.copy()
    print("Початкова матриця для методу Гауса: ", "\n", A)
    iter = 1
    i_max = -1
    j_max = -1
    max_element = -9999
    n = 0
    while n < np.shape(A)[0]:
        for i in range(n, np.shape(A)[0]):
            for j in range(n, np.shape(A)[1]):
                if A[i][j] > max_element:
                    max_element = A[i][j]
                    i_max = i
                    j_max = j
            P = np.eye(A.shape[0])
            M = np.eye(A.shape[0])
            tmp = P[:, 0].copy()
            P[:, 0] = P[:, i_max].copy()
            P[:, i_max] = tmp
            print(f"Крок: {iter}")
            print("P: ", "\n", P)
            A_cur = P @ A_cur
            tmp = P[0, :].copy()
            P[0, :] = P[j_max, :].copy()
            P[j_max, :] = tmp
            A_cur = P @ A_cur
            for j in range(i, A.shape[0]):
                if i == j:
                    M[i][i] = 1 / A_cur[i, i]
                else:
                    M[j, i] = -A_cur[j, i] / A_cur[i, i]
            A_cur = np.round(M @ A_cur, 6)
            print("M: ", "\n", M)
            iter += 1
            n += 1
    print("Розв'язок системи: ", calc_solution(A_cur))

#Gauss_method(A_gauss)

A = np.array([[7, 2, 3, 1], [2, 5, 1, 2], [3, 1, 9, 4], [1, 2, 4, 8]])
b = np.array([1, 8, 6, 9])

def Seidel_method(A, b):
    print("\n")
    A = np.array(A, dtype='float32')
    b = np.array(b, dtype='float32')
    print("Початкова матриця для методу Зейделя: ", "\n", A)
    x = np.zeros_like(b)
    eps = 0.0001
    iter = 0
    stop_condition = False
    while not stop_condition:
        x_new = np.zeros_like(x)
        for i in range(len(A)):
            S1 = sum(A[i][j] * x_new[j] for j in range(i))
            S2 = sum(A[i][j] * x[j] for j in range(i+1, len(A)))
            x_new[i] = (b[i] - S1 - S2) / A[i][i]
        x_norm = max(abs(x_new[i] - x[i]) for i in range(len(A)))  # неперервна (кубічна) норма векторів
        stop_condition = x_norm <= eps
        if stop_condition:
            print("Умова припинення виконалась!")
            break
        x = x_new
        iter += 1
        print(f'Ітерація {iter}: x_{iter} = {x_new}')
        print('Норма векторів = ', x_norm)
        print('==========')
    print("Кількість ітерацій: ", iter)
    print(f"Розв'язок системи з точністю {eps}: ", x_new)

def Cond(A):
    print('\n')
    if np.linalg.det(A) != 0:
        A_inv = np.linalg.inv(A)
        print("Обернена матриця: ", "\n", A_inv)
        norma_A = np.linalg.norm(A, np.inf)              # рядок
        print("Норма початкової матриці = ", norma_A)
        norma_A_inv = np.linalg.norm(A_inv, np.inf)      # рядок
        print("Норма оберненої матриці = ", norma_A_inv)
        cond = norma_A * norma_A_inv
        print("Число обумовленості = ", cond)

Seidel_method(A, b)
Cond(A)


