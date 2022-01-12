import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


A = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
B = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
C = [0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100]
D = [0, 1, -2, 3, -4, 5, -6, 7, -8, 9, -10]


plt.plot(A, B, label='1')
plt.plot(A, C, label='2')
plt.plot(A, D, label='3')
plt.legend()
plt.show()
