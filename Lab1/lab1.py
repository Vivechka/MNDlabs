import random as rd
import numpy as np
import prettytable as pt
import time as t

start_time = t.time()
X1 = np.array([rd.randint(1, 20) for i in range(8)])
X2 = np.array([rd.randint(1, 20) for i in range(8)])
X3 = np.array([rd.randint(1, 20) for i in range(8)])

a0 = 2
a1 = 3
a2 = 7
a3 = 4

Y = np.array([a0+a1*X1[i]+ a2*X2[i]+a3*X3[i] for i in range(8)])
Imax = np.where(Y == max(Y))[0][0]


x01 = (max(X1)+min(X1))/2
x02 = (max(X2)+min(X2))/2
x03 = (max(X3)+min(X3))/2

dx1 = x01 - min(X1)
dx2 = x02 - min(X2)
dx3 = x03 - min(X3)

Xn1 = np.array([(X1[i]-x01)/dx1 for i in range(8)])
Xn2 = np.array([(X2[i]-x02)/dx2 for i in range(8)])
Xn3 = np.array([(X3[i]-x03)/dx3 for i in range(8)])

Yet = a0+a1*x01+ a2*x02+a3*x03

table = pt.PrettyTable()
table.field_names = ['№', 'x1', 'x2', 'x3', 'Y', ' ', 'XN1', 'XN2', 'XN3']

for i in range(8):
    table.add_row([i+1, X1[i], X2[i], X3[i], Y[i], ' ', round(Xn1[i], 3), round(Xn2[i], 3), round(Xn3[i], 3)])
table.add_row(['x0', x01, x02, x03, Yet, ' ', '-', '-', '-'])
table.add_row(['dx', dx1, dx2, dx3, '-', ' ', '-', '-', '-'])


print(table)
print('Точка плану, що задовольняє заданому критерію оптимальості - (', X1[Imax],';', X2[Imax], ';', X3[Imax],')' )

print("Час роботи програми: %s секунд" % (t.time() - start_time))