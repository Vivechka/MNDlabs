import random
from prettytable import PrettyTable
from functools import partial
from scipy.stats import f, t as ts
    
m = 3
n = 8

x1_min = 20
x1_max = 70

x2_min = -20
x2_max = 40

x3_min = 70
x3_max = 80

y_min = 200 + (x1_min + x2_min + x3_min) / 3
y_max = 200 + (x1_max + x2_max + x3_max) / 3

plan_matrix_normal = [[-1, -1, -1, 1, 1, 1, -1],
       [-1, -1, 1, 1, -1, -1, 1],
       [-1, 1, -1, -1, 1, -1, 1],
       [-1, 1, 1, -1, -1, 1, -1],
       [1, -1, -1, -1, -1, 1, 1],
       [1, -1, 1, -1, 1, -1, -1],
       [1, 1, -1, 1, -1, -1, -1],
       [1, 1, 1, 1, 1, 1, 1]]

plan_matrix = [[x1_min, x2_min, x3_min, x1_min * x2_min, x1_min * x3_min, x2_min * x3_min, x1_min * x2_min * x3_min],
              [x1_min, x2_min, x3_max, x1_min * x2_min, x1_min * x3_max, x2_min * x3_max, x1_min * x2_min * x3_max],
              [x1_min, x2_max, x3_min, x1_min * x2_max, x1_min * x3_min, x2_max * x3_min, x1_min * x2_max * x3_min],
              [x1_min, x2_max, x3_max, x1_min * x2_max, x1_min * x3_max, x2_max * x3_max, x1_min * x2_max * x3_max],
              [x1_max, x2_min, x3_min, x1_max * x2_min, x1_max * x3_min, x2_min * x3_min, x1_max * x2_min * x3_min],
              [x1_max, x2_min, x3_max, x1_max * x2_min, x1_max * x3_max, x2_min * x3_max, x1_max * x2_min * x3_max],
              [x1_max, x2_max, x3_min, x1_max * x2_max, x1_max * x3_min, x2_max * x3_min, x1_max * x2_max * x3_min],
              [x1_max, x2_max, x3_max, x1_max * x2_max, x1_max * x3_max, x2_max * x3_max, x1_max * x2_max * x3_max]]
f1 = 0
f2 = 0  
q = 0.05            
while True:
    y_matrix = [[random.randint(int(y_min), int(y_max)) for _ in range(m)] for _ in range(n)]

    average_y = [round(sum(i) / len(i), 3) for i in y_matrix]

    b0 = sum(average_y) / n
    b1 = sum([average_y[i] * plan_matrix_normal[i][0] for i in range(n)]) / n
    b2 = sum([average_y[i] * plan_matrix_normal[i][1] for i in range(n)]) / n
    b3 = sum([average_y[i] * plan_matrix_normal[i][2] for i in range(n)]) / n
    b12 = sum([average_y[i] * plan_matrix_normal[i][0] * plan_matrix_normal[i][1] for i in range(n)]) / n
    b13 = sum([average_y[i] * plan_matrix_normal[i][0] * plan_matrix_normal[i][2] for i in range(n)]) / n
    b23 = sum([average_y[i] * plan_matrix_normal[i][1] * plan_matrix_normal[i][2] for i in range(n)]) / n
    b123 = sum([average_y[i] * plan_matrix_normal[i][0] * plan_matrix_normal[i][1] * plan_matrix_normal[i][2] for i in range(n)]) / n

    s = [sum([(y_matrix[j][i] - average_y[i]) ** 2 for i in range(m)]) / m for j in range(n)]
    gp = max(s) / sum(s)
    f1 = m-1
    f2 = n 
    q1 = q / f1
    fisher_value = f.ppf(q=1 - q1, dfn=f2, dfd=(f1 - 1) * f2)
    gt = fisher_value / (fisher_value + f1 - 1)

    if gp > gt:
        m += 1
    else: break

d = 8

sb = sum(s) / n
s_beta_2 = sb / (n * m)
s_beta = s_beta_2 ** (1 / 2)

bb = [b0, b1, b2, b3, b12, b13, b23, b123]
t = [abs(bb[i]) / s_beta for i in range(n)]
student = partial(ts.ppf, q=1 - q)
f3 = f1*f2
tt = student(df=f3)
blist = [b0, b1, b2, b3, b12, b13, b23, b123]

for i in range(n):
    if t[i] < tt:
        blist[i] = 0
        d -= 1

y_reg = [b0 + b1 * plan_matrix[i][0] + b2 * plan_matrix[i][1] + b3 * plan_matrix[i][2] +
             b12 * plan_matrix[i][3] + b13 * plan_matrix[i][4] + b23 * plan_matrix[i][5] +
             b123 * plan_matrix[i][6] for i in range(n)]
sad = (m / (n - d)) * int(sum([(y_reg[i] - average_y[i]) ** 2 for i in range(n)]))
fp = sad / sb
f4 = n-d
fisher = partial(f.ppf, q=0.95)
ft = fisher(dfn=f4, dfd=f3)
if fp > ft:
    toPrint = 'неадекватно оригіналу при рівні значимості 0.05'
else:
    toPrint = 'адекватно оригіналу при рівні значимості 0.05'

table = PrettyTable()

headers_x = ['X{}'.format(i) for i in range(0, m+1)]
headers_x.extend(['X12', 'X13', 'X23', 'X123'])
headers_y = ['Y{}'.format(i) for i in range(1, m+1)]
headers_y.extend(['av_Y', 'S^2'])

table.field_names = [*headers_x, *headers_y]
x0 = [[1] for _ in range(n)]
for i in range(n):
    table.add_row([*x0[i], *plan_matrix_normal[i], *y_matrix[i], average_y[i], s[i]])

print(table)
print("Дисперсія однорідна за критерієм Кохрена")
print("Кількість значущих коефіцієнтів за критерієм Стьюдента: ", d)
print("За критерієм Фішера рівняння регресії ", toPrint)
print('Рівняння')
print('y = {} + {} * x1 + {} * x2 + {} * x3 + {} * x1x2 + {} * x1x3 + {} * x2x3 + {} * x1x2x3'
      .format(round(b0, 3), round(b1, 3), round(b2, 3),
      round(b3, 3), round(b12, 3), round(b13, 3),
      round(b23, 3), round(b123, 3)))