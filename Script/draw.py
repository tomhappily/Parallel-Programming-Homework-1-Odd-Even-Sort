import matplotlib.pyplot as plt
import numpy as np
import json


categories = [4, 12, 24, 48]

process4_total_time = []
process12_total_time = []
process24_total_time = []
process48_total_time = []

process4_io_time = []
process12_io_time = []
process24_io_time = []
process48_io_time = []

process4_comm_time = []
process12_comm_time = []
process24_comm_time = []
process48_comm_time = []

total_time = []
io_time = []
comm_time = []


for i in range(4):
    with open('./time/4process/rank'+str(i)+'.json', 'r') as file:
        data = json.load(file)
        process4_total_time.append(data["total_time"])
        process4_io_time.append(data["io_time"])
        process4_comm_time.append(data["comm_time"])

for i in range(12):
    with open('./time/12process/rank'+str(i)+'.json', 'r') as file:
        data = json.load(file)
        process12_total_time.append(data["total_time"])
        process12_io_time.append(data["io_time"])
        process12_comm_time.append(data["comm_time"])

for i in range(24):
    with open('./time/24process/rank'+str(i)+'.json', 'r') as file:
        data = json.load(file)
        process24_total_time.append(data["total_time"])
        process24_io_time.append(data["io_time"])
        process24_comm_time.append(data["comm_time"])

for i in range(48):
    with open('./time/48process/rank'+str(i)+'.json', 'r') as file:
        data = json.load(file)
        process48_total_time.append(data["total_time"])
        process48_io_time.append(data["io_time"])
        process48_comm_time.append(data["comm_time"])

total_time.append(max(process4_total_time))
total_time.append(max(process12_total_time))
total_time.append(max(process24_total_time))
total_time.append(max(process48_total_time))

speedUp = []

for i in range(4):
    speedUp.append(total_time[0]/total_time[i])

plt.plot(categories, speedUp, label='speedup', color='blue', linestyle='-', linewidth=2)
x_indexes = np.arange(len(categories))
# 设置标签和标题
plt.xlabel('processes')
plt.ylabel('runtime(seconds)')
plt.xticks(categories)
plt.legend()

# 显示图形
plt.show()

# io_time.append(max(process4_io_time))
# io_time.append(max(process12_io_time))
# io_time.append(max(process24_io_time))
# io_time.append(max(process48_io_time))
#
# comm_time.append(max(process4_comm_time))
# comm_time.append(max(process12_comm_time))
# comm_time.append(max(process24_comm_time))
# comm_time.append(max(process48_comm_time))
#
# cpu_time = [a - b - c for a, b, c in zip(total_time, io_time, comm_time)]
#
# # 设置柱的宽度
# bar_width = 0.5
# x_indexes = np.arange(len(categories))
# print(x_indexes)
#
# # 绘制堆叠柱状图
# plt.bar(x_indexes, cpu_time, width=bar_width, label='CPU', color='skyblue')
# plt.bar(x_indexes, io_time, width=bar_width, bottom=cpu_time, label='I/O', color='salmon')
# plt.bar(x_indexes, comm_time, width=bar_width, bottom=np.add(cpu_time, io_time), label='Comm', color='lightgreen')
#
# # 设置标签和标题
# plt.xlabel('processes')
# plt.ylabel('runtime(seconds)')
# plt.xticks(x_indexes, categories)
# plt.legend()
#
# # 显示图形
# plt.show()
