import gzip, sys
import matplotlib.pyplot as plt
import numpy as np

# 1. 定义7个区域
REGIONS = [
    (-200000, -10000),
    (-10000, -1000),
    (-1000, -50),
    (-50, 50),
    (50, 1000),
    (1000, 10000),
    (10000, 200000)
]

# 2. 解析VCF.GZ
vcf_gz_path = sys.argv[1]
region_data = [[] for _ in range(7)]  # 7个空列表存储各区域长度值

with gzip.open(vcf_gz_path, 'rt') as f:
    for line in f:
        if line.startswith('#'): 
            continue
        parts = line.split('\t')
        ref_len = len(parts[3])
        alt_len = len(parts[4].split(',')[0])
        length = alt_len - ref_len
        
        for i, (low, high) in enumerate(REGIONS):
            if low <= length < high or (i == 3 and low <= length <= high):
                region_data[i].append(length)
                break

# 3. 创建图表
fig, axes = plt.subplots(1, 7, figsize=(21, 7), sharey=True)

TITLES = [
    "DEL≥10,000",
    "1,000<DEL≤10,000",
    "50≤DEL<1,000",
    "SNV and Indel",
    "50≤INS<1,000",
    "1,000≤INS/DUP<10,000",
    "INS/DUP≥10,000"
]

# 4. 为每个区域绘制详细分布
for i, ax in enumerate(axes):
    data = region_data[i]
    if data:
        # 计算每个长度值的频次
        unique, counts = np.unique(data, return_counts=True)
        
        # 绘制条形图（每个长度值一个条形）
        ax.bar(unique, counts, width=1, color='#1f77b4')
        
        # 设置区域边界
        low, high = REGIONS[i]
        ax.set_xlim(low, high)
    
    # 设置区域标题
    ax.set_title(f"{TITLES[i]}", fontsize=20)

    step = (high - low) // 2
    # 设置主要刻度
    ax.set_xticks(np.arange(low, high, step))
    
    # 禁用科学计数法
    ax.ticklabel_format(axis='x', style='plain', useOffset=False)
    
    # 旋转刻度标签45度
    ax.tick_params(axis='x', rotation=45)

# 5. 设置对数坐标轴
plt.yscale('log')
axes[0].set_ylabel("Count", fontsize = 20)
plt.yticks([10**i for i in range(0, 8)])

# 6. 输出PNG
plt.tight_layout()
plt.savefig("pangenie_var_plot.png", dpi=150)
plt.close()
