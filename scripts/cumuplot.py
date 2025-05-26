import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

def parse_hapblock(file_path):
    spans = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('BLOCK'):
                parts = line.strip().split()
                span = float(parts[-3])
                spans.append(span / 1e6)  # 转换为Mbp
    return sorted(spans, reverse=True)  # 按长度降序排列

spans = parse_hapblock(sys.argv[1])

total_span = sum(spans)
cumulative_span = np.cumsum(spans)
cumulative_coverage = cumulative_span / total_span

n50_index = np.argmax(cumulative_coverage >= 0.5)
n50_value = spans[n50_index]


plt.figure(figsize=(10, 6))
plt.plot(cumulative_coverage, spans, color='#2c7bb6', linewidth=2, label='Phase Blocks')

# 标注N50
plt.axvline(x=0.5, color='gray', linestyle='--', linewidth=1.5, 
            label=f'N50 = {n50_value:.2f} Mbp')
plt.text(0.5, max(spans)*0.8, f'N50: {n50_value:.2f} Mbp', 
         verticalalignment='top', color='gray')

# 坐标轴美化
plt.xlabel('Cumulative Coverage', fontsize=12)
plt.ylabel('Phase Block Size (Mbp)', fontsize=12)
# plt.title('Cumulative Coverage vs. Phase Block Size (HapCut2)', fontsize=14)
# plt.grid(True, linestyle='--', alpha=0.7)
plt.xlim(0, 1)
plt.legend()
plt.tight_layout()
plt.savefig('cumulative_coverage_plot.png', dpi=300)

