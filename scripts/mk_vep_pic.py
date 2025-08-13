import re, sys, ast, csv
import matplotlib.pyplot as plt

def mk_csv(js_line, output):
    # print(js_line)
    p = r"arrayToDataTable\((.*?)\)"
    data_str = re.search(p, js_line).group(1)
    data_list = ast.literal_eval(data_str)

    with open(output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter = '\t')
        for row in data_list:
            writer.writerow(row)

def custom_autopct(pct):
    return f'{pct:.1f}%' if pct >= 5 else ''

def mk_pie(js_line, output):
    p = r"arrayToDataTable\((.*?)\)"
    data_str = re.search(p, js_line).group(1)
    data_list = ast.literal_eval(data_str)

    labels = [row[0] for row in data_list[1:]]
    sizes = [row[1] for row in data_list[1:]]

    # threshold = 1  # 设置阈值
    # other_size = sum(s for s in sizes if s/sum(sizes)*100 < threshold)
    # sizes = [s for s in sizes if s/sum(sizes)*100 >= threshold] + [other_size]
    # print(labels, sizes)
    # labels = [l for i,l in enumerate(labels) if sizes[i]>=threshold] + ['Other']

    # 1. 按数值大小降序排序（核心步骤）
    sorted_data = sorted(zip(sizes, labels), reverse=True)
    sorted_sizes, sorted_labels = zip(*sorted_data)

    # 2. 创建带数值的图例标签[[5]]
    # legend_labels = [f"{label} ({size})" for size, label in sorted_data]
    legend_labels = sorted_labels

    fig, ax = plt.subplots(figsize=(12, 8))
    wedges, texts, autotexts = ax.pie(
        sorted_sizes,
        labels=None,  # 禁用直接标签
        # autopct='%1.1f%%',
        autopct=custom_autopct,
        startangle=90,
        pctdistance=0.85,
        wedgeprops={'edgecolor': 'w', 'linewidth': 1.5}
    )

    # 4. 设置百分比文本样式[[2]]
    for autotext in autotexts:
        if autotext.get_text():  # 只处理显示的百分比
            autotext.set_size(10)
            autotext.set_color('white')
            autotext.set_weight('bold')

    # 5. 添加右侧图例[[6]][[11]]
    ax.legend(
        wedges, 
        legend_labels,
        loc="center left",
        bbox_to_anchor=(1, 0.5),  # 定位到右侧
        frameon=False,
        fontsize=10
    )

    plt.savefig(output, dpi=300, bbox_inches='tight')

def mk_bar(js_line, output):
    print(js_line)
    p = r"arrayToDataTable\((.*?)\)"
    data_str = re.search(p, js_line).group(1)
    data_list = ast.literal_eval(data_str)

    chromosomes = [row[0] for row in data_list[1:]]
    counts = [row[1] for row in data_list[1:]]

    def chromosome_key(chrom):
        chrom_num = chrom.replace("chr", "")
        if chrom_num.isdigit():
            return int(chrom_num)
        else:
            # chrX -> 23, chrY -> 24
            return {"X": 23, "Y": 24}.get(chrom_num, 99)

    sorted_indices = sorted(range(len(chromosomes)), key=lambda i: chromosome_key(chromosomes[i]))
    chromosomes = [chromosomes[i] for i in sorted_indices]
    counts = [counts[i] for i in sorted_indices]

    
    plt.figure(figsize=(14, 7))
    bars = plt.bar(chromosomes, counts, color='#1f77b4')

    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.ylim(0, max(counts) * 1.15)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.savefig(output, dpi=300, bbox_inches='tight')

js_lines = open(sys.argv[1]).readlines()

mk_csv(js_lines[12], 'var_class.csv')
mk_csv(js_lines[21], 'cons_type_severe.csv')
mk_csv(js_lines[30], 'cons_type_all.csv')
mk_csv(js_lines[39], 'coding_cons_type.csv')

mk_pie(js_lines[12], 'var_class.png')
mk_pie(js_lines[21], 'cons_type_severe.png')
mk_pie(js_lines[30], 'cons_type_all.png')
mk_pie(js_lines[39], 'coding_cons_type.png')

mk_bar(js_lines[47], 'var_chrom.png')
