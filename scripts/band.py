import sys, re

hapblock = sys.argv[1]
MIN_SPAN = 5000  # 最小SPAN阈值

def hapcut2_to_ideogram(input_file, output_file):
    # 初始化变量
    tmp_blocks = []
    x_blocks = []

    current_block = {"chr": None, "positions": [], "value": 1, "span": 0}
    block_counter = 0  # 用于交替颜色
    
    with open(input_file, 'r') as f_in:        
        for line in f_in:
            line = line.strip()
            
            # 检测BLOCK起始行
            if line.startswith("BLOCK:"):
                # 如果上一个区块未处理则先处理（理论上不会出现）
                if current_block["positions"] and current_block["span"] >= MIN_SPAN:
                    tmp_blocks.append(current_block)
                    if current_block['chr'] == 'X':
                        x_blocks.append(current_block)

                # 解析SPAN值
                span_match = re.search(r'SPAN:\s*(\d+)', line)
                span = int(span_match.group(1)) if span_match else 0
                
                # 初始化新block
                current_block = {"chr": None, "positions": [], "value": block_counter % 2 + 1, "span": span}
                block_counter += 1  # 颜色交替
                
            # 检测区块结束分隔符
            elif line.startswith("********"):
                if current_block["positions"] and current_block["span"] >= MIN_SPAN:
                    tmp_blocks.append(current_block)
                    if current_block['chr'] == 'X':
                        x_blocks.append(current_block)
                current_block = {"chr": None, "positions": [], "value": None, "span": 0}
                
            # 处理数据行
            elif line and not line.startswith("#"):  # 忽略注释行
                parts = re.split(r'\s+', line)  # 兼容制表符和空格分隔
                if len(parts) >= 5:
                    chr = parts[3].replace("chr", "")  # 移除chr前缀（根据karyotype文件调整）
                    pos = int(parts[4])
                    if not current_block["chr"]:
                        current_block["chr"] = chr
                    # 确保同一区块内染色体一致
                    if current_block["chr"] == chr:
                        current_block["positions"].append(pos)
                    else:
                        print(f"Warning: Chromosome mismatch in block {block_counter}")
        
        # 处理最后一个区块
        if current_block["positions"] and current_block["span"] >= MIN_SPAN:
            tmp_blocks.append(current_block)
            if current_block['chr'] == 'X':
                x_blocks.append(current_block)

    # deal with X 
    if len(x_blocks) > 80:
        x_st = min(b['positions'][0] for b in x_blocks)
        x_ed = max(b['positions'][0] for b in x_blocks)
        x_val = x_blocks[-1]['value']

        tmp_blocks = [b for b in tmp_blocks if b['chr'] != 'X']
        tmp_blocks.append({
            'chr':'X',
            'positions': [x_st, x_ed],
            'value': 0,
            'span': x_ed - x_st
        })

    # output
    with open(output_file, 'w') as fo:
        fo.write("Chr\tStart\tEnd\tValue\n")
        for bl in tmp_blocks:
            if bl['positions']:
                st = min(bl['positions'])
                ed = max(bl['positions'])
                fo.write(f"{bl['chr']}\t{st}\t{ed}\t{bl['value']}\n")


# 调用示例
hapcut2_to_ideogram(hapblock, "contigs.tsv")
# print(f"Filtered blocks with SPAN >= {MIN_SPAN} bp")
