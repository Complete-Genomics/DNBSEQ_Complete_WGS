import pysam, sys
from multiprocessing import Pool

bam, cpu, output, li = sys.argv[1:]

id_map = {}
f = open(li)
for line in f:
    line = line.strip()
    id = line.split('#')[-1]
    id_map[line] = id
f.close()

def process_chunk(chunk):
    # 处理一个数据块（自定义逻辑）
    for read in chunk:
        read.query_name = id_map.get(read.query_name, read.query_name)
    return chunk

# 分块读取 BAM 文件
with pysam.AlignmentFile(bam, "rb") as bam_in:
    chunks = [list(bam_in.fetch(until_eof=True))]  # 示例简单分块，实际需优化分块策略

# 使用多进程池
with Pool(processes=cpu) as pool:  # 指定进程数
    processed_chunks = pool.map(process_chunk, chunks)

# 合并结果并写入文件
with pysam.AlignmentFile(output, "wb", header=bam_in.header) as bam_out:
    for chunk in processed_chunks:
        for read in chunk:
            bam_out.write(read)