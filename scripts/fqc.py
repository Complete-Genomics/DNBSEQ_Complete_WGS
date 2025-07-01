import math

def fspec(num):
    if num >= 1_000_000:
        value = num / 1_000_000
        rounded = round(value, 1)  # 保留1位小数
        return f">{rounded:.1f} million"  # 强制显示小数位（如2.0 → 2.0）
    else:
        value_k = num / 1_000
        return f">{int(value_k)}k"

f = open('report.csv')
for line in f:
    l = line.strip().split(',')
    pref, nums = l[0], l[1:]
    if pref == 'Sample':
        spec = 'summary'
    elif pref == 'cWGS bam avg depth': 
        depths = [float(num) for num in nums]
        nums = [round(float(num),1) for num in nums]
        spec = f">{int(min(nums))}X"
        nums = [str(num) for num in nums]

    elif pref.startswith('Total') or pref == 'Average fragment length' or pref.startswith('Phased contig N50'):
        nums = [float(num) for num in nums]
        spec = fspec(min(nums))
        if pref.startswith('Phased contig N50'):
            spec += 'b'
            spec.replace('millionb','Mb')

        nums = [str("{0:.1E}".format(num)) for num in nums]

    elif 'mapping rate' in pref or 'genome covered' in pref:
        _nums = [float(num.replace('%','')) for num in nums]
        spec = f">{int(min(_nums))}%"

    print(f"{pref},{spec}," + ",".join(nums))
    # elif i == 14: # phase block
    #     _nums = [int(num) for num in nums]
    #     a = max(_nums)
    #     b = int(math.ceil(a * 1000) / 1000)
    #     spec = f"<{b}"
    #     print(f"{pref},{spec}," + ",".join(nums))
f.close()


avg_depth = sum(depths) / len(depths)
nums = [str(round(d/avg_depth, 2)) for d in depths] if avg_depth else ['0' for d in depths]
a, b = round(min(depths),1), round(max(depths),1)
spec = f"between {a} and {b}"
print(f"sample variation from average,{spec}," + ','.join(nums))