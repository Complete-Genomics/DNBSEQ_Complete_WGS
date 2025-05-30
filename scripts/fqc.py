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
i = 0
for line in f:
    i += 1
    l = line.strip().split(',')
    pref, nums = l[0], l[1:]
    if i == 1:
        print('Sample,summary,' + ','.join(nums))
    elif i == 2: # avg depth
        depths = [float(num) for num in nums]
        nums = [round(float(num),1) for num in nums]
        spec = f">{int(min(nums))}X"
        nums = [str(num) for num in nums]
        print(f"stLFR bam avg depth,{spec}," + ','.join(nums))

    elif i in [3,4,5,6,8,9,12,13,15,16]: #snp, indel,lfr cnt, lfr avg len, n50, pc per sample
        nums = [float(num) for num in nums]
        spec = fspec(min(nums))
        if i in [9,15]:
            spec += 'b'
            spec.replace('millionb','Mb')
        elif i == 16:
            spec = fspec(max(nums))
            spec.replace('>','<')

        nums = [str("{0:.1E}".format(num)) for num in nums]
        print(f"{pref},{spec}," + ",".join(nums))

    elif i in [7,10,11]: #bc split rate, map rate, cov
        _nums = [float(num.replace('%','')) for num in nums]
        spec = f">{int(min(_nums))}%"
        print(f"{pref},{spec}," + ",".join(nums))

    elif i == 14: # phase block
        _nums = [int(num) for num in nums]
        a = max(_nums)
        b = int(math.ceil(a * 1000) / 1000)
        spec = f"<{b}"
        print(f"{pref},{spec}," + ",".join(nums))
f.close()


avg_depth = sum(depths) / len(depths)
nums = [str(round(d/avg_depth, 2)) for d in depths]
a, b = round(min(depths),1), round(max(depths),1)
spec = f"between {a} and {b}"
print(f"sample variation from average,{spec}," + ','.join(nums))