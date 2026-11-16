import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend("agg")

# data in depth.txt follow "depth,count" format.
# depth start from 0 
count_table = np.loadtxt("03.depth.txt",delimiter=' ',dtype=int)
total_pos=np.sum(count_table[:,1])
count_table_101 = np.zeros(101)
for item in count_table:
    if item[0] < 100 :
        count_table_101[item[0]] += item[1]
    else:
        count_table_101[100] += item[1]

freq_table_101 = count_table_101 * 100 / total_pos

plt.figure()
plt.bar(np.arange(0,101),freq_table_101,width=0.3,color=['b']*101)
plt.xlim(0,101)
plt.ylabel('Fraction of basese(%)')
plt.xlabel('Sequencing depth')
plt.xticks([0,20,40,60,80,100],[0,20,40,60,80,'100+'])
plt.savefig('03.depth_freq.png')

cumulative_table_101=np.zeros(101)
for i in range(101):
    cumulative_table_101[i] = np.sum(freq_table_101[i:])
plt.figure()
plt.plot(np.arange(0,101),cumulative_table_101,color='r')
plt.xlim(0,101)
plt.ylabel('Fraction of basese(%)')
plt.xlabel('Cumulative Sequencing depth')
plt.xticks([0,20,40,60,80,100],[0,20,40,60,80,'100+'])
plt.savefig('03.cumulative_depth_freq.png')

average = np.sum(count_table[:,0]*count_table[:,1])/total_pos
print("Average sequencing depth %.6f (only unique mapped reads)"%average)
print("Coverage %.6f %%"%(100.0-freq_table_101[0]))
print("Coverage at least 4X \t%.6f"%cumulative_table_101[3]);
print("Coverage at least 10X\t%.6f"%cumulative_table_101[9]);
print("Coverage at least 20X\t%.6f"%cumulative_table_101[19]);
print("Coverage at least 30X\t%.6f"%cumulative_table_101[29]);

print("depth\tbases\tpercent\tcumulative_percent\n")
for i in range(101):
    print("%d\t%d\t%.6f\t%.6f"%(i,count_table_101[i],
          freq_table_101[i],cumulative_table_101[i]))
