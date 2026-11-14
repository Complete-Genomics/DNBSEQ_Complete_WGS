import argparse
import numpy as np
#import matplotlib as mpl
#mpl.use('agg')
#import matplotlib.pyplot as plt

###########################################################
# parse args & set global variables
###########################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='input pos file', required=True, type=str)
#parser.add_argument('-o', help='output prefix', required=True, type=str)
parser.add_argument('-b', help='break gap threshold, bigger gap break barcode into several LFR fragments', required=True, type=int)
parser.add_argument('--max-chr-id', help='skip reads with ref_id greater than this value; default 24 = chr1-22,X,Y,M in GRCh38', type=int, default=24, dest='max_chr_id')
#parser.add_argument('-p', help='read-pair threshold of high quality LFR', required=True, type=int)
#parser.add_argument('-l', help='length threshold of high quality LFR', required=True, type=int)
args = parser.parse_args()

file_name = args.i
#prefix = args.o
gap_threshold=args.b
#pair_threshold=args.p
#length_threshold=args.l

###########################################################
# utils class
###########################################################

# frequence & distribution class
class Freq:
    def __init__(self):
        self.datas={}
    
    def Touch( self , key , value=1 ):
        if key in self.datas:
            self.datas[key] += value
        else :
            self.datas[key] =value

    def Print( self , promt ):
        for ( k , v ) in self.datas.items():
            print('%s\t%d\t%d'%(promt,k,v))
    
    def Draw( self ,filename, title , xlabel , ylabel ,bins=None,xlog=False,ylog=False,xlim=None,ylim=None):
        plt.switch_backend('agg')
        if bins == None :
            plt.figure(figsize=(8,6))
        else :
            plt.figure(figsize=(15,7))
        if bins != None:
            plt.subplot(121)
            plt.scatter(self.datas.keys(),self.datas.values())
            if xlog :  plt.xscale('log')
            if ylog :  plt.yscale('log')
            if xlim :  plt.xlim(xlim)
            if ylim :  plt.ylim(ylim)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.subplot(122)
            bin_freq=Freq()
            for b in bins:
                bin_freq.Touch(0)
            for( k , v ) in self.datas.items() :
                start = 0 
                for step in bins :
                    if step <= k :
                        start = step
                    else : break
                bin_freq.Touch(step,v)
            X = bin_freq.datas.keys()
            Y = bin_freq.datas.values()
            plt.bar(np.arange(len(X))+0.5,Y,width=0.9)
            for x,y in zip(np.arange(len(X)),Y):
                plt.text(x+0.5,y+0.05,'%d' %y, ha='center',va='bottom')
            plt.xticks(np.arange(len(X)),bins)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
        else :
            plt.scatter(self.datas.keys(),self.datas.values())
            if xlog :  plt.xscale('log')
            if ylog :  plt.yscale('log')
            if xlim :  plt.xlim(xlim)
            if ylim :  plt.ylim(ylim)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
        plt.suptitle(title)
        #fig.subplots_adjust(top=0.85)
        #plt.tight_layout()
        plt.savefig(filename)

# information of 1 LFR.
class LFRInfo:
    def __init__( self ):
        self.chr_id = -1
        self.start = -1
        self.end = -1
        self.read_num = 0
        self.readpair_num = 0 
        self.read_ids = {}
    
    # same LFR means same chr_id && distance little then gap_threshold
    def isInSameLFR( self , chr_id , pos):
        if self.chr_id == -1 :
            return True
        else :
            if chr_id != self.chr_id :
                return False ;
            if pos < self.end :
                raise Exception("input position is not sorted ! chr_id = %d , pos = %d " % (chr_id , pos ))
            if pos - self.end + 1 >= gap_threshold:
                return False
            else:
                return True
    
    # make sure you call isInSameLFR and return value is True before call addPos.
    # otherwise incorrect data will be generated.
    def addPos( self , read_id , chr_id , pos):
        if self.chr_id == -1 :
            self.chr_id = chr_id
            self.start = pos
        self.end = pos
        self.read_num +=1
        if not read_id in self.read_ids :
            self.read_ids[read_id]=pos
            self.readpair_num += 1
        
    def Length( self ):
        return self.end - self.start + 1
    
# information of 1 barcode
class BarcodeInfo:
    def __init__(self , barcode_id):
        self.id = id
        self.lfrs = []
        self.curr_lfr = LFRInfo()
        self.read_num = 0
    
    def addMappingInfo( self , read_id , chr_id , pos ):
        if self.curr_lfr.isInSameLFR( chr_id,pos ):
            self.curr_lfr.addPos( read_id , chr_id , pos )
        else:
            self.saveLFR()
            self.addMappingInfo( read_id , chr_id , pos )
        self.read_num += 1
    
    def saveLFR( self ):
        self.lfrs.append(self.curr_lfr)
        self.curr_lfr=LFRInfo()
        
    def LFRNum( self ):
        return len(self.lfrs)
    def ReadPairNum(self):
        rp=0
        for lfr in self.lfrs :
            rp += lfr.readpair_num
        return rp

###########################################################
# Load & Parseing input data
###########################################################

# all_data contain 4 columns : read_id  barcode_id  ref_id  map_pos
# data in all_data should be :
#       first grouped by ref_id
#       then  sorted by map_pos
all_data = np.loadtxt(file_name, delimiter='\t',dtype=int)
if args.max_chr_id >= 0:
    all_data = all_data[all_data[:, 2] <= args.max_chr_id]
barcodes={}
for i in all_data:
    if not i[1] in barcodes :
        barcodes[i[1]] = BarcodeInfo(i[1])
    barcodes[i[1]].addMappingInfo(i[0],i[2],i[3])
for ( _ , bs ) in barcodes.items():
    bs.saveLFR()
all_data =  None

############################################################
## Collection metrics
############################################################
#barcode_num = len(barcodes) # total barcode num 
#barcode_read_num_freq = Freq()
#barcode_lfr_freq = Freq()
#lfr_num=0
#lfr_readpair_num_freq = Freq()
#lfr_read_num_freq = Freq()
#lfr_length_freq = Freq()
#high_quality_lfr_num = 0
#high_quality_lfr_fac = 0.0
#
#for ( barcode_id , barcode_info ) in barcodes.items():
#    barcode_read_num_freq.Touch(barcode_info.read_num)
#    lfr_num += barcode_info.LFRNum()
#    barcode_lfr_freq.Touch(barcode_info.LFRNum())
#    for lfr in barcode_info.lfrs :
#        lfr_readpair_num_freq.Touch(lfr.readpair_num)
#        lfr_read_num_freq.Touch(lfr.read_num)
#        lfr_length_freq.Touch(lfr.Length())
#        if lfr.readpair_num >= pair_threshold  and lfr.Length() >= length_threshold :
#            high_quality_lfr_num+=1
#
#high_quality_lfr_fac = (high_quality_lfr_num * 1.0) / (lfr_num * 1.0 )
###########################################################
# Print barcode details
###########################################################
f=open('05.barcode.details.txt', mode='w', encoding='utf-8')
f.write('barcode-id\tread-pair-num\tread-num\tlfr-num\n')
for ( barcode_id , barcode_info ) in barcodes.items():
    f.write("%d\t%d\t%d\t%d\n" % (barcode_id,
             barcode_info.ReadPairNum(),
             barcode_info.read_num,
             barcode_info.LFRNum()))
f.close()
f=open('05.long_fragment.details.txt', mode='w', encoding='utf-8')
f.write('barcode-id\tref-name\tlfr-readpair-num\tlfr-read-num\tlfr-length\n')
for ( barcode_id , barcode_info ) in barcodes.items():
    for lfr in barcode_info.lfrs:
        f.write("%d\t%s\t%d\t%d\t%d\n" % (barcode_id,
                                  lfr.chr_id,
                                  lfr.readpair_num,
                                  lfr.read_num,
                                  lfr.Length()))
f.close()
############################################################
## Generate reports
############################################################
#print('############ REPORT START #################################')
#print('')
#print('')
#print('------------General report ---------------------')
#print('The used gap_threshold : %d ' % gap_threshold)
#print('Total barcode number : %d ' % barcode_num )
#print('Total LFR number : %d ' % lfr_num)
#print('The used read-pair threshold of high-quality-LFR : %d ' % pair_threshold)
#print('The used length threshold of high-quality-LFR : %d ' % length_threshold)
#print('Total high-quality-LFR number %d ' % high_quality_lfr_num)
#print('Total high-quality-LFR proportion %.3f ' % high_quality_lfr_fac)
#print('------------General report end -----------------')
#print('')
#print('')
#print('------------Barcode read number freq -----------')
#print('------------USE `grep ^BRSN` to filter ---------')
#barcode_read_num_freq.Print('BRSN')
##barcode_read_num_freq.Draw(filename="%s_barcode_read_num_freq.png" % prefix ,
##                           title ='mapped read number freq of barcode',
##                           xlabel='mapped read number',
##                           ylabel='count',
##                           bins=[0,10,20,30,50,100,150,200,300,400,500,1000],
##                           ylog=True)
#print('------------Barcode read number end ------------')
#print('')
#print('')
#print('------------Barcode lfr number freq ------------')
#print('------------USE `grep ^BLFR` to filter ---------')
#barcode_lfr_freq.Print('BLFR')
##barcode_lfr_freq.Draw(filename="%s_barcode_lfr_freq.png" % prefix ,
##                           title ='LFR number freq of barcode',
##                           xlabel='LFR number',
##                           ylabel='count')
#print('------------Barcode lfr number end -------------')
#print('')
#print('')
#print('------------LFR readpair number freq -----------')
#print('------------USE `grep ^LRPN` to filter ---------')
#lfr_readpair_num_freq.Print('LRPN')
##lfr_readpair_num_freq.Draw(filename="%s_lfr_readpair_num_freq.png" % prefix ,
##                           title ='mapped readpair number of LFR',
##                           xlabel='readpair number',
##                           ylabel='count',
##                           bins=[0,10,20,30,50,100,150,200],
##                           ylog=True)
#print('------------LFR readpair number end ------------')
#print('')
#print('')
#print('------------LFR read number freq ---------------')
#print('------------USE `grep ^LRSN` to filter ---------')
#lfr_read_num_freq.Print('LRSN')
##lfr_read_num_freq.Draw(filename="%s_lfr_read_num_freq.png" % prefix ,
##                           title ='mapped read number of LFR',
##                           xlabel='read number',
##                           ylabel='count',
##                           bins=[0,10,20,30,50,100,150,200,300,400,500,1000],
##                           ylog=True)
#print('------------LFR read number end ----------------')
#print('')
#print('')
#print('------------LFR length freq --------------------')
#print('------------USE `grep ^LFRL` to filter ---------')
#lfr_length_freq.Print('LFRL')
##lfr_length_freq.Draw(filename="%s_lfr_length_freq.png" % prefix ,
##                           title ='length of LFR',
##                           xlabel='LFR length',
##                           ylabel='count',
##                           bins=[0,100,500,1000,5000,10000,20000,30000,50000,100000,200000],
##                           xlog=True,
##                           ylog=True)
#print('------------LFR length freq end ----------------')
#print('')
#print('')
#print('############ REPORT END ###################################')
