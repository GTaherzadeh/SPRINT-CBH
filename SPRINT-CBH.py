###################################
#!/usr/bin/env python
#coding=utf-8
# ./SPRINT-CBH.py pdbid
# Ghazaleh Taherzadeh
###################################
"""predict Carbohydrate-binding residues for a new chain"""
import sys
import string
import math
import os
import numpy
import scipy.io

# function
def encode_restype(res):
    AAs = 'ARNDCQEGHILKMFPSTWYV'
    code = []
    for a in AAs:
        if res == a:
            code.append('1')
        else:
            code.append('0')    
    return code

def compute_entropy(dis_list):
    if sum(dis_list) == 0:# zero-list: (0,0,...,0)
        return 0.0
    prob_list = map(lambda x:(x+0.0)/sum(dis_list),dis_list)
    ent = 0.0
    for prob in prob_list:
        if prob != 0:
            ent -= prob*math.log(prob,len(dis_list))
    return ent

def compute_pcc(x,y):
    mean_x = (sum(x)+0.0)/len(x)
    dev_x = sum([(i-mean_x)**2 for i in x])
    mean_y = (sum(y)+0.0)/len(y)
    dev_y = sum([(i-mean_y)**2 for i in y])
    if dev_x == 0 or dev_y == 0:
        return 0.0
    ret = 0.0
    for i in xrange(len(x)):
        ret += (x[i]-mean_x)*(y[i]-mean_y)
    return ret/math.sqrt(dev_x*dev_y)
#/*************************************************/
pid = sys.argv[1] 
win = 4 
base_path = 'path_to_your_sequence_dir'
rsa_path = base_path+'/Spider2/'
pssm_path = base_path+'/PSSM/'
fasta_path = base_path+'/fasta/'
fea_path = base_path+'/fea/'
res_path = base_path+'Results/'
error_file = base_path+'error.log'
#/*************************************************/
fin = file(rsa_path+pid+'.spd3','r')
sp = fin.readlines()[1:]
fin.close()
rsa_res = ''.join([x.split()[1] for x in sp])
rsa_pre = [string.atof(x.split()[6]) for x in sp]
fin = file(pssm_path+pid+'.pssm','r')
pssm = fin.readlines()
fin.close()
if len(pssm[-6].split()) != 0 or pssm[3].split()[0] != '1': 
    print 'error on reading pssm, line -6 is not a spare line;\
     or line 3 is not the first line'
    sys.exit(1)
pssm = pssm[3:-6]
fin = file(fasta_path+pid+'.seq','r')
ann = fin.readlines()
fin.close()
if len(ann) != 2:
    print 'check sequence',pid
    sys.exit(1)
fastaseq = ann[1].split()[0]
if not fastaseq == rsa_res:
    print 'Sequence inconsistent!'
    print 'fasta: ',fastaseq
    print '  rsa: ',rsa_res
    exit(1)
fout = file(fea_path+pid+'.info','w')
fout.write('>%s\n' %pid)
pos = 0
for i in xrange(len(fastaseq)):
    res = fastaseq[i]
    fout.write('%5d%5s%5s'%(i+1,res,res))
    if pssm[pos].split()[1] == res:# check for residue type
        for p_e in pssm[pos].split()[2:22]:
            fout.write(':%2s' %p_e)
        for p_e in pssm[pos].split()[22:42]:
            fout.write(':%3s' %p_e)
        fout.write(':%5s' %pssm[pos].split()[42])
    else:
        print 'Error reading pssm file!'
        flog = file(error_file,'a')
        flog.write(pid+': error on writing pssm, %s:%s\n' \
        %(pssm[pos].split()[1],res))
        flog.close()
        sys.exit(1)
    if rsa_res[pos] == res:
        fout.write(':%5.3f' %rsa_pre[pos])
    else:
        print 'Error reading rsa file!'
        flog = file(error_file,'a')
        flog.write(pid+': error on writing rsa, %s:%s\n' %(rsa_res[pos],res))
        flog.close()
        sys.exit(1)
    pos += 1
    fout.write('\n')
fout.close()
#/*************************************************/
# build feature file
fin = file(fea_path+'%s.info'%pid,'r')
info = fin.readlines()[1:]
fin.close()
output = file(fea_path+'%s.fea'%pid,'w')
seq_len = len(info)
# initial
out_list = []
b = []
for i in xrange(len(info)):
    out_list.append([])
for i in xrange(len(info)):
    b.append([])
pssm = map(lambda x:map(lambda y:'%7.5f' %(1/(1+math.pow(math.e,-string.atoi(y)))),x.split(':')[1:21]),info)
pssm_t = []
for i in xrange(20):
    pssm_t.append('%7.5f' %(1/(1+math.e**0)))
for i in xrange(win):
    pssm.insert(0,pssm_t)
    pssm.append(pssm_t)
for i in xrange(win,seq_len+win):
    for j in xrange(i-win,i+win+1):
        out_list[i-win].append(':'.join(pssm[j]))
wop = ['%7.5f' %compute_entropy(z) for z in map(lambda x:map(lambda y:string.atoi(y),x.split(':')[21:41]),info)]
for i in xrange(win):
    wop.insert(0,'%7.5f' %(0))
    wop.append('%7.5f' %(0))
for i in xrange(win,seq_len+win):
    for j in xrange(i-win,i+win+1):
        out_list[i-win].append(wop[j])
pssm = [[string.atoi(y) for y in x.split(':')[1:21]] for x in info]
pssm_t = []
for i in xrange(20):
    pssm_t.append(0)
for i in xrange(win):
    pssm.insert(0,pssm_t)
    pssm.append(pssm_t)
for i in xrange(win,seq_len+win):
    for j in xrange(i-win,i+win+1):
        if j != i:
            out_list[i-win].append('%.4f'%(compute_pcc(pssm[i],pssm[j])))
rt = map(lambda x:x.split(':')[0][-1],info)
for i in xrange(win):
    rt.insert(0,'X')
    rt.append('X') 
for i in xrange(win,seq_len+win):
    for j in xrange(i-win,i+win+1):
        out_list[i-win].append(':'.join(encode_restype(rt[j])))
for i in xrange(0,len(fastaseq)):
    out_list[i].append(str(len(fastaseq)/float(1000)))
rsa = map(lambda x:x.split(':')[42].split()[0],info)
for i in xrange(win):
    rsa.insert(0,'1.000')
    rsa.append('1.000')
for i in xrange(win,seq_len+win):
    for j in xrange(i-win,i+win+1):
        out_list[i-win].append(rsa[j])
rsa = [string.atof(x.split(':')[42]) for x in info]
for i in xrange(win):
    rsa.insert(0,1.0)
    rsa.append(1.0)
for i in xrange(win,seq_len+win):
    for j in xrange(1,win+1):
        out_list[i-win].append('%.4f'%(sum(rsa[i-j:i+j+1])/(2*j+1)))
temp_ter = []
ter = map(lambda x:x.split(':')[0][-1],info)
for i in xrange(0,len(fastaseq)):
    if i == 0 or i == 1 or i == 2 or i == len(fastaseq)-3 or i == len(fastaseq)-2 or i == len(fastaseq)-1:
       temp_ter.append(1) 
    else: 
       temp_ter.append(0)  
for i in xrange(win):
    temp_ter.insert(0,'0')
    temp_ter.append('0') 
for i in xrange(win,seq_len+win):
    for j in xrange(i-win,i+win+1):
        out_list[i-win].append(str(temp_ter[j]))
# output to fea file
myList = []
for i in xrange(len(out_list)):
    output.write('-1')
    output.write('\t')
    out_list1 = []
    out_list1.append(':'.join(out_list[i]))
    myList = [t.split(':') for t in out_list1]
    for k in myList:
        for p in enumerate(k):
            D = (p)
	    output.write(':'.join([str(D[0]+1),D[1]]))
	    output.write('\t')
    output.write('\n')
output.close()
#/*************************************************/
# Classifier
outfile = file(res_path+'%s.prob'%pid,'w')
sys.path.append ('path_to_libsvm/python')
import svm
from svm import *
from svmutil import *
#y, x = svm_read_problem('./All_Data.txt')
#m = svm_train(y, x, '-s 0 -t 2 -g 0.05 -c 1 -b 1')
#svm_save_model('Fullmodel.model', m)
m = svm_load_model(base_path+'Fullmodel.model')
a, b = svm_read_problem(fea_path+'%s.fea'%pid)
p_label, p_acc, p_val = svm_predict(a, b, m,'-b 1')
for item in p_val:
  outfile.write("%s\n" % item[0])
outfile.close()
