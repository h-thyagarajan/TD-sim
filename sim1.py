import random 
import numpy as np
import matplotlib.pyplot as plt

gene1="ATGC"
gene2="TAAAACGTCTCGATCGCTTGCGCAACTTGT"
gene3="GAAGTGTCTACCATCCCTAAGCCCATTTCC"
gene4="CGCATATTAACCCCTGATTGTATCCGCATC"

w_gene=[gene1,gene2,gene3,gene4]
def score_match(subject, query):
    score = 0
    for i in range(0,len(subject)):
        subject_base = subject[i]
        query_base = query[i]

        if subject_base == query_base:
            score = score + 1
        else:
            score = score
    return score/len(subject)
def non3_fun(gene):
    nt_listgene=[]
    for i in range(len(gene)):
        if i%3!=0 or i==0:
            nt_listgene+=[gene[i]]
    nt_gene=''.join(nt_listgene)
    return nt_gene

def fitness(gene, w_gene):
    if gene.find('AAAAA')>=1 or gene.find('GGGGG')>=1 or gene.find('TTTTT')>=1 or gene.find('CCCCC')>=1:
        w=0
    nt_gene = non3_fun(gene)
    nt_w_gene = non3_fun(w_gene)
    score=score_match(nt_gene,nt_w_gene)
    if score>=0.75:
        w=score
    if 0.75>score_match(nt_gene,nt_w_gene)>=0.5:
        w=0.5*score
    if 0.5>score_match(nt_gene,nt_w_gene):
        w=0
    return w

def tot_fitness(gene_l1, gene_l2, w_gene):
    gene1a=gene_l1[:4]
    gene1b=gene_l2[:4]
    gene2a=gene_l1[34:64]
    gene3a=gene_l1[94:124]
    gene4a=gene_l1[154:184]
    gene2b=gene_l2[34:64]
    gene3b=gene_l2[94:124]
    gene4b=gene_l2[154:184]
    w2a=fitness(gene2a,w_gene[1])
    w3a=fitness(gene3a,w_gene[2])
    w4a=fitness(gene4a,w_gene[3])
    w2b=fitness(gene2b,w_gene[1])
    w3b=fitness(gene3b,w_gene[2])
    w4b=fitness(gene4b,w_gene[3])
    lw=[np.average([w2a,w2b]),np.average([w3a,w3b]),np.average([w4a,w4b])]
    tot_w=np.average(lw) 
    if gene1a != 'ATAT':
        tot_w=0
    if gene1b !='ATGC' and gene1b !='ATAT':
        tot_w=0
    return lw, tot_w

def mutate(dna, mu):
    for j in range(mu):
        dna_list = list(dna)
        mutation_site = random.randint(0, len(dna_list) - 1)
        dna_list[mutation_site] = np.random.choice(list('ATCG'))
        dna = ''.join(dna_list)
    return dna

def intiate_pop(n,mu,w_gene):

    l1_gene1, l1_gene2, l1_gene3, l1_gene4 =[], [], [], []
    l2_gene1, l2_gene2, l2_gene3, l2_gene4 =[], [], [], []
    for i in range(n):
        l1_gene1+=['ATAT']    
        l1_gene2+=[mutate(gene2,mu)]
        l1_gene3+=[mutate(gene3,mu)]
        l1_gene4+=[mutate(gene4,mu)]
        l2_gene2+=[mutate(gene2,mu)]
        l2_gene3+=[mutate(gene3,mu)]
        l2_gene4+=[mutate(gene4,mu)]
        if i<n//2:
            l2_gene1+=['ATGC']
        else:
            l2_gene1+=['ATAT']
  
    def igdna(length):
        dna1,dna2,dna3 = "", "", ""
        for count in range(length):
            dna1+=random.choice("CGTA")
            dna2+=random.choice("CGTA")
            dna3+=random.choice("CGTA")
        return dna1,dna2,dna3
    
    l1,l2,ltot_w,lw=[],[],[],[]
    
    for i in range(n):
        igdna1=igdna(30)
        l1+= [l1_gene1[i]+igdna1[0]+l1_gene2[i]+igdna1[1]+l1_gene3[i]+igdna1[2]+l1_gene4[i]]
        igdna2=igdna(30)
        l2+= [l2_gene1[i]+igdna2[0]+l2_gene2[i]+igdna2[1]+l2_gene3[i]+igdna2[2]+l2_gene4[i]]
        w=tot_fitness(l1[i],l2[i],w_gene)
        ltot_w+= [w[1]]
        lw+= [w[0]]
    return l1, l2, ltot_w, lw

n=500
mu=10
l=intiate_pop(n,mu,w_gene)
female_index=[]
male_index=[]
for i in range(len(l[1])):
    if l[1][i][:4]=="ATAT":
        female_index+=[i]
    else:
        male_index+=[i]
        
def surv_adult(l):
    lw=[]
    for i in range(len(l[3])):
        lw+=[l[3][i][0]]
    lw=np.array(lw)
    lindex=np.arange(len(lw))
    ladults=np.delete(lindex, np.where(lw<0.5))
    return ladults  
s_index=surv_adult(l)

def n_matings(l,s):
    lw=[]
    for i in range(len(l[3])):
        lw+=[l[3][i][1]]
    lw=np.array(lw)
    #add male specific part
    nmates=np.zeros(len(lw))
    for i in range(len(lw)):
        if i in s_index:
            if lw[i]>0.85:
                nmates[i]=3
            elif lw[i]>0.7:
                nmates[i]=2
            elif lw[i]>0.55:
                nmates[i]=1
            else:
                nmates[i]=0
        else:
            nmates[i]=0
    return nmates
n_mate=n_matings(l,s_index)
        
def mating_pairs(n_matings, female_index,male_index):
    l_pair=[]
    for i in female_index:
        f_mates=n_matings[i]
        weight=[]
        for j in male_index:
            weight+=[n_matings[j]]
        p_weight=np.array(weight)/sum(weight)
        male=np.random.choice(np.array(male_index),int(f_mates), p=p_weight,replace='False')
        l_pair+=[[i,male]]
    return l_pair
l_pair=mating_pairs(n_mate, female_index,male_index)

def fecundity(l_pair, lw):
    tot_fec=[]
    for i in range(len(l_pair)):
        female=l_pair[i][0]
        males=l_pair[i][1]
        n_mates=len(males)
        l_fec=[]
        if n_mates==0:
            tot_fec+=[[0]]
        else:
            for j in range(n_mates):
                fec=int((10*lw[female][0]-n_mates)*(np.log10(lw[female][2]+1))*3*lw[males[j]][0])
                l_fec+=[fec]
            tot_fec+=[l_fec]
    return tot_fec
eggs=fecundity(l_pair,l[3])

def final_eggs(eggs,n): 
    f_concate=np.concatenate(eggs)
    #print(f_concate, sum(f_concate))
    while sum(f_concate)>n:
        i=np.random.choice(len(f_concate))
        if f_concate[i]>0:
            f_concate[i]=f_concate[i]-1
    #print(f_concate, sum(f_concate))
    final_eggs=[]
    
    counter=0
    for i in range(len(eggs)):
        lf=f_concate[counter:counter+len(eggs[i])]
        counter+=len(eggs[i])
        final_eggs+=[lf]
    return np.array(final_eggs)
final_eggs=final_eggs(eggs,n)

#def sires_i(final_eggs,l_pair):
#    mated_males=[]
#    for i in range(len(l_pair)):
#        if len(l_pair[i][1])>0:
#            mated_males+=[l_pair[i][1]]
#    sires=[]
#    for i in range(len(mated_males)):
#        for j in range(len(final_eggs[i])):
#            if final_eggs[i][j]>0:
#                sires+=[mated_males[i][j]]
#    sires_index=np.unique(np.array(sires))
#    return sires_index
#sires_i=sires_i(final_eggs,l_pair)

def recombination(genome1,genome2,nrecomb):
    n=len(genome1)
    rec_bp=np.sort(np.random.choice(range(n),nrecomb))
    for i in range(len(rec_bp)):
        genome1_new=genome1[:rec_bp[i]]+genome2[rec_bp[i]:]
        genome2_new=genome2[:rec_bp[i]]+genome1[rec_bp[i]:]
        genome1=genome1_new
        genome2=genome2_new
    return genome1, genome2

rf=0.8
lrecomb=[]
for i in range(n):
    if np.random.uniform()<0.8:
        nrecomb=1
        n_rand=np.random.uniform()
        if n_rand<0.6:
            nrecomb=2
        if n_rand<0.3:
            nrecomb=3
        if n_rand<0.1:
            nrecomb=4
    else:
        nrecomb=0
    lrecomb+=[nrecomb]

l_offs_parents=[]
for i in range(len(final_eggs)):
    for j in range(len(final_eggs[i])):
        for k in range(final_eggs[i][j]):
            l_offs_parents+=[[l_pair[i][0],l_pair[i][1][j]]]
new_l=[]
l1,l2,lw,ltot_w=[],[],[],[]
for i in range(len(l_offs_parents)):
    dam_i=l_offs_parents[i][0]
    sire_i=l_offs_parents[i][1]
    nrecomb=lrecomb[i]
    if nrecomb>0:
        m1=np.random.choice([0,1])
        m2=np.random.choice([0,1])
        m3=np.random.choice([dam_i,sire_i])
        l1+=[str(recombination(l[m1][dam_i],l[m2][sire_i],nrecomb))]
        if m3==dam_i:
            l2+=[str(l[abs(m1-1)][dam_i])]
        else:
            l2+=[str(l[abs(m2-1)][sire_i])]
    else:
        m1=np.random.choice([0,1])
        m2=np.random.choice([0,1])
        m3=np.random.choice([dam_i,sire_i])
        if m3==dam_i:
            l2+=[str(l[abs(m1-1)][dam_i])]
            l1+=[str(l[abs(m2)][sire_i])]
        else:
            l2+=[str(l[abs(m2-1)][sire_i])]
            l1+=[str(l[abs(m1)][dam_i])]
            
for i in range(n):
    w=tot_fitness(l1[i],l2[i],w_gene)
    ltot_w+=[w[1]]
    lw+= [w[0]]
new_l=[l1,l2,ltot_w,lw]