import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
#yoooooo
def igdna(length):
    dna1,dna2,dna3 = "", "", ""
    for count in range(length):
        dna1+=random.choice("CGTA")
        dna2+=random.choice("CGTA")
        dna3+=random.choice("CGTA")
    return dna1,dna2,dna3

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
    gene=str(gene)
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

def mutate(dna, mu):
    for j in range(mu):
        dna_list = list(dna)
        mutation_site = np.random.randint(0, len(dna_list) - 1)
        dna_list[mutation_site] = np.random.choice(list('ATCG'))
        dna = ''.join(dna_list)
    return dna

def recombination(genome1,genome2,nrecomb):
    rec_bp=np.sort(np.random.choice(range(len(genome1)),nrecomb))
    g1=genome1
    g2=genome2
    for i in range(len(rec_bp)):
        genome1_new=genome1[:rec_bp[i]]+genome2[rec_bp[i]:]
        genome2_new=genome2[:rec_bp[i]]+genome1[rec_bp[i]:]
        g1=genome1_new
        g2=genome2_new
    return g1, g2

def tot_fitness(gene_l1, gene_l2, w_gene):
    gene1a=gene_l1[:4]
    gene1b=gene_l1[184:364]
    gene1c=gene_l1[544:724]
    gene1d=gene_l1[904:1084]
    gene2a=gene_l2[:4]
    gene2b=gene_l2[184:364]
    gene2c=gene_l2[544:724]
    gene2d=gene_l2[904:1084]
    w1b=fitness(gene1b,w_gene[1])
    w1c=fitness(gene1c,w_gene[2])
    w1d=fitness(gene1d,w_gene[3])
    w2b=fitness(gene2b,w_gene[1])
    w2c=fitness(gene2c,w_gene[2])
    w2d=fitness(gene2d,w_gene[3])
    lw=[np.average([w1b,w2b]),np.average([w1c,w2c]),np.average([w1d,w2d])]
    tot_w=np.average(lw)
    if gene1a != 'ATGC' and gene1a != 'ATAT':
        tot_w=0
    if gene2a != 'ATGC' and gene2a != 'ATAT':
        tot_w=0
    return lw, tot_w

def intiate_pop(n,mu,w_gene):
    l1_gene1, l1_gene2, l1_gene3, l1_gene4 =[], [], [], []
    l2_gene1, l2_gene2, l2_gene3, l2_gene4 =[], [], [], []
    for i in range(n):
        l1_gene1+=['ATAT']
        l1_gene2+=[mutate(w_gene[1],mu)]
        l1_gene3+=[mutate(w_gene[2],mu)]
        l1_gene4+=[mutate(w_gene[3],mu)]
        l2_gene2+=[mutate(w_gene[1],mu)]
        l2_gene3+=[mutate(w_gene[2],mu)]
        l2_gene4+=[mutate(w_gene[3],mu)]
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
        igdna1=igdna(180)
        l1+= [l1_gene1[i]+igdna1[0]+l1_gene2[i]+igdna1[1]+l1_gene3[i]+igdna1[2]+l1_gene4[i]]
        igdna2=igdna(180)
        l2+= [l2_gene1[i]+igdna2[0]+l2_gene2[i]+igdna2[1]+l2_gene3[i]+igdna2[2]+l2_gene4[i]]
        w=tot_fitness(l1[i],l2[i],w_gene)
        lw+= [w[0]]
        ltot_w+= [w[1]]
#    print(lw)
    return l1, l2, lw, ltot_w

def surv_adult(l):
    lw=[]
    for i in range(len(l[2])):
        lw+=[l[2][i][0]]
    lw=np.array(lw)
    lindex=np.arange(len(lw))
    ladults=np.delete(lindex, np.where(lw<0.5))
    l_final_adults=[]
    for i in range(len(ladults)):
        if l[3][ladults[i]]>0:
            l_final_adults+=[ladults[i]]
    return np.array(l_final_adults)

def n_matings(l,s):
    lw=[]
    for i in range(len(l[2])):
        lw+=[l[2][i][1]]
    lw=np.array(lw)
    #add male specific part
    nmates=np.zeros(len(lw))
    for i in range(len(lw)):
        if i in s:
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
                fec=int((8*lw[female][0]/n_mates)*(np.log10(lw[female][2]+1))*3*lw[males[j]][0])
                l_fec+=[fec]
            tot_fec+=[l_fec]
    return tot_fec

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

def timeseries(n,imu,gen,new_rf,mut,w_gene):
    l=intiate_pop(n,imu,w_gene)
    lw2_out_mean,lw3_out_mean,lw4_out_mean=[],[],[]
    lw2_out_std,lw3_out_std,lw4_out_std=[],[],[]
    l_out=[]
    for t in range(gen):
        if t<=25:
            rf=0.8
        else:
            rf=new_rf
        female_index=[]
        male_index=[]
        for i in range(len(l[1])):
            if l[1][i][:4]=="ATAT":
                female_index+=[i]
            else:
                male_index+=[i]
        s_index=surv_adult(l)
        n_mate=n_matings(l,s_index)
        l_pair=mating_pairs(n_mate, female_index,male_index)
        eggs=fecundity(l_pair,l[2])
        f_eggs=final_eggs(eggs,n)

        lrecomb=[]
        for i in range(n):
            if np.random.uniform()<rf:
                nrecomb=1
                n_rand=np.random.uniform()
                if n_rand<rf*0.8:
                    nrecomb=2
                if n_rand<rf*0.4:
                    nrecomb=3
                if n_rand<rf*0.2:
                    nrecomb=4
            else:
                nrecomb=0
            lrecomb+=[nrecomb]
        l_offs_parents=[]
        for i in range(len(f_eggs)):
            for j in range(len(f_eggs[i])):
                for k in range(f_eggs[i][j]):
                    l_offs_parents+=[[l_pair[i][0],l_pair[i][1][j]]]

        l1,l2,lw,ltot_w=[],[],[],[]
        for i in range(len(l_offs_parents)):
            dam_i=l_offs_parents[i][0]
            sire_i=l_offs_parents[i][1]
            nrecomb=lrecomb[i]
            og1=np.random.choice([0,1])
            og2=np.random.choice([0,1])
            recg=np.random.choice([0,1])
            g1,g2='',''
            if nrecomb==0:
                g1=str(l[og1][dam_i])
                g2=str(l[og2][sire_i])
            else:
                if np.random.uniform()>0.5:
                    g1=recombination(str(l[og1][dam_i]),str(l[og2][sire_i]),nrecomb)[recg]
                    g2=str(l[abs(og2-1)][sire_i])
                else:
                    g1=str(l[abs(og1-1)][dam_i])
                    g2=recombination(str(l[og1][dam_i]),str(l[og2][sire_i]),nrecomb)[recg]
            l1+=[mutate(g1,int(mut*len(g1)))]
            l2+=[mutate(g2,int(mut*len(g2)))]

        for i in range(n):
            w=tot_fitness(l1[i],l2[i],w_gene)
            ltot_w+=[w[1]]
            lw+= [w[0]]
        l=[l1,l2,lw,ltot_w]
        l_out+=[l]
        lw2,lw3,lw4=[],[],[]
        for i in range(len(lw)):
            lw2+=[lw[i][0]]
            lw3+=[lw[i][1]]
            lw4+=[lw[i][2]]

        lw2_out_mean+=[np.average([lw2])]
        lw3_out_mean+=[np.average([lw3])]
        lw4_out_mean+=[np.average([lw4])]
        lw2_out_std+=[np.std([lw2])]
        lw3_out_std+=[np.std([lw3])]
        lw4_out_std+=[np.std([lw4])]
        print(t,new_rf, rf)
    lw_out_mean=[lw2_out_mean,lw3_out_mean,lw4_out_mean]
    lw_out_std=[lw2_out_std,lw3_out_std,lw4_out_std]
    return l_out, lw_out_mean, lw_out_std

n=500
imu=30
gen=50
mut=0.001
wgene1="ATGC"
wgene2="AATTTCGCCGACGTGATGACATTCCAGGCAGTGCCTCTGCCGCCGGACCCCTCTCGTGATTGGGTAGCTGGACATGCCCTTGTAAGATATAACAAGAGCCTGCCTGTCTAATGATCTCACGGCGAAAGTCGGGGAGACAGCAGCGGCTGCAGACATTATACCGCAACAACACTAAGGTGA"
wgene3="GATAACTCCGTAATTGACTACGCGTTCCTCTAGACCTTACTTGACCGGATACAGTGTCTTTGACACGTTTATGGGTTACAGCAATCACATCCAAGACTGGCTATGCACGAAGCAACTCTTGAGTGTTAAAATGTTGACCCCTGTATTTGGGATGCGGGTAGTAGATGAGTGCAGGGACTC"
wgene4="CGAGGTCAAGTACATTACCCTCTCATAGGGGGCGTTCTAGATCACGTTACCACCATATCATTCGAGCATGACACCATCTCCGCTGTGCCCATCCTAGTAGTCATTATTCCTATCACGCTTTCGAGTGTCTGGTGGCGGATATCCCCCACGAATGAAAATGTTTTTCGCTGACAGTCATAT"
w_gene=[wgene1,wgene2,wgene3,wgene4]

for r in np.arange(0,1.1,0.1):
    rf=np.round(r,1)
    l_out=timeseries(n,imu,gen,rf,mut,w_gene)

    def plot_w_score(l_out,gen):
        w_score_df=pd.DataFrame(dict(g2_mean=l_out[1][0],g3_mean=l_out[1][1],g4_mean=l_out[1][2],g2_std=l_out[2][0],g3_std=l_out[2][1],g4_std=l_out[2][2],gen=range(gen)))
        sns.set_style("darkgrid")
        ax1 = sns.lineplot(x="gen", y="g2_mean", data=w_score_df)
        ax1.fill_between(w_score_df["gen"], y1=w_score_df["g2_mean"] - (w_score_df["g2_std"]/np.sqrt(n)), y2=w_score_df["g2_mean"] + (w_score_df["g2_std"]/np.sqrt(n)), alpha=.3)
        ax2 = sns.lineplot(x="gen", y="g3_mean", data=w_score_df)
        ax2.fill_between(w_score_df["gen"], y1=w_score_df["g3_mean"] - (w_score_df["g3_std"]/np.sqrt(n)), y2=w_score_df["g3_mean"] + (w_score_df["g3_std"]/np.sqrt(n)), alpha=.3)
        ax3 = sns.lineplot(x="gen", y="g4_mean", data=w_score_df)
        ax3.fill_between(w_score_df["gen"], y1=w_score_df["g4_mean"] - (w_score_df["g4_std"]/np.sqrt(n)), y2=w_score_df["g4_mean"] + (w_score_df["g4_std"]/np.sqrt(n)), alpha=.3)
        plt.xlabel("Generations")
        plt.ylabel("w_scoreu "+u"\u00B1"+" SE")
        plt.legend(labels=["gene2","gene3","gene4"])
        plt.tight_layout()
        plt.savefig('H:\sim\w_score1_rf_'+str(rf)+'.png')
        plt.close()

    def genomes_txt_raw(l_out):
        genome1_data, genome2_data=[],[]
        for i in range(gen):
            genome1_data+=[l_out[0][i][0]]
            genome2_data+=[l_out[0][i][1]]

        with open('genome1_rf_'+str(rf)+'.txt', 'w') as f:
            for i in range(len(genome1_data)):
                if i<10:
                    f.write('gen'+str(0)+str(i)+'\n'+str(genome1_data[i])+'\n')
                else:
                    f.write('gen'+str(i)+'\n'+str(genome1_data[i])+'\n')

        with open('genome2_rf_'+str(rf)+'.txt', 'w') as g:
            for i in range(len(genome2_data)):
                if i<10:
                    g.write('gen'+str(0)+str(i)+'\n'+str(genome2_data[i])+'\n')
                else:
                    g.write('gen'+str(i)+'\n'+str(genome2_data[i])+'\n')

    def fasta_format_all(gen, n, l_out):
        l1_gen=l_out[0][gen][0]
        l2_gen=l_out[0][gen][1]

        with open('genome_'+str(gen)+'_all_rf_'+str(rf)+'.fas', 'w') as f:
            for i in range(n):
                f.write('>No'+str(i)+'a\n'+l1_gen[i]+'\n')
                f.write('>No'+str(i)+'b\n'+l2_gen[i]+'\n')

    def fasta_format_gene2(gen, n, l_out):
        l1_gen=l_out[0][gen][0]
        l2_gen=l_out[0][gen][1]

        with open('genome_'+str(gen)+'_gene2_rf_'+str(rf)+'.fas', 'w') as f:
            for i in range(n):
                f.write('>No'+str(i)+'a\n'+l1_gen[i][184:364]+'\n')

                f.write('>No'+str(i)+'b\n'+l2_gen[i][184:364]+'\n')

    def fasta_format_gene3(gen, n, l_out):
        l1_gen=l_out[0][gen][0]
        l2_gen=l_out[0][gen][1]

        with open('genome_'+str(gen)+'_gene3_rf_'+str(rf)+'.fas', 'w') as f:
            for i in range(n):
                f.write('>No'+str(i)+'a\n'+l1_gen[i][544:724]+'\n')
                f.write('>No'+str(i)+'b\n'+l2_gen[i][544:724]+'\n')

    def fasta_format_gene4(gen, n, l_out):
        l1_gen=l_out[0][gen][0]
        l2_gen=l_out[0][gen][1]

        with open('genome_'+str(gen)+'_gene4_rf_'+str(rf)+'.fas', 'w') as f:
            for i in range(n):
                f.write('>No'+str(i)+'a\n'+l1_gen[i][904:1084]+'\n')
                f.write('>No'+str(i)+'b\n'+l2_gen[i][904:1084]+'\n')

    def fasta_format_intergene(gen, n, l_out):
        l1_gen=l_out[0][gen][0]
        l2_gen=l_out[0][gen][1]

        with open('genome_'+str(gen)+'_intergene_rf_'+str(rf)+'.fas', 'w') as f:
            for i in range(n):
                f.write('>No'+str(i)+'a\n'+l1_gen[i][4:184]+l1_gen[i][364:544]+l1_gen[i][724:904]+'\n')
                f.write('>No'+str(i)+'b\n'+l2_gen[i][4:184]+l2_gen[i][364:544]+l2_gen[i][724:904]+'\n')

    plot_w_score(l_out,gen)
    genomes_txt_raw(l_out)
    fasta_format_all(49,n,l_out)
    fasta_format_gene2(49,n,l_out)
    fasta_format_gene3(49,n,l_out)
    fasta_format_gene4(49,n,l_out)
    fasta_format_intergene(49,n,l_out)
