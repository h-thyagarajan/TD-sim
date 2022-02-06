def fasta_format_all(gen, n):
    g1 = open("genome1.txt", "r")
    g2 = open("genome2.txt", "r")
    l1=g1.read()
    l2=g2.read()
    
    gen=49
    n=500
    l1_gen=l1[6+544007*gen:6+544007*gen+544000]
    l2_gen=l2[6+544007*gen:6+544007*gen+544000]
    
    with open('genome_'+str(gen)+'.fas', 'w') as f:
        for i in range(n):
            f.write('>No'+str(i)+'\n'+l1_gen[2+1088*i:2+1088*i+1084]+'\n')
            f.write('>No'+str(i+n)+'\n'+l2_gen[2+1088*i:2+1088*i+1084]+'\n')