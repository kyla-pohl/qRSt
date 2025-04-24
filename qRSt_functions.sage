Sym = SymmetricFunctions(QQ)
s = Sym.schur()





def a(llambda,i,j):
   try:
       return Partition(llambda).arm_lengths()[i][j]
   except IndexError:
       raise IndexError("Aieee! {}, {}, {}".format(llambda, i, j))






def l(llambda,i,j):
    return Partition(llambda).leg_lengths()[i][j]








def alpha(nu, llambda, Q, T):
    
    assert (llambda <= nu), "uh-oh, llambda is not inside nu:"+str(llambda)+str(nu)
    
    listllambda = list(llambda)
    listllambda.append(0)
    for i in range(len(nu)):
        if listllambda[i] != nu[i]:
            r = i
    
    R = []
    for(i,j) in Partition(llambda).cells():
        if i == r:
            R.append((i,j))
            
    listllambdaconj = list(llambda.conjugate())
    listllambdaconj.append(0)
    for i in range(len(nu.conjugate())):
        if listllambdaconj[i] != nu.conjugate()[i]:
            c = i
    
    C = []
    for(i,j) in Partition(llambda.conjugate()).cells():
        if i == c:
            C.append((i,j))
    
    
    result = prod( (1-(Q^a(llambda,i,j) * T^(1+l(llambda,i,j))))/(1-(Q^a(nu,i,j) * T^(1+l(nu,i,j)))) for (i,j) in R )
    result *= prod( (1-Q^(a(llambda,j,i)+1) * T^l(llambda,j,i)) / (1-Q^(1+a(nu,j,i)) * T^l(nu,j,i)) for (i,j) in C)
    
    
    return result






def n(nu, llambda):
    return sum([sum(nu.leg_lengths()[i]) for i in range(len(nu))]) - sum([sum(llambda.leg_lengths()[i]) for i in range(len(llambda))])






def beta(llambda, mu, Q, T):
    return 1 / (alpha(llambda, mu, Q, T))






def nprime(llambda, mu):
    return sum([sum(llambda.arm_lengths()[i]) for i in range(len(llambda))]) - sum([sum(mu.arm_lengths()[i]) for i in range(len(mu))])







def gamma(nu, llambda, mu, Q, T):
    return ((1-Q^ (nprime(llambda, mu) - nprime(nu, llambda)) * T^ (n(nu, llambda) - n(llambda, mu)) ) * (1-Q^ (nprime(llambda, mu) - nprime(nu, llambda) + 1) * T^ (n(nu, llambda) - n(llambda, mu) - 1) )) / ((1-Q) * (1-T))





def F_QT(llambda, mu, Q, T):
    llambda = Partition(llambda)
    mu = Partition(mu)
    
    size = llambda.size()
        
    
    # We start with the case mu in D(llambda).
    L = Partitions(size - 1, outer = llambda)
    for i in range(len(L)):
        if mu == L[i]:

            # This is the list of all partitions nu that are one bigger than llambda and contain llambda.
            B = Partitions(size + 1, inner = llambda)


            probB = [0]


            # For each of these partitions, we need to define the probability of ending up there.
            for b in B:
                probb = T^(n(b,llambda) - n(llambda, mu) - 1) * alpha(b, llambda, Q, T) * beta(llambda, mu, Q, T) / gamma(b, llambda, mu, Q, T)

                # This list will be our "buckets".
                probB.append(probB[len(probB)-1] + probb)

            #assert probB[len(probB)-1] == 1, "uh oh -- the probabilities for mu do not add up to one!"

            rand = random()

            # Which "bucket" does rand fall into? That determines the partition F should return.
            for j in range(len(probB)-1):
                if rand > probB[j]:
                    if rand < probB[j+1]:
                        return B[j]
                    
    
    # Second, we have the mu = llambda case.
    else:
        assert mu == llambda, "uh oh -- mu and llambda are different!"
        
        # This is the list of all partitions nu that are one bigger than llambda and contain llambda.
        B = Partitions(size + 1, inner = llambda)
        
        
        probB = [0]

        
        # For each of these partitions, we need to define the probability of ending up there.
        for b in B:
            probb = T^(n(b,llambda)) * alpha(b, llambda, Q, T)
            
            # This list will be our "buckets".
            probB.append(probB[len(probB)-1] + probb)
        
        #assert probB[len(probB)-1] == 1, "uh oh -- the probabilities for llambda do not add up to one!"
        
        rand = random()
        
        # Which "bucket" does rand fall into? That determines the partition F should return.
        for j in range(len(probB)-1):
            if rand > probB[j]:
                if rand < probB[j+1]:
                    return B[j]










def qRSt(sigma,Q,T):
    
    M = matrix(Permutation(sigma))
    
    for i in range(100):
        if g in SymmetricGroup(i):
            if g not in SymmetricGroup(i-1):
                m = i
    
    N = [ [1 for j in range(m+1)] for i in range(m+1)]
    for i in range(m+1):
        N[0][i] = Partitions(0)[0]
    for i in range(1,m+1):
        N[i][0] = Partitions(0)[0]
    for i in range(1,m+1):
        for j in range(1,m+1):
            N[i][j] = None


    for u in range(m+1):
        for v in range(m+1):
            if N[u][v] == None:

                k = sum(M[i,j] for i in range(u) for j in range (v))

                p = list(N[u-1][v])
                q = list(N[u][v-1])
                if len(q) - len(p) > 0:
                    for i in range(len(q) - len(p)):
                        p.append(0)

                if len(p) - len(q) > 0:
                    for i in range(len(p) - len(q)):
                        q.append(0)    

                Plist = list(Partitions(k, inner = [max(p[i], q[i]) for i in range(len(p))]))


                if len(Plist)==1:
                    N[u][v] = Plist[0]
                else:
                    assert p==q, "uh-oh - p isn't equal to q:"+str(p)+str(q)
                    r = N[u-1][v-1]
                    N[u][v] = F_QT(p,r,Q,T)
    
    MatrixTableauP = matrix({SkewPartition([ N[i][m], N[i-1][m]]).cells()[0]: (i) for i in range(1,m+1)})

    ListTableauP = [list(MatrixTableauP[i]) for i in range(MatrixTableauP.nrows())]

    ZerolessListTableauP = [[i for i in ListTableauP[j] if i != 0] for j in range(len(ListTableauP))]

    TableauP = StandardTableau(ZerolessListTableauP)

    
    
    MatrixTableauQ = matrix({SkewPartition([ N[m][i], N[m][i-1]]).cells()[0]: (i) for i in range(1,m+1)})

    ListTableauQ = [list(MatrixTableauQ[i]) for i in range(MatrixTableauQ.nrows())]

    ZerolessListTableauQ = [[i for i in ListTableauQ[j] if i != 0] for j in range(len(ListTableauQ))]

    TableauQ = StandardTableau(ZerolessListTableauQ)

    
    #print('P =',TableauP)
    #print('Q =',TableauQ)
    
    return TableauP, TableauQ





def qRSt_output(g,q,t):
    Output = qRSt(g, q, t)
    Output[0].pp()
    print('-----------')
    Output[1].pp()

