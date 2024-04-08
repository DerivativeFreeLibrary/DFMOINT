import numpy as np
import pandas as pd
import ctypes

class cache_factory():
    '''
    n : numero di variabili
    m : numero di vincoli
    q : numero di funzioni obiettivo
    '''

    n        = 0
    m        = 0
    q        = 0
    soglia   = 0
    ncache_hits = 0
    cur_dim  = 0
    size     = 0
    use_pandas = True
    store = 0
    store_pd = 0

    def __init__(self,n,m,q):
        self.n = n
        self.m = m
        self.q = q
        if self.use_pandas:
            self.soglia   = np.finfo(np.float32).eps/np.sqrt(n)
        else:
            self.soglia   = 1.e-6

        df    = pd.DataFrame(columns=['id','x','fob','constr'])
        self.store_pd = pd.DataFrame(df,index=df.index.copy())
        #self.store_pd = pd.DataFrame(columns=['id','x','fob','constr'])
        self.store = []
        #self.store = [dict() for k in range(30000)]
        self.ncache_hits = 0
        self.cur_dim  = 0
        self.size = 0

    def find(self, x):
        '''
        cerca nella cache il punto x
        ritorna true e fob + constr se lo trova
        altrimenti torna false e [] []
        '''

        self.ncache_hits += 1
        #if self.ncache_hits > 350:
        #    print('xxxxxx : ',x)


        # if not trovato and not i.empty:
        #     print('beccato (1)!!!!!')
        #     input()
        # if trovato and i.empty:
        #     print('beccato (2)!!!!!')
        #     print('NO PANDAS: trovato')
        #     print('   PANDAS: non trovato')
        #     print(x,hash(np.float16(x).tostring()),hash(np.float32(x).tostring()))
        #     print(x1,hash(np.float16(x1).tostring()),hash(np.float32(x1).tostring()))


        if self.use_pandas:

            #print('id == "'+str(hash(np.float32(x).tostring()))+'"')
            #print(self.store_pd[['id','x','fob','constr']].query('id == "'+str(hash(np.float32(x).tostring()))+'"'))
            #input()
            #i = self.store_pd[['id','x','fob','constr']].query('id == "'+str(hash(np.float32(x).tostring()))+'"')
            i = self.store_pd[self.store_pd.id == hash(np.float32(x).tostring())]
            if i.empty:
                return False,np.nan,np.empty(0),np.empty(0)
            else:
                return True,i.id.values[0],i.fob.values[0],i.constr.values[0]
            #ind = ctypes.c_size_t(hash(x.tostring())).value
            #if ind in self.store_pd.index:
            #    return True, ind, self.store_pd.loc[ind].fob, self.store_pd.loc[ind].constr
            #else:
            #    return False,np.nan,np.empty(0),np.empty(0)

        else:
            id = np.nan
            x1 = np.array(0)
            fob  = np.array(0)
            constr = np.array(0)
            trovato = False
            for ind,c in enumerate(self.store):
                if c:
                    #print('.',end='')
                    if  np.all(np.abs(x-c['x']) <= self.soglia):
                        id = ind
                        fob = np.copy(c['fob'])
                        constr = np.copy(c['constr'])
                        x1 = np.copy(c['x'])
                        trovato = True
                        #print('CACHE: x=',x)
                        #print('CACHE: z=',c['x'])
                        #print('xxxxxx : ',ind,c['x'])
            return trovato,id,fob,constr

        if False:
            for ind,c in enumerate(self.store):
                if c:
                    #print('.',end='')
                    if  np.all(np.abs(x-c['x']) <= self.soglia):
                        #print('CACHE: x=',x)
                        #print('CACHE: z=',c['x'])
                        #print('xxxxxx : ',ind,c['x'])
                        return True,ind,c['fob'],c['constr']
            #print()
            return False,np.nan,np.empty(0),np.empty(0)

    def insert(self,x,fob,constr):
        '''
        inserisce in store la tripla x,fob,constr
        se non c'è già
        '''
        flag,ind,fobz,constrz = self.find(x)
        if not(flag):
            if self.use_pandas:
                #self.store_pd = self.store_pd.append({'id':hash(np.float32(x).tostring()),
                #    'x':np.copy(x), 'fob':np.copy(fob),'constr':np.copy(constr)},ignore_index=True)
                self.store_pd = pd.concat([self.store_pd, pd.DataFrame([{'id':hash(np.float32(x).tostring()),
                    'x':np.copy(x), 'fob':np.copy(fob),'constr':np.copy(constr)}])])
                #ind = ctypes.c_size_t(hash(x.tostring())).value
                #self.store_pd.loc[ind] = {'id':ind,
                #    'x':np.copy(x), 'fob':np.copy(fob),'constr':np.copy(constr)}
            else:
                self.store.append({'x': np.copy(x),'fob': np.copy(fob),'constr': np.copy(constr)})
            #if self.size < 30000:
            #    self.size += 1
            #if self.cur_dim >= 30000:
            #    self.cur_dim = 0
            #self.store[self.cur_dim] = {'x': np.copy(x),'fob': np.copy(fob),'constr': np.copy(constr)}
            #self.cur_dim += 1
        return

    def print(self):
        for ind,c in enumerate(self.store):
            print('CACHE: x=',c['x'])
