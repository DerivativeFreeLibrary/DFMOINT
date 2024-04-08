import numpy as np

class filter_elem():
    ''' singolo elemento del filtro. Le info salvate sono:
        1- il punto (ovviamente)
        2- il vettore delle funzioni (di penalità)
        3- il passo usato lungo le dir dense (continue)
        4- i passi usati lungo le dir coordinate (continue)
        5- la xi per il decremento negli spostamenti lungo le discrete
        6- i passi lo spostamento lungo le direzioni primitive
        N.B. quando le primitive aumentano, deve aumentare la dimensione
        del vettore alfa_tilde per TUTTI i punti nel filtro
    '''
    x = np.empty(0,dtype=float)
    fobs = np.empty(0,dtype=float)
    alfa = 0.0
    alfa_c = np.empty(0,dtype=float)
    alfa_tilde = np.empty(0,dtype=float)
    xi = 0.0

    def __init__(self,x,f,alfa,alfa_c,alfa_tilde,xi,allones):
        self.x = np.array(x,dtype=float)
        self.fobs = np.array(f,dtype=float)
        self.alfa = np.double(alfa)
        self.xi   = np.double(xi)
        self.alfa_c = np.array(alfa_c,dtype=float)
        self.alfa_tilde = np.array(alfa_tilde,dtype=float)
        self.flag_allones = allones

    def print(self,flagx=True):
        if flagx:
            print(self.x,' ',self.fobs,' ',self.alfa,' ',self.alfa_c,' ',self.xi)
        else:
            print(self.fobs,' ',self.alfa,' ',self.alfa_c,' ',self.xi)
