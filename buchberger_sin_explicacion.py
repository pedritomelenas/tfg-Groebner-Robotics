#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sympy import *


def s_polynomial(f, g):
    return expand(lcm(LM(f), LM(g)) * (1 / LT(f) * f - 1 / LT(g) * g))

def buchberger(F):
    G = list(F) #hacemos una copia de la lista F
    pairs = [] #lista de pares

    for i, f1 in enumerate(F): #i es el indice y f1 el polinomio correspondiente
        for f2 in F[i + 1:]:
            pairs.append((f1, f2)) #metemos pares en la lista pairs

    while pairs: #recorremos la lista de pares hasta que se vacie
        f1, f2 = pairs.pop(0) #almacenamos en f1 y f2 la ultima pareja y la quitamos
        s = s_polynomial(f1, f2) #calculamos su s-polinomio
        h = reduced(s, G)[1] #reducimos s con respecto a G

        if h != 0: #si su resto es distinto de 0
            for g in G:
                pairs.append((g, h)) #Se añaden los pares nuevos a comprobar

            G.append(h)#Se añade h a la lista de generadores
    return G
