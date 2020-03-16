#!/usr/bin/python

#Auteurs: Alexandre Dufour (p1054564), Jean-François Blanchette (p1099987)

#Instruction d'exécution:
#chmod a+x plast.py
#python plst.py

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
from math import pow
import argparse

MATCH = 5
MISMATCH = -4
E=0
ss=0
n=1
L = 0.192
K = 0.176
db_path=""
seed=""
sequence=""

def getdata():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-i")
    parser.add_argument("-db")
    parser.add_argument("-E", default=-4, type=int)
    parser.add_argument("-ss", default=1 * pow(10, -3), type=float)
    parser.add_argument("-seed")

    args = parser.parse_args()

    global sequence
    sequence = args.i
    global E
    E = -args.E
    global db_path
    db_path = args.db
    global ss
    ss = args.ss
    global seed
    seed = len(args.seed)
    global n
    n = len(sequence)

    print(seed)

#Extrait les séquences d'un fichier fasta
def read_fasta(path):

    fragtable = []
    with open(path, 'rt') as f:
        sequence = ''
        for line in f:
            newline = line.strip() #Retire les espaces vides et les '\n'
            if(newline.startswith('>')):
                if(sequence != ''):
                    fragtable.append(sequence)
                    sequence = ''
            else:
                sequence += newline
        fragtable.append(sequence)
    f.close()
    return fragtable


#Trouve le max de 3 valeurs
def max3(p1, p2, p3):
    if p1 < p2 and (p2 > p3) and (p2 > 0):
        return p2
    elif p2 < p3 and (p1 < p3) and (p3 > 0):
        return p3
    elif p1 > 0 and (p1 > p2) and (p1 > p3):
        return p1
    else:
        return 0

def kmer(seq, g):
    #On conserve les éléments sous la forme d'un dictionnaire
    seqk = []
    for i in range(len(seq)-g+1):
        seqk.append(seq[i:i+g])
    return seqk


#Recherche exacte par algorithme naif
def research(kmers, bank):
    maxscore = [] #Tableau qui contient les scores maximales entre les HSP et la séquence de la banque de données apparenté
    e_value = [] #Tableau qui va contenir toutes les e-value
    for i in range(len(bank)):
        fusion = [] #Tableau qui va contenir les informations de tous les HSP par séquence
        for j in range(len(kmers)):
            for k in range(len(bank[i])-len(kmers[j])):
                l = 0
                while((l < len(kmers[j])) and (k+l < len(bank[i])) and (bank[i][k+l] == kmers[j][l])):
                    l = l+1
                    if l == len(kmers[j]):
                        value = prolongement(kmers, i, j, k, bank)
                        fusion.append((kmers, j, k, bank[i], value[0], value[1]))
                        if(fusionHSPs(fusion) != None ):
                            maxscore.append(fusionHSPs(fusion))
        if(len(maxscore) > 0):
            s = max(maxscore)
            e = evalue(len(bank[i]), n, s)
            if(e < ss): #on s'assure que la valeur de la e_value soit plus significative que le second seuil
                e_value.append((e, i))
            output('unknown.fasta', s, 'ident', value[2])

'''
kmers: listes de kmer produit par la liste d'origine
seq: sequence de la banque de données
kmer: numéro du kmer d'origine dans la liste à étendre
posseq: position dans la séquence de la banque de donnée où commence le chevauchement
bank: la banque de données
'''
def prolongement(kmers, seq, kmer, posseq, bank):
    #F = forward, B = Backward
    F = 0
    B = 0
    Forward = True
    Backward = True
    score_max = 0
    value = len(kmers[0])*MATCH #le score d'origine est le score entre le kmer et la séquence
    maxvalue = 0 #score maximum atteint jusqu'à présent
    while(score_max > E and (Forward == True or Backward == True)):
        options = [0, 0, 0]  # options[0] = deux côté, options[1] = à gauche, options[2] = à droite
        #On avance vers l'avant dans la séquence de la banque
        if(Forward == True):
            F = F+1
            if(F+kmer < len(kmers) and (posseq+F+len(kmers[0]) < len(bank[seq]))):
                if(bank[seq][posseq+F+len(kmers[0])] == kmers[kmer+F][len(kmers[kmer+F])-1]): #Accès au dernier élément du prochain kmer
                    options[0] += MATCH
                    options[2] += MATCH
                else:
                    options[0] += MISMATCH
                    options[2] += MISMATCH
            else:
                Forward = False #On a atteint la fin soit des kmers, soit de la séquence de la banque de donnée
                F = F-1
        else:
            options[2] = -1000 #Comme forward est false, on veut veut pas considérer cette case dans l'évaluation du maximum
        #On avance vers l'arrière de la séquence de la bank
        if(Backward == True):
            B = B+1
            if(B <= kmer and posseq-B > 0):
                if (bank[seq][posseq-B] == kmers[kmer-B][0]):#Accès au premier élément du kmer précédent
                    options[0] += MATCH
                    options[1] += MATCH
                else:
                    options[0] += MISMATCH
                    options[1] += MISMATCH
            else:
                Backward = False
                B = B-1
        else:
            options[1] = -1000 #Comme backward est false, on ne veut pas considérer cette case du tableau
        score_max = max(options)
        if(value + score_max >= maxvalue):
            value += score_max
            maxvalue = value
            if (options.index(max(options)) == 1):
                F = F-1
            if (options.index(max(options)) == 2):
                B = B-1
    return(F, B, maxvalue)

#   0          1                           2             3      4  5
#(kmers, numéro du kmer dans kmers, position dans seq, bank[i], F, B)
#fusion_list contient toute les informations pertinente aux HSPs d'une sequence de la banque
#Cette fonction permet d'analyser et de fusionner tout les HSP d'une seule séquence de la banque de donnée
def fusionHSPs(fusion_list):
    score = []
    for i in range(len(fusion_list)):
        for j in range(len(fusion_list)):
            seq = ''
            if(i < j):
                if(fusion_list[i][1] + len(fusion_list[i][0][0]) + fusion_list[i][4] > #On vérifie si la fin d'un HSP d'un kmer chevauche le début d'un HSP d'un autre kmer
                   fusion_list[j][1]-fusion_list[i][5]):
                    dataseq = fusion_list[i][3][fusion_list[i][2]-fusion_list[i][5]:fusion_list[j][2]+fusion_list[j][4]] #Extraction de la séquence local dans la banque de donnée
                    for k in range(fusion_list[i][1] - fusion_list[i][5], fusion_list[j][1] + fusion_list[j][4]): #On parcours les kmers du début du premier à la fin du deuxième pour former la séquence HSP
                        seq += fusion_list[i][0][k][0]
                    for a in pairwise2.align.globalmx(dataseq, seq, MATCH, MISMATCH):
                        score.append(a[2]) # On stock toutes les valeurs de score obtenues, et on gardera la plus grande
    if(len(score) != 0 ):
        return max(score) #On renvoit la valeur maximale du score pour les HSP avec cette séquence de la banque de donné


def bitscore(S):
    return round((L*S-np.log(K)/np.log(2)))


def evalue(m, n, score):
    B = bitscore(score)
    e = m*n*pow(2, (-1)*B)
    return e

def output(path):
    ids = []
    sequence = read_fasta(path)

    with open(path, 'rt') as f:
        for line in f:
            newline = line.strip()  # Retire les espaces vides et les '\n'
            if (newline.startswith('>')):
                ids.append(newline)

'''
    print(ids[0]+ " score: "+ str(score) + " ident " + str(ident)+"\n"
            + "# Best HSP score:" +str(HSPscore) + ', ' + "bitscore" + str(bitscore(HSPscore)) +
            "e-value: " +str(evalue(1,1,bitscore(HSPscore))) + '\n')

    for i in ids:

        print(ids[i]+ " score: "+ str(score) + " ident " + str(ident)+"\n"
                + "# Best HSP score: " +str(HSPscore) + ', ' + "bitscore " + str(bitscore(HSPscore)) +
                "e-value: " +str(evalue(1,1,bitscore(HSPscore))))
    f.close()
'''
#Fonction principale
def main():
    getdata()
    #On retire les séquences de la banque de données
    bank = read_fasta(db_path)
    #Extraction de chaque kmers de la séquence en input. La positoon du kmer est sa position dans le tableau
    kmers = kmer(sequence, seed)
    research(kmers, bank)
    print('end')




#Execution de la fonction principale

main()
output("unknown.fasta",1,1,1)

#main()
output("unknown.fasta")