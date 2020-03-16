#!/usr/bin/python

#Auteurs: Alexandre Dufour (p1054564), Jean-François Blanchette (p1099987)

#Instruction d'exécution:
#chmod a+x AlignementMultiple.py
#python AlignementMultiple.py

from Bio.SubsMat.MatrixInfo import blosum62
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


A = -10
B = -1
INFINITY = -100000000
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

#Construit le tableau d'alignement de deux séquences
def buildmatrix(sequenceI, sequenceJ):
    scorematrix = [[0 for x in range((len(sequenceJ)+1))] for y in range(3*(len(sequenceI)+1))]
    #Initialisation de la première ligne et de la première colonne
    for i in range(len(sequenceJ)):
        scorematrix[0][i] = INFINITY
        scorematrix[1][i] = INFINITY
        scorematrix[2][i] = A + (i * B)
    for i in range(3*len(sequenceI)):
        if(i % 3 == 0):
            scorematrix[i][0] = INFINITY
        if (i % 3 == 1):
            scorematrix[i][0] = A + (i * B)
        if (i % 3 == 2):
            scorematrix[i][0] = INFINITY
    #Début de l'algorithme d'alignement à partir de la seconde rangé, seconde colonne
    for i in range(3, 3*len(sequenceI)):
        for j in range(1, len(sequenceJ)):
            #Recherche du score dans blosum62
            sequI = sequenceI[i//3]
            sequJ = sequenceJ[j-1]
            try:
                scorematch = blosum62[(sequI, sequJ)]
            except KeyError:
                scorematch = blosum62[(sequJ, sequI)]

            #Matrice M
            if(i % 3 == 0):
                scorematrix[i][j] = max3(scorematrix[i - 3][j - 1] + scorematch, scorematrix[i - 2][j - 1] + scorematch, scorematrix[i - 1][j - 1] + scorematch)
            #Matrice X
            if(i % 3 == 1):
                scorematrix[i][j] = max2(scorematrix[i-3][j] + B, scorematrix[i-4][j] + A + B)
            #Matrice Y
            if(i % 3 == 2):
                scorematrix[i][j] = max2(scorematrix[i][j-1] + B, scorematrix[i-2][j-1] + A + B)
    return scorematrix

def max2(p1, p2):
    if (p1 < p2):
        return p2
    return p1


#Trouve la séquence Sc ayant le meilleur score d'alignement avec les autres séquences
def sc(alinementmatrix):
    scmax = 0
    for i in range(len(alinementmatrix)):
        scoresc = 0
        for j in range(len(alinementmatrix[i])):
            scoresc += alinementmatrix[i][j]
        if(scoresc > scmax):
            scmax = scoresc
            sequenceI = i+1
    print('scoresc = ' + str(scmax))
    print('sequence = ' + str(sequenceI))

#Permet de toruver la séquence consensus
def matconsensus(seq1, seq2, seq3, seq4, seq5):
    consensus = [[0 for x in range(len(seq1))] for y in range(21)]
    alignement = ""
    for i in range((len(seq5))):
        char1 = seq1[i]
        char2 = ''
        char3 = ''
        char4 = ''
        char5 = ''
        x1 = 1
        x2 = 0
        x3 = 0
        x4 = 0
        x5 = 0
        if(seq2[i] == char1):
            x1 += 1
        else:
            char2 = seq2[i]
            x2 = 1
        if(seq3[i] == char1):
            x1 += 1
        else:
            if (seq3[i] == char2):
                x2 += 1
            else:
                char3 = seq3[i]
                x3 = 1
        if (seq4[i] == char1):
            x1 += 1
        else:
            if (seq4[i] == char2):
                x2 +=1
            else:
                if( seq4[i] == char3):
                    x3 += 1
                else:
                    char4 = seq4[i]
                    x4 = 1
        if (seq5[i] == char1):
            x1 += 1
        else:
            if (seq5[i] == char2):
                x2 +=1
            else:
                if( seq5[i] == char3):
                    x3 += 1
                else:
                    if (seq5[i] == char4):
                        x4 += 1
                    else:
                        char5 = seq5[i]
                        x5 = 1
        letter = max(x1, x2, x3, x4, x5)
        if(letter == x1):
            alignement += char1
        else:
            if (letter == x2):
                alignement += char2
            else:
                if(letter == x3):
                    alignement += char3
                else:
                    if(letter == x4):
                        alignement += char4
                    else:
                        alignement += char5
    print("\nAlignement Consensus:\n" + alignement + "\n")

    print("Pourcentage d'identité:")
    similaire = 0
    for i in range(len(alignement)):
        if (seq1[i] == alignement[i]):
            similaire += 1
    pourcentage = similaire * 100 / len(alignement)
    print('seq2 = ' + str(pourcentage))

    similaire = 0
    for i in range(len(alignement)):
        if (seq2[i] == alignement[i]):
            similaire += 1
    pourcentage = similaire * 100 / len(alignement)
    print('seq4 = ' + str(pourcentage))

    similaire = 0
    for i in range(len(alignement)):
        if (seq3[i] == alignement[i]):
            similaire += 1
    pourcentage = similaire * 100 / len(alignement)
    print('seq3 = ' + str(pourcentage))

    similaire = 0
    for i in range(len(alignement)):
        if (seq4[i] == alignement[i]):
            similaire += 1
    pourcentage = similaire * 100 / len(alignement)
    print('seq5 = ' + str(pourcentage))

    similaire = 0
    for i in range(len(alignement)):
        if (seq5[i] == alignement[i]):
            similaire += 1
    pourcentage = similaire * 100 / len(alignement)
    print('seq1 = ' + str(pourcentage))

#inspiré de https://stackoverflow.com/questions/5686211/is-there-a-function-that-can-calculate-a-score-for-aligned-sequences-given-the-a?fbclid=IwAR1BGxd50Ek_ftVbkeX-8sAyMQcTLfR16on1GdEX8vCo8BkRapeZ2CZoe78
#Permet de retrouver le score de deux séquence déjà aligné
def score_pairwise(seq1, seq2, gap_s, gap_e):
    score = 0
    gap = False
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if(seq1[i] == '-' and seq2[i] == '-'): #Si deux gaps match ensemble, on ne modifie pas le score
            score += 0
            gap = False
        else:
            if not gap: #Si on n'est pas déjà dans un gap
                if '-' in pair:
                    gap = True
                    score += gap_s
                else:
                    try:
                        score = blosum62[(seq1[i], seq2[i])]
                    except KeyError:
                        score = blosum62[(seq2[i], seq1[i])]
            else: #Si on est déjà dans un gap
                if '-' not in pair: #Si on n'est pas dans un gap
                    gap = False
                    try:
                        score = blosum62[(seq1[i], seq2[i])]
                    except KeyError:
                        score = blosum62[(seq2[i], seq1[i])]
                else:
                    score += gap_e
    return score


#Fonction principale
def main():
    fragtable = read_fasta('sequences.fasta')
    alinementmatrix = [[0 for x in range(len(fragtable))] for y in range(len(fragtable))]
    for i in range(len(fragtable)):
        for j in range(len(fragtable)):
            if(i != j):
                matrix = buildmatrix(fragtable[i], fragtable[j])
                #Va rechercher la valeur maximal dans le tableau de score
                alinementmatrix[i][j] = max3(3*matrix[len(fragtable[i])-1][len(fragtable[j])-1], matrix[3*len(fragtable[i])-2][len(fragtable[j])-1], matrix[3*len(fragtable[i])-3][len(fragtable[j])-1])
    for k in range(len(alinementmatrix)):
        print(alinementmatrix[k])
    sc(alinementmatrix)

    blosum = matlist.blosum62

    #Permet de trouver l'alignement entre deux séquences
    alignement = pairwise2.align.globalds(fragtable[1],fragtable[0],blosum, A, B)
    print(alignement)

    #Permet de retrouver l'alignement consensus
    matconsensus("MEKVPGEMEIERRERSEELSEAERKAVQATWARLYANCEDVGVAILVRFFVNFPSAKQYFS----QFKHMEEPLEMER-SPQLRKHACRVMGALNTVVENL—HDPEKV---SSVLSLVGKAHALKHKVEPVYFKILSGV---ILEVIAEEFANDFPPETQRAWAKLRGLIYSHVTAAYKEVGWVQQVPNATTPPATLPSSGP",
                 "M----------------GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDK----FKHLKSEDEMKA-SEDLKKHGATVLTALGGILKKKGHHEAE--------IKPLAQSHATKHKIPVKYLEFISEC---IIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGF----------------QG-",
                 "M----------------VLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFP----HF-------DLSHGSAQIKGHGKKVVAALIEAANH--IDD---IAG—TLSKLSDLHAHKLRVDPVNF---KLLGQCFLVVVAIHHPAALTPEVHASLDKFLCAVGTVLTAKYR-----------------------",
                 "MERLESEL------------------IRQSWRAVSRSPLEHGTVLFSRLFALEPSLLPLFQYNGRQFSSPEDCLS---SPEFLDHIRKVMLVIDAA--VTNVEDLSSLEEY-LATLGR----KHRAVGVRLSSFST---VGESLLYMLEKCLGPDFTPATRTAWSQLYGAVVQAMSR-----GW----------------DGE",
                 "MG------EIGFTEKQEAL-------VKESWEILKQDIPKYSLHFFSQILEIAPAAKGLFS----FLRDSDE—VPHNN-PKLKAHAVKVFKMTCETAIQLR—EEGKVVVA---DTTLQYLGSIHLKSGVIDP-HFEVVKEALLRTLKEGLGEKYNEEVEGAWSQ----AYDHLALAIK---------------TEMKQEES")

    #Permet de retrouver le score d'un alignement
    print("\nScore d'alignement = " + str(score_pairwise("MERLESEL------------------IRQSWRAVSRSPLEHGTVLFSRLFALEPSLLPLFQYNGRQFSSPEDCLS---SPEFLDHIRKVMLVIDAA--VTNVEDLSSLEEY-LATLGR----KHRAVGVRLSSFST---VGESLLYMLEKCLGPDFTPATRTAWSQLYGAVVQAMSR-----GW----------------DGE",
                   "MG------EIGFTEKQEAL-------VKESWEILKQDIPKYSLHFFSQILEIAPAAKGLFS----FLRDSDE--VPHNN-PKLKAHAVKVFKMTCETAIQLR--EEGKVVVA---DTTLQYLGSIHLKSGVIDP-HFEVVKEALLRTLKEGLGEKYNEEVEGAWSQ----AYDHLALAIK---------------TEMKQEES",
                   A, B)))

main()
	