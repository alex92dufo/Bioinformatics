#!/usr/bin/python

#Auteurs: Alexandre Dufour (p1054564), Jean-François Blanchette (p1099987)

#Instruction d'exécution:
#chmod a+x main.py
#python main.py


# Constante de point alloué
INDEL = -8
MISMATCH = -4
MATCH = 4

# Récupération des données de fragments dans fastq
def read_fastq(path):

    fragtable = []
    with open(path, 'rt') as f:
        for line in f:
            sequence = f.readline().rstrip()
            _ = f.readline()
            quality = f.readline().rstrip()
            fragtable.append(sequence)
    return fragtable


# Renvoie du maximum entre 3 valeurs
def max(p1, p2, p3):
    if p1 < p2 and (p2 > p3) and (p2 > 0):
        return p2
    elif p2 < p3 and (p1 < p3) and (p3 > 0):
        return p3
    elif p1 > 0 and (p1 > p2) and (p1 > p3):
        return p1
    else:
        return 0


#Produit la matrice d'alignement entre deux séquences
def buildmatrix(sequenceI, sequenceJ, path):
    scoreMatrix = [[0 for x in range(len(sequenceI)+1)] for y in range(len(sequenceJ)+1)]
    for i in range(len(scoreMatrix)):
        scoreMatrix[i][0] = i*INDEL #initialisation de la première colonne de la matrice
    #Début de l'algorithme d'alignement à partir de la seconde rangé, seconde colonne
    for i in range(len(sequenceI)):
        for j in range(len(sequenceJ)):
            aline = MATCH
            if(sequenceI[i] != sequenceJ[j]):
                aline = MISMATCH
            scoreMatrix[i+1][j+1] = max(scoreMatrix[i][j] + aline, scoreMatrix[i][j+1] + INDEL, scoreMatrix[i+1][j] + INDEL)
            #Écriture du chemin dans la matrice path
            if(scoreMatrix[i+1][j+1] == (scoreMatrix[i][j] + aline)):
                path[i+1][j+1] = 'aline'
            else:
                if (scoreMatrix[i + 1][j + 1] == (scoreMatrix[i+1][j] + INDEL)):
                    path[i+1][j+1] = 'seqIIndel'
                else:
                    path[i+1][j+1] = 'seqJIndel'
    return scoreMatrix


#Reproduit l'alignement entre deux séquences selon la matrice d'alignement et l'emplacement où le maximum se trouve dans la matrice.
def traceback(index, path, seqI, seqJ):
    seqIaline = ''
    seqJaline = ''

    m = index
    n = len(seqJ)
    longueur = 0 #Détermine la longueur du chevauchement
    for i in range(m, len(seqJ)):
        seqJaline += '-'
        seqIaline += seqI[i]
    while(m != 1):
        longueur += 1
        if(path[m][n] == 'aline'):
            seqIaline = seqI[m-1] + seqIaline
            seqJaline = seqJ[n-1] + seqJaline
            m = m-1
            n = n-1
        else:
            if (path[m][n] == 'seqJIndel'):
                seqIaline = seqI[m - 1] + seqIaline
                seqJaline = '-' + seqJaline
                m = m-1
            else:
                if (path[m][n] == 'seqIIndel'):
                    seqIaline = '-' + seqIaline
                    seqJaline = seqJ[n - 1] + seqJaline
                    n = n-1
    for i in range(len(seqJ) - n, 0, -1):
        seqJaline = seqJ[i] + seqJaline
        seqIaline = '-' + seqIaline
    print(longueur)
    print(seqJaline + '\n' + seqIaline + '\n')



#Produit la matrice de tous les score
def matricealignement(fragtable):
    #Initilialisation de la matrice des scores entre alignements
    alinementMatrix = [[0 for x in range(len(fragtable))] for y in range(len(fragtable))]
    #Comparaison des séquences entre elles
    for i in range(len(fragtable)):
        for j in range(len(fragtable)):
            highestValue = 0
            index = 0
            path = [[0 for x in range(len(fragtable[i])+1)] for y in range(len(fragtable[j])+1)] #Permet de garder en mémoire le chemin du chevauchement
            if (i != j):
                #Construction de la matrice d'alignement d'une paire de séquence
                alinement = buildmatrix(fragtable[i], fragtable[j], path)
                #Détection de la valeur la plus élevé dans la dernière colonne de la matrice d'alignement
                for m in range(len(alinement)):
                        if (alinement[m][len(alinement[m])-1] > highestValue):
                            highestValue = alinement[m][len(alinement[m])-1]
                            index = m
                alinementMatrix[i][j] = highestValue
                traceback(index, path, fragtable[i], fragtable[j])
    return alinementMatrix


#Fonction permettant de trouver l'inverse complémentaire d'une séquence
def complementaryreverse(sequences):
    reversesequence = ['' for x in range(len(sequences))]
    print(len(reversesequence[0]))
    for i in range(len(sequences)):
        for j in range(len(sequences[i])):
            if(sequences[i][j] == 'A'):
                reversesequence[i] = 'T' + reversesequence[i]
            else:
                if (sequences[i][j] == 'T'):
                    reversesequence[i] = 'A' + reversesequence[i]
                else:
                    if (sequences[i][j] == 'G'):
                        reversesequence[i] = 'C' + reversesequence[i]
                    else:
                        if (sequences[i][j] == 'C'):
                            reversesequence[i] = 'G' + reversesequence[i]
    return reversesequence

#Imprime tous les chevauchements dans le format seqX seqY score
def printall(matriceScore):
    for i in range(len(matriceScore)):
        for j in range(len(matriceScore)):
            score = matriceScore[i][j]
            if(score > 0):
                print('seq' + str(j + 1) + '\t' + 'seq' + str(i + 1) + '\t' + str(score))

#Imprime les chevauchements dont le score est supérieur à 80 dans le format seqX seqY score
def printsup(matriceScore):
    for i in range(len(matriceScore)):
        for j in range(len(matriceScore)):
            score = matriceScore[i][j]  # Détermine s'il y a un chevauchement entre les les séquences i et j
            if (score >= 80):
                print('seq' + str(j+1) + '\t' + 'seq' + str(i+1) + '\t' + str(score))

#Fonction principale
def main():
    print('Quel est le nom du fichier fastq que vous désirez analyser?')
    data = input()
    fragreads = read_fastq(data)
    print('Voici les séquences lues:')
    for i in range(len(fragreads)):
        print(fragreads[i])
    task = ''
    while (task != '0'):
        print('\n' + 'Que désirez-vous faire?' + '\n' + '1. Inverser la séquence et trouver son complémentaire' + '\n' + '2. Trouver les scores de chevauchement des séquences' + '\n' + '0. Quitter')
        task = input()
        if(task == '1'): #Produit le complémentaire inverse des séquences passées en paramètre.
            complementary = complementaryreverse(fragreads)
            print(complementary)
            print('\n')
            print('Quel est le nom du fichier fastq possédant les séquences forward?')
            data = input()
            fragreads = read_fastq(data) #On reprend la séquence

            #On échange les séquences des reads reverse trouvé dans le fichier reads.fq par leur inverse complémentaire
            fragreads[5] = complementary[0]
            fragreads[6] = complementary[1]
            fragreads[8] = complementary[2]
            fragreads[9] = complementary[3]
            fragreads[11] = complementary[4]
            fragreads[14] = complementary[5]
            fragreads[17] = complementary[6]

            print("Les séquences ont été enregistrées en mémoire.")

        if(task == '2'): #On mesure les scores entre les séquences paires à paires.
            print('Veuillez patienter, nous traitons les données.')
            matricescore = matricealignement(fragreads)

            print("Que désirez-vous voir?" + '\n' + "1. Les paires de séquence ayant un chevauchement" + '\n' + "2. Les paires de séquence ayant un chevauchement dont le score est minimalement de 80." + '\n' + "3. La matrice de score")
            task = input()

            if(task == '1'):
                printall(matricescore)

            if(task == '2'):
                printsup(matricescore)

            if(task == '3'):
                for i in range(len(matricescore)):
                    print(matricescore[i])


main()

