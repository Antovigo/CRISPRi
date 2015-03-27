#!/usr/bin/python
# -*- coding:utf-8 -*-

def main ():
    dist_aug = 10 # Distance minimale par rapport à l'AUG initiateur
    nbmmA = 3 # Nombre de mismatchs
    nbmmB = 3
    seqA = read_sequence ('sequences/gfp-superfolder.txt')
    seqB = read_sequence ('sequences/mCherryE.txt')
    
    for nbmm in range (0, 15, 3):
        #dual (dist_aug, seqA, nbmm, seqB, nbmm)
        generator (seqB, dist_aug, nbmm, ('ctgttttgaatggtcccaaaac', 'g'), ('aaaac', 'gttttgggaccattcaa'))

def dual (dist_aug, seqA, nbmmA, seqB, nbmmB):
    '''Deux spacers séparés par une répétition'''
    print ('GFP superfolder :')
    generator (seqA, dist_aug, nbmmA, ('aaac', 'gttttagagctatg'), ('aacagcatagctctaaaac', ''))
    print ('mCherry :')
    generator (seqB, dist_aug, nbmmB, ('ctgttttgaatggtcccaaaac', 'g'), ('aaaac', 'gttttgggaccattcaa'))

def generator (seq, dist_aug, nbmm, ohA, ohB):
    '''Génération des oligos'''
    # CHOIX DE LA CIBLE
    
    debut = pam (seq, dist_aug)
    target = seq [debut + 3 : debut + 33]
    
    # GÉNÉRER LES TÉMOINS NÉG
    # target = 'actcttaatgtcgcgaccaggatgtgcgat'
    
    # MISMATCHS
    # overhangs
    FA, FB = ohA
    RA, RB = ohB
    
    target = target [:20 - nbmm] + com (target [20 - nbmm : 20]) + target [20:]
    forward = FA + com (target [::-1]).upper() + FB # Reverse complement
    reverse = RA + target.upper() + RB
    
    print ('\033[1;31m{0} nucleotides downstream from AUG, {1} mismatchs\033[1;m'.format (debut, nbmm))
    
    print_oligos (forward, reverse)
    print_assembly (forward, reverse, FA, RB)
    print_interference (seq, target, nbmm, debut, FA, RB)
    
    print ('')
    
def read_sequence (url):
    '''Lire un fichier. Il faudra rajouter des parsers pour les différents types'''
    return open (url).read().replace(' ','').replace('\n','').lower()
    
def pam (seq, dist_aug):
    '''Trouver le PAM le plus proche de la distance voulue'''
    return seq.find ('cc', dist_aug)

def emph (seq, debut, fin):
    '''Met en emphase une partie d'une séquence en majuscules et en rouge'''
    return seq [:debut] + '\033[1;34m' + seq [debut:fin].upper() + '\033[1;m' + seq [fin:]

def base_com (base):
    '''Donne la complémentaire d'une base'''
    return {'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a'}[base]

def com (seq):
    '''Donne le complementaire de la séquence'''
    return ''.join ([base_com (i) for i in seq])
    
def print_oligos (forward, reverse):    
    print ("F : 5' " + forward + " 3'")
    print ("R : 3' " + reverse + " 5'")
    
def print_assembly (forward, reverse, FA, RB):  
    print ('')
    #print ('\033[1;31mAssemblage :\033[1;m')
    print ("5' " + " "* len (RB) + forward + " 3'")
    print ("3' " + " " * len (FA) + reverse [::-1] + " 5'")
    
def print_interference (seq, target, nbmm, debut, FA, RB):
    print ('')
    #print ('\033[1;31mInterférence :\033[1;m')
    print ('\033[1;33m ' * (28 - nbmm) + com(target [20 - nbmm:20]).replace ('t', 'u'))
    print (' ' * 8 + com(target [:20 - nbmm]).replace ('t', 'u') + ' ' * nbmm + com (target [20:]).replace ('t', 'u') + '\033[1;m')
    print (seq [debut - 5 : debut + 40])
    print (emph (com(seq [debut - 5 : debut + 3]) + ' ' * 30 + com (seq [debut + 38 : debut + 45]), 5, 7))
    print (' ' * 8 + com (seq [debut + 3 : debut + 33]))
    print ('')

main()
