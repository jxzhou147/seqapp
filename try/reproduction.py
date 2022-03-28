#!/usr/bin/python3
from tkinter import E
import pandas as pd
import numpy as np
# Herceptin (trastuzumab) VH sequence information
Her_VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEW"+\
        "VARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVY"+\
        "YCSRWGGDGFYAMDYWGQGTLVTVSS"
Her_VH_mCDR3 = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQA"+\
        "PGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAY"+\
        "LQMNSLRAEDTAVYYGQGTLVTVSS"
Her_CDRH3 = "CSRWGGDGFYAMDYW"

# Herceptin (trastuzumab) VK sequence information
Her_VK = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLI"+\
        "YSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPP"+\
         "TFGQGTKVEIK"
Her_CDRL1 = "QDVNTA"
Her_CDRL3 = "QQHYTTPPT"

# pH for net charge calculation (Sharma et. al)
pH = 5.5


def net_charge(aa_seq):
  # The following function calculates the net charge of an input amino acid sequence
  # at a given pH.
  # aa_seq: an amino acid sequence

  # The following code reads in the pKa values taken from:
  # http://homepage.smc.edu/kline_peggy/Organic/Amino_Acid_pKa.pdf
  pKas = pd.read_csv("~/Documents/python/seqapp/try/data/pKas.csv",index_col = 0,header=None)

  # The following code reads in the Eisenberg hydrophobicity scales taken from:
  # Eisenberg et al. 1984
  Eisenberg = pd.read_csv("~/Documents/python/seqapp/try/data/Eisenberg.csv",index_col = 0,header=None)

  # Identify the positions of hydrophobic and hydrophilic residues
  phobic = [1,2,5,8,10,13,18,19,20]
  philic = [3,4,6,7,9,11,12,14,15,16,17]
  phobic = list(map(lambda x:x-1,phobic))
  philic = list(map(lambda x:x-1,philic))

  # Create integer list of amino acid counts
  AA_counts = list(map(aa_seq.count,pKas.index)) 
  
  # TODO: What does this calculate?
  
  dod = 1/(10**(pKas - pH) + 1)
  dod[1] = dod[1].map(lambda x: -1*x)
  dod[2] = dod[2].map(lambda x: 1-x)
  dod.iloc[[1,2,3,19],[2]] = -1*dod.iloc[[1,2,3,19],[2]]
  dod.iloc[[6,8,14],[2]] = 1-dod.iloc[[6,8,14],[2]]
  dod["row_sum"] = dod.apply(lambda x:x.sum(),axis =1)
  dod["row_sum"] = AA_counts*dod["row_sum"] 


  # Multiply integer list with WHATEVER DOD IS, and sum up to create net charge
  net_charge = sum(dod["row_sum"])
  return(net_charge)

# net_charge function test of one AA sequence
# print(net_charge('EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEW'))

def HIndex(aa_seq):
  # The following function calculates the hydrophobicity index of an input amino acid sequence.
  # aa_seq: an amino acid sequence

  # The following code reads in the pKa values taken from:
  # http://homepage.smc.edu/kline_peggy/Organic/Amino_Acid_pKa.pdf
  pKas = pd.read_csv("~/Documents/python/seqapp/try/data/pKas.csv",index_col = 0,header=None)

  # The following code reads in the Eisenberg hydrophobicity scales taken from:
  # Eisenberg et al. 1984
  Eisenberg = pd.read_csv("~/Documents/python/seqapp/try/data/Eisenberg.csv",index_col = 0,header=None)

  # Identify the positions of hydrophobic and hydrophilic residues
  phobic = [1,2,5,8,10,13,18,19,20]
  philic = [3,4,6,7,9,11,12,14,15,16,17]
  phobic = list(map(lambda x:x-1,phobic))
  philic = list(map(lambda x:x-1,philic))

  # Create integer list of amino acid counts
  AA_counts = list(map(aa_seq.count,Eisenberg.index))

  # Calculate the hydrophobicity index
  Eisenberg[1] =AA_counts*Eisenberg[1]
  phobic_sum = Eisenberg.iloc[phobic,[0]]
  philic_sum = Eisenberg.iloc[philic,[0]]
  HIndex = -phobic_sum[1].sum()/philic_sum[1].sum()
  
  return(HIndex)

# HIndex function test of one AA sequence
# print(HIndex("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEW")+HIndex(Her_CDRL1),HIndex(Her_CDRL3))



def run(x):

  # Calculate the net charge of one amino acid sequence and add the net charge of the entire VH sequence minus CDRH3
  NET_CHARGE = net_charge(x)

  VHNetCharge = net_charge(x) + net_charge(Her_VH_mCDR3)
  FabNetCharge = VHNetCharge + net_charge(Her_VK)
  FvCSP = VHNetCharge * net_charge(Her_VK)


  # Calculate the hydrophobicity index of one amino acid sequence
  CDH3_HI = HIndex(x)
  HISum = CDH3_HI + HIndex(Her_CDRL1) + HIndex(Her_CDRL3)

  return NET_CHARGE, CDH3_HI, FabNetCharge, FvCSP, HISum

