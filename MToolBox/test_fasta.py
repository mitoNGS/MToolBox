#!/usr/bin/env python
# encoding=utf8

import os

def check_fasta(inhandle):
    """
		Check if a file is a fasta, returns Boolean (True/False)
		"""
    a = open(inhandle, 'r')
    i = a.readline()
    while i.strip() == "" and i != "":
        i = a.readline()
    return i.strip().startswith('>')


def fasta_list_maker(folder):
    fasta_list = []
    dirList = os.listdir(folder)
    for d in dirList:
        if os.path.isdir(os.path.join(folder, d)) == False:
            if check_fasta(os.path.join(folder,d)) == True:
                fasta_list.append(os.path.join(folder,d))
    return fasta_list

if __name__ == "__main__":
	ff = fasta_list_maker(os.getcwd())
	for f in ff:
		print os.path.split(f)[1]