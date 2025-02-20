"""
Protein Sequencing Project
Name:
Roll Number:
"""

from os import read
import numpy as np
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    with open(filename) as file1:
        txt = file1.read().splitlines()
    return ''.join(txt)


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    temp = []
    for i in range(startIndex,len(dna),3):
        temp.append(dna[i:i+3])
        if dna[i:i+3]== 'TAG' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TGA' :
            break
    rna = [string.replace('T' , 'U') for string in temp] 
    return rna


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f = open(filename, "r")
    dict1 = json.load(f)
    dict2 = dict((i,k) for k,v in dict1.items() for i in v)
    dict3 = { k.replace('T', 'U'): v for k, v in dict2.items() }
    return dict3


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein = []
    if codons[0] == 'AUG' :
        protein.append('Start')
    for i in range(1,len(codons)):
        if i == 'UAG' or i== 'UAA' or i == 'UGA' :
            protein.append('Stop')
        else:
            protein.append(codonD[codons[i]])
    return protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    txt = readFile(dnaFilename)
    dict1 = makeCodonDictionary(codonFilename)
    finallst = []
    count = 0
    i = 0
    while i < len(txt):
        if txt[i:i+3] == 'ATG' :
            seq = dnaToRna(txt, i)
            proteinlst = generateProtein(seq, dict1)
            finallst.append(proteinlst)
            i = i + 3 * len(seq)
            
        else :
            count = count + 1
            i = i + 1
    print("total bases: ", len(txt))
    print("unused bases: ", count)
    print("synthesized proteins length: ",len(finallst))
    return finallst


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    res = []
    for i in range(len(proteinList1)):
            if proteinList1[i] in proteinList2:
                if proteinList1[i] not in res:
                    res.append(proteinList1[i])
    return res


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    res = []
    for i in proteinList:
        for j in i:
            res.append(j)
    return res


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    resdict = {}
    for i in aaList:
        if i not in resdict:
            resdict[i] = 1
        else:
            resdict[i] +=1
    return resdict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    lst1, lst2 = combineProteins(proteinList1), combineProteins(proteinList2)
    dict1, dict2 = aminoAcidDictionary(lst1), aminoAcidDictionary(lst2)
    freq1, freq2 = {}, {}
    temp = []
    aminodiffs = []  
    for i in dict1:
        freq1[i] = dict1[i] / len(lst1)
        if i not in temp and i != 'Start' and i != 'Stop' :
            temp.append(i)
    for j in dict2:
        freq2[j] = dict2[j] / len(lst2)
        if j not in temp and j != 'Start' and j != 'Stop' :
            temp.append(j)
    for acid in temp:
        frequency1 = 0
        frequency2 = 0
        if acid in freq1:
            frequency1 = freq1[acid]
        if acid in freq2:
            frequency2 = freq2[acid]
        diff = frequency2 - frequency1
        if diff > cutoff or diff < -cutoff:
           difflst = [acid, frequency1, frequency2]
           aminodiffs.append(difflst)
    return aminodiffs


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences:")
    for i in commonalities:
        cmmnproteins = ""
        lst = i [1:(len(i)-1)]
        count = 0
        for j in lst:
            cmmnproteins += j
            count += 1
            if count != len(lst):
                cmmnproteins += "-"
        if len(cmmnproteins) != 0:
            print(cmmnproteins)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for i in differences:
        acid = i[0]
        freq1 = round(i[1]*100, 2)
        freq2 = round(i[2]*100, 2)
        print(str(acid) + " "+ str(freq1)+ " % in Seq1" +", "+  str(freq2) + " % in Seq2")
    return 


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    res = []
    plst1, plst2 = combineProteins(proteinList1), combineProteins(proteinList2)
    acids1, acids2 = aminoAcidDictionary(plst1), aminoAcidDictionary(plst2)
    for i in acids1:
        if i not in res:
            res.append(i)
    for i in acids2:
        if i not in res:
            res.append(i)
    res.sort()
    return res


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    lst1 = combineProteins(proteinList)
    dict1 = aminoAcidDictionary(lst1)
    res = []
    for i in labels:
        if i in dict1:
            res.append(dict1[i] / len(lst1))
        else:
            res.append(0)
    return res


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    w = 0.4
    xValues = np.arange(len(xLabels))
    plt.bar(xValues, freqList1, width=-w, align='edge', label= label1, edgecolor = edgeList)
    plt.bar(xValues, freqList2, width= w, align='edge', label= label2, edgecolor = edgeList)
    plt.xticks(ticks=list(range(len(xLabels))), labels=xLabels, rotation = "horizontal")
    plt.legend()
    plt.title("Comparison Of Frequencies")
    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    edgelst = []
    words = []
    for i in range(len(biggestDiffs)):
        words.append(biggestDiffs[i][0])
    for i in range(len(labels)):
        if labels[i] in words:
            edgelst.append('black')
        else:
            edgelst.append('white')
    
    return edgelst


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    humanproteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    eleproteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")
    cmmnproteins = commonProteins(humanproteins, eleproteins)
    diffs1 = findAminoAcidDifferences(humanproteins, eleproteins, 0.005)
    displayTextResults(cmmnproteins, diffs1)
    labels = makeAminoAcidLabels(humanproteins, eleproteins)
    f1 = setupChartData(labels, humanproteins)
    f2 = setupChartData(labels, eleproteins)
    edges1 = makeEdgeList(labels, diffs1)
    createChart(labels, f1, "Human" , f2, "Elephant", edgeList= edges1)
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    #test.testFindAminoAcidDifferences()
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()

    ## Uncomment these for Week 2 ##
   
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
   
