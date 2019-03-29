#Usage: $ python predDTI_Enzyme.py Enzyme_Target.fasta Drug_selFV.csv > predDTI_Enzyme.out

import sys
import numpy as np
from sklearn.ensemble import RandomForestClassifier

protFile = sys.argv[1]
DrugFVfile = sys.argv[2]
selDipep = ['CN', 'CQ', 'WG', 'TH', 'GL', 'FL', 'GM', 'HM', 'PM', 'GF', 'VF', 'WS', 'DT', 'GT', 'RW', 'CW', 'GW', 'HW', 'WW']

protFV = []
protF = open(protFile, 'r')
allSeqVal = protF.read().rstrip().split('>')
protF.close()
#print allSeqVal

allSeqName = []
for p in allSeqVal[1:]:
  seq = ''
  seqVal = p.split('\n')
  allSeqName.append(seqVal[0])
  for seqFrag in seqVal[1:]:
    seq = seq + seqFrag
  #print seq
  tmpFV = []
  for dipep in selDipep:
    tmpFV.append(p.count(dipep))
  #print tmpFV
  tmpFV = [p1/400.0 for p1 in tmpFV]
  protFV.append(tmpFV)
#print protFV
#print len(allSeqName)
#print allSeqName


drugF = open(DrugFVfile, 'r')
allDrugVal = drugF.read().split('\n')[1:-1]
drugF.close()
#print allDrugVal

queryX = []
drugName = []
p2 = 0
while p2 < len(allSeqName):
  tmpQueryX = []
  for d2 in allDrugVal:
    tVal = d2.rstrip().split(',')
    if tVal[0] not in drugName:
      drugName.append(tVal[0])
    qX = tVal[1:] 
    qX = [float(j) for j in qX]
    tmpQueryX.append(protFV[p2]+ qX)
  queryX.append(tmpQueryX)
  p2 = p2 + 1
#print queryX
#print len(drugName)
#print drugName

RFpram = ['Fold1.csv', 'Fold2.csv', 'Fold3.csv', 'Fold4.csv', 'Fold5.csv']

if len(allSeqName) == len(queryX):
  predY = []
  for i in RFpram:
    nameF = i.split('_')[0]
    eVal = int(i.split('_')[1])

    inpF = open(nameF, 'r')
    allVal = inpF.read().split('\n')[1:-1]
    inpF.close()

    X = []
    Y = []
    for k in allVal:
      tmpVal = k.rstrip().split(',')[2:]
      tmpX = tmpVal[:-1]
      tmpX = [float(l) for l in tmpX]
      X.append(tmpX)

      if tmpVal[-1] == 'Interactive':
        Y.append(1)
      else:
        Y.append(-1)

    clf = RandomForestClassifier(n_estimators=eVal)
    clf.fit(np.array(X),np.array(Y))

    TpredY = []
    for q in queryX:
      tmpY = clf.predict(q)
      TpredY.append(tmpY)
    predY.append(TpredY)
  #print predY

  if len(predY[0]) == len(predY[1]) == len(predY[2]) == len(predY[3]) == len(predY[4]):
    n1 = 0
    while n1 < len(allSeqName):
      print 'Query:'
      print allSeqName[n1]
      print 'Interacts (potential) with the following drugs:'
      n2 = 0
      while n2 < len(drugName):
        #print n1, n2, predY[0][n1][n2], predY[1][n1][n2], predY[2][n1][n2], predY[3][n1][n2], predY[4][n1][n2]
        tmpPredY = [predY[0][n1][n2], predY[1][n1][n2], predY[2][n1][n2], predY[3][n1][n2], predY[4][n1][n2]]
        if tmpPredY.count(1) > 3:
          print drugName[n2], '[predicted', tmpPredY.count(1), 'times out of 5]'
        n2 += 1  
      print
      n1 += 1
