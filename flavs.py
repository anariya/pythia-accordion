from random import random

flavs = [1, 2]
probs = [0.5, 0.5] # nt used yet

numHadrons = 19
numEvents = 10000000
numDDbar = 0
numUUbar = 0
numUDbar = 0
numDUbar = 0
numDDbarjoin = 0
numUUbarjoin = 0
numUDbarjoin = 0
numDUbarjoin = 0

posOld = 1
negOld = -1
firstPos = True
firstNeg = True

for _ in range(numEvents):
    for i in range(numHadrons - 1):
        fromPos = random() < 0.5
        if fromPos:
            posNew = -1 if (random() < 0.02) else -2
            if firstPos:
                firstPos = False
            elif (posNew == -1 and posOld == 1):
                numDDbar += 1
            elif (posNew == -2 and posOld == 2):
                numUUbar += 1
            elif (posNew == -1 and posOld == 2):
                numUDbar += 1
            else:
                numDUbar += 1
            posOld = -posNew
        else:
            negNew = 1 if (random() < 0.02) else 2
            if firstNeg:
                firstNeg = False
            elif (negNew == 1 and negOld == -1):
                numDDbar += 1
            elif (negNew == 2 and negOld == -2):
                numUUbar += 1
            elif (negNew == 1 and negOld == -2):
                numDUbar += 1
            else:
                numUDbar += 1
            negOld = -negNew

    if (negOld == -1 and posOld == 1):
        numDDbarjoin += 1
    elif (negOld == -2 and posOld == 2):
        numUUbarjoin += 1
    elif (negOld == -1 and posOld == 2):
        numUDbarjoin += 1
    else:
        numDUbarjoin += 1

totalreg = numDDbar + numUUbar + numUDbar + numDUbar
totaljoin = numDDbarjoin + numUUbarjoin + numUDbarjoin + numDUbarjoin

ratioDDbar = (numDDbarjoin / totaljoin) / (numDDbar / totalreg)
ratioUUbar = (numUUbarjoin / totaljoin) / (numUUbar / totalreg)
ratioUDbar = (numUDbarjoin / totaljoin) / (numUDbar / totalreg)
ratioDUbar = (numDUbarjoin / totaljoin) / (numDUbar / totalreg)

print(f"ratioDDbar = {ratioDDbar}")
print(f"ratioUUbar = {ratioUUbar}")
print(f"ratioUDbar = {ratioUDbar}")
print(f"ratioDUbar = {ratioDUbar}")

    

