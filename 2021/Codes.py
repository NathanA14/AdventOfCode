# -*- coding: utf-8 -*-

############## DAY 1 ##############

File1 = "InputD1.txt"
f1 = open(File1,"r")
S1 = (f1.read()).split("\n")
L1 =[]
for s in S1:
    if s != "":
        L1.append(int(s))

def D1_1(L):
    c = 0
    for i in range(1,len(L)):
        if L[i] > L[i-1]:
            c+=1
    return (c)

def D1_2(L):
    L_sum = []
    for i in range (2,len(L)):
        L_sum.append(L[i-2] + L[i-1]+L[i])
    return L_sum

#print(D1_1(D1_2(L1)))

############### DAY 2 ###############

File2 = "InputD2.txt"
f2 = open(File2,"r")
S2 = (f2.read()).split("\n")
L2 = []
for s in S2:
    L2.append(s)

def D2_1(L):
    Dict = {"forward" : 1, "down" : 1, "up" : -1}
    c = {"forward" : 0, "down" : 0, "up" : 0}
    for instruction in L:
        if instruction != "":
            Inst = instruction.split(" ")
            c[Inst[0]] += Dict[Inst[0]] * int(Inst[1])
    return (c["forward"] * (c["down"] + c["up"]))

def D2_2(L):
    #L = ["forward 5", "down 5", "forward 8", "up 3", "down 8", "forward 2"]
    Dict = {"forward" : 1, "down" : 1, "up" : -1}
    c = {"aim" : 0, "depth" : 0, "forward" : 0}
    for instruction in L:
        if instruction != "":
            Inst = instruction.split(" ")
            if Inst[0] == "forward":
                c["forward"] += int(Inst[1])
                c["depth"] += c["aim"] * int(Inst[1])
            else:
                c["aim"] += Dict[Inst[0]] * int(Inst[1])
    return (c["depth"] * c["forward"])

#print (D2_2(L2))

############# DAY 3 ############
File3 = "InputD3.txt"
f3 = open(File3,"r")
S3 = f3.read().split("\n")
L3 = []
for i in range(len(S3)):
    if S3[i] == "":
        break
    L3.append([])
    for j in range(len(S3[i])):
        L3[i].append(int(S3[i][j]))

def bitStrToInt(S):
    I = 0
    for i in range(len(S)):
        I += int(S[i]) * (2 ** (len(S) - 1 - i))
    return (I)

def CountNb01(L):
    if len(L) <=1:
        return [[None],[None]]
    C = [[0] * len(L[1]),[0]*len(L[1])]
    for i in range(len(L)):
        for j in range(len(L[i])):
            C[L[i][j]][j] += 1
    return C

def D3_1(L):
    C = CountNb01(L)
    Gamma = Epsilon = ""
    for i in range(len(C[1])):
        Gamma += (str(bin(C[1][i]>C[0][i])).split("b")[1])
        Epsilon += (str(bin(C[1][i]<C[0][i])).split("b")[1])
    Gamma = bitStrToInt(Gamma)
    Epsilon = bitStrToInt(Epsilon)
    return (Gamma,Epsilon, Gamma * Epsilon)

def D3_2(L):
    O2 = L[:]
    CO2 = L[:]
    for j in range(len(L[0])):
        C_O2 = CountNb01(O2)
        C_CO2 = CountNb01(CO2)
        if len(O2) > 1:
            if C_O2[1][j] >= C_O2[0][j]:
                for b in O2[:]:
                    if b[j] == 0:
                        O2.remove(b)
            else:
                for b in O2[:]:
                    if b[j] == 1:
                        O2.remove(b)

        if len(CO2) > 1:
            if C_CO2[0][j] <= C_CO2[1][j]:
                for b in CO2[:]:
                    if b[j] == 1:
                        CO2.remove(b)
            else:
                for b in CO2[:]:
                    if b[j] == 0:
                        CO2.remove(b)
    NO2 = NCO2 = ""
    for i in range(len(O2[0])):
        NO2 += (str(bin(O2[0][i] == 1)).split("b")[1])
        NCO2 += (str(bin(CO2[0][i] == 1)).split("b")[1])
    NO2 = bitStrToInt(NO2)
    NCO2 = bitStrToInt(NCO2)

    return(O2,NO2,"\n",CO2,NCO2,"\n",NO2*NCO2)

print(D3_2(L3))
