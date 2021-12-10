# -*- coding: utf-8 -*-
import copy
import numpy as np

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

#print(D3_2(L3))

################# DAY 4 #######################
File4 = "InputD4.txt"
f4 = open(File4,"r")
S4 = f4.read().split("\n")
L4_nb = S4[0].split(",")
for i in range(len(L4_nb)):
    L4_nb[i] = int(L4_nb[i])
#print(L4_nb)
L4_mtrx = []
for i in range(1,len(S4)):
    if S4[i] != "":
        mtrx_line=[]
        for n in S4[i].split(" "):
            if n != "":
                mtrx_line.append(int(n))
        L4_mtrx[len(L4_mtrx) - 1].append(mtrx_line)
    else:
        L4_mtrx.append([])
#print(L4_mtrx[0])

def Check_bingo(Checked):
    for n in range(len(Checked)):
        for j in range(len(Checked[n][0])):
            Check = True
            for i in range(len(Checked[n])):
                if not(False in Checked[n][i]):
                    return(True,n)
                if Checked[n][i][j] and Check:
                    if i == len(Checked[n])-1:
                        return(True,n)
                else:
                    Check = False
    return (False,None)
#print(Check_bingo([[[True,False,False,False,False],[False,False,False,False,False],[True,False,False,False,False],[True,False,False,False,False],[True,False,False,False,False]]]))

def D4_1(L_mtrx, L_nb):
    Checked = []
    for i in L_mtrx:
        Checked.append([])
        for j in L_mtrx[0]:
            Checked[len(Checked)-1].append([False] * len(L_mtrx[0][0]))
    for nb in L_nb:
        for n in range(len(L_mtrx)):
            for i in range(len(L_mtrx[n])):
                for j in range(len(L_mtrx[n][i])):
                    if not(Checked[n][i][j]) and L_mtrx[n][i][j] == nb:
                        Checked[n][i][j] = True
        Check = Check_bingo(Checked)
        if Check[0]:
            mtrx = L_mtrx[Check[1]]
            sum = 0
            for i in range(len(mtrx)):
                for j in range(len(mtrx[i])):
                    if not(Checked[Check[1]][i][j]):
                        sum += mtrx[i][j]
            return(sum, nb, sum * nb, mtrx)
    return('no bingo winner')

def D4_2(L_mtrx,L_nb):
    notWinned = copy.deepcopy(L_mtrx)
    while len(notWinned) > 1:
        notWinned.remove(D4_1(notWinned,L_nb)[3])
        if [] in notWinned:
            notWinned.remove([])
    return (D4_1(notWinned,L_nb))

#print(D4_2(L4_mtrx,L4_nb))

############### DAY 5 ################@@

File5 = "InputD5.txt"
#File5="test.txt"
f5 = open(File5,"r")
S5 = f5.read().split("\n")
L5 = []
for line in S5:
    if line == "":
        break
    s = line.split(" -> ")
    L5.append([])
    for nb in s:
        nb = nb.split(",")
        L5[len(L5)-1].append( ( int(nb[0]), int(nb[1]) ) )

def visualise(L):
    s=""
    for line in L:
        for n in line:
            if n == 0:
                s+="."
            else:
                s+=str(n)
        s+="\n"
    return s


def D5_12(L):
    overlap = []
    if File5 == "test.txt" : N = 10
    else : N = 1000
    for i in range(N):
        overlap.append([0] * N)
    for line in L:
        if line[0][0] == line[1][0]:
            jmin = min(line[0][1], line[1][1])
            jmax = max(line[0][1], line[1][1])
            i = line[0][0]
            for j in range(jmin, jmax + 1):
                overlap[j][i] += 1
        elif line[0][1] == line[1][1]:
            imin = min(line[0][0], line[1][0])
            imax = max(line[0][0], line[1][0])
            j = line[1][1]
            for i in range(imin, imax + 1):
                overlap[j][i] += 1
        else:
            x = line[0][0]
            y = line[0][1]
            compare = line[1][0] + 1
            while x != compare:
                overlap[y][x] += 1
                if line[0][0] < line[1][0]:
                    x+=1
                    compare = line[1][0] + 1
                else:
                    x-=1
                    compare = line[1][0] - 1
                if line[0][1] < line[1][1]:
                    y+=1
                else:
                    y-=1
        #print(line)
        #print(visualise(overlap))

    c = 0
    for line in overlap:
        for nb in line:
            if nb > 1:
                c+=1
    #print(visualise(overlap))
    return c

#print(D5_12(L5))


############## DAY 6 ##############
File6 = "InputD6.txt"
# File6="test.txt"
f6 = open(File6,"r")
S6 = f6.read().split(",")
L6 = []
for s in S6:
    L6.append(int(s))

def D6_1(Fishs, Days):
    for d in range(Days):
        print(d)
        for i in range(len(Fishs)):
            if Fishs[i] == 0:
                Fishs.append(8)
                Fishs[i] = 6
            else:
                Fishs[i] -= 1
        #print(Fishs)
    return len(Fishs)

def D6_2(Fishs,Days):
    maxDays = 9
    L = [0] * maxDays
    for i in range(len(Fishs)):
        L[Fishs[i]] += 1
    for day in range(Days):
        L[(7 + day) % maxDays] += L[day % maxDays]

    # for i in range(len(Fishs)):
    #     Fishs[i] = ( Fishs[i] + Days%7 - 1) % 7
    return (L,sum(L))

#print(D6_2(L6,256))

################### DAY 7 #################

File7 = "InputD7.txt"
# File7="test.txt"
f7 = open(File7,"r")
S7 = f7.read().split(",")
L7 = []
for s in S7:
    L7.append(int(s))

def D7_1(L):
    L.sort()
    if len(L) < 1:
        return None
    if len(L) % 2 == 0 :
        M = L[(len(L) // 2 ) + 1 - 1]
    else:
        M = L[(len(L)-1)//2 - 1]

    sum = 0
    for n in L:
        sum += abs(M - n)
    return sum

def D7_2(L):
    Avrg1 = int(sum(L)/len(L) - 0.5)
    Avrg2 = int(sum(L)/len(L) + 0.5)

    Sum1 = 0
    Sum2 = 0
    for n in L:
        N1 = abs(Avrg1 - n)
        N2 = abs(Avrg2 - n)
        Sum1 += ( N1 * (N1+1) ) // 2
        Sum2 += ( N2 * (N2+1) ) // 2

    return min(Sum1,Sum2)
# print(D7_2(L7))

####################### DAY 8 ##################
File8 = "InputD8.txt"
# File8="test.txt"
f8 = open(File8,"r")
S8 = f8.read().split("\n")
L8 = []
for entry in S8:
    if entry != "":
        L8.append(["",""])
        line = entry.split("|")
        L8[len(L8) - 1][0] = line[0].split(" ")
        L8[len(L8) - 1][1] = line[1].split(" ")
while L8.count("") > 0:
    L8.remove("")
for line in L8:
    while line[0].count("") > 0:
        line[0].remove("")
    while line[1].count("") > 0:
        line[1].remove("")

def Check1478(digit):
    if len(digit) == 2 or len(digit) == 3 or len(digit) == 4 or len(digit) == 7:
        return True
    return False

def DecodeOutput(line):
    code = line[0][:]
    D = {}
    order = [""]*7 # haut, hautGauche, hautDroite, milieu, bas gauche, bas droite, bas

    #Etape 1
    x0 = ""
    x2x5 = ""
    x1x3 =""
    x4x6 = ""
    while x0 == "":
        for digit in code:
            if len(digit) == 3:
                Num = digit
            elif len(digit) == 2:
                Den = digit
        for i in Num:
            if not(i in Den):
                x0 = i
                order[0] = i
            else:
                x2x5 += i
    while x1x3 == "":
        for digit in code:
            if len(digit) == 4:
                for i in digit:
                    if not(i in x2x5):
                        x1x3 += i
    while x4x6 == "":
        for digit in code:
            if len(digit) == 7:
                for i in digit:
                    if not(i in x2x5 or i in x0 or i in x1x3):
                        x4x6 += i
    #Etape 2
    x2x3x4=""
    les3_6=[]
    while x2x3x4 == "":
        for digit in code:
            if len(digit) == 6:
                les3_6.append(digit)
        for i in "abcdefg":
            if not(i in les3_6[0] and i in les3_6[1] and i in les3_6[2]):
                x2x3x4 += i

    for i in x2x3x4:
        if i in x2x5:
            order[2] = i
            order[5] = x2x5[x2x5.index(i) - 1] #if is 1=> 0, if is 0=>-1
        elif i in x1x3:
            order[3] = i
            order[1] = x1x3[x1x3.index(i) - 1]
        elif i in x4x6:
            order[4] = i
            order[6] = x4x6[x4x6.index(i) - 1]
    return order




def D8_1(L):
    c = 0
    for line in L:
        output = line[1]
        for digit in output:
            if Check1478(digit):
                c += 1
    return c

def D8_2(L):
    c = 0
    Numbers = []
    for line in L:
        decode = DecodeOutput(line)
        Nb = [
        "".join(sorted(decode[0]+decode[1]+decode[2]+decode[4]+decode[5]+decode[6])),
        "".join(sorted(decode[2]+decode[5])),
        "".join(sorted(decode[0]+decode[2]+decode[3]+decode[4]+decode[6])),
        "".join(sorted(decode[0]+decode[2]+decode[3]+decode[5]+decode[6])),
        "".join(sorted(decode[1]+decode[2]+decode[3]+decode[5])),
        "".join(sorted(decode[0]+decode[1]+decode[3]+decode[5]+decode[6])),
        "".join(sorted(decode[0]+decode[1]+decode[3]+decode[4]+decode[5]+decode[6])),
        "".join(sorted(decode[0]+decode[2]+decode[5])),
        "".join(sorted(decode[0]+decode[1]+decode[2]+decode[3]+decode[4]+decode[5]+decode[6])),
        "".join(sorted(decode[0]+decode[1]+decode[2]+decode[3]+decode[5]+decode[6]))]
        output = line[1]
        for i in range(len(output)):
            output[i] = "".join(sorted(output[i]))
        Code = {}
        for digit in line[0]:
            for i in range(len(Nb)):
                if Nb[i] == "".join(sorted(digit)):
                    Code["".join(sorted(digit))] = i
        Numbers.append(int(str(Code[output[0]])+str(Code[output[1]])+str(Code[output[2]])+str(Code[output[3]])))

    return sum(Numbers)

# print(D8_2(L8))


############## DAY 9 #################

File9 = "InputD9.txt"
# File9="test.txt"
f9 = open(File9,"r")
S9 = f9.read().split("\n")
L9 = []
for line in S9:
    L9.append([])
    for i in line:
        L9[len(L9)-1].append(int(i))
while [] in L9:
    L9.remove([])

def TestNeighbors(L,i,j):
    if L[i][j] == 9:
        return 0
    elif L[i][j] == 0:
        return 1
    elif i!=len(L)-1 and i!=0:
        if j != len(L[0]) - 1 and j != 0 and L[i-1][j] > L[i][j] and L[i+1][j] > L[i][j] and L[i][j-1] > L[i][j] and L[i][j+1] > L[i][j]:
                return(L[i][j] + 1)
        elif j == len(L[0]) - 1 and L[i-1][j] > L[i][j] and L[i+1][j] > L[i][j] and L[i][j-1] > L[i][j]:
            return (L[i][j] + 1)
        elif j == 0 and L[i-1][j] > L[i][j] and L[i+1][j] > L[i][j] and L[i][j+1] > L[i][j]:
            return(L[i][j] + 1)
    elif i == len(L) - 1:
        if j != len(L[0]) - 1 and j != 0 and L[i-1][j] > L[i][j] and L[i][j-1] > L[i][j] and L[i][j+1] > L[i][j]:
            return(L[i][j] + 1)
        elif j == len(L[0]) - 1 and L[i-1][j] > L[i][j] and L[i][j-1] > L[i][j]:
            return(L[i][j]+1)
        elif j == 0 and L[i-1][j] > L[i][j] and L[i][j+1] > L[i][j]:
            return(L[i][j]+1)
    elif i == 0:
        if j != len(L[0]) - 1 and j != 0 and L[i+1][j] > L[i][j] and L[i][j-1] > L[i][j] and L[i][j+1] > L[i][j]:
            return(L[i][j] + 1)
        elif j == len(L[0]) - 1 and L[i+1][j] > L[i][j] and L[i][j-1] > L[i][j]:
            return(L[i][j]+1)
        elif j == 0 and L[i+1][j] > L[i][j] and L[i][j+1] > L[i][j]:
            return(L[i][j]+1)
    return 0

def D9_1(L):
    c = 0
    for i in range(len(L)):
        for j in range(len(L[i])):
            c += TestNeighbors(L,i,j)
    return c

def D9_2(L):
    basins = []
    for i in range(len(L)):
        basins.append([False]*len(L[0]))
    Added = copy.deepcopy(basins)
    for i in range(len(L)):
        for j in range(len(L[0])):
            if not(L[i][j] == 9):
                basins[i][j] = True
    AdjPos=[(0,0),(0,-1),(-1,0),(0,1),(1,0)]
    NbBasins = []
    for i in range(len(basins)):
        for j in range(len(basins[i])):
            ToTrack=[]
            if basins[i][j] and not(Added[i][j]):
                NbBasins.append(0)
                ToTrack = [(i,j)]
            while len(ToTrack)>0:
                # print(ToTrack)
                Tracking = ToTrack[len(ToTrack)-1]
                for pos in AdjPos:
                    try:
                        TrackAfter = tuple(np.add(Tracking,pos))
                        if -1 not in TrackAfter and not(Added[TrackAfter[0]][TrackAfter[1]]) and basins[TrackAfter[0]][TrackAfter[1]]:
                            ToTrack.append(TrackAfter)
                    except:
                        pass
                if not(Added[Tracking[0]][Tracking[1]]):
                    NbBasins[len(NbBasins)-1] += 1
                    Added[Tracking[0]][Tracking[1]] = True
                ToTrack.remove(Tracking)
    Maxis_3 = []
    while len(Maxis_3) !=3:
        M = max(NbBasins)
        Maxis_3.append(M)
        NbBasins.remove(M)
    return(Maxis_3[0]*Maxis_3[1]*Maxis_3[2])

# print(D9_2(L9))

################# DAY 10 ############

File10 = "InputD10.txt"
# File10="test.txt"
f10 = open(File10,"r")
S10 = f10.read()
L10 = S10.split("\n")

def CheckLine(line):
    opened = []
    C_dict = {')' : 3, ']' : 57, '}' : 1197, '>' : 25137}
    Chunk_dict = {'(' : ')', '[' : ']', '{' : '}', '<' : '>'}
    for s in line:
        if s in "({[<":
            opened.append(s)
        elif s != Chunk_dict[opened[len(opened)-1]]:
            return C_dict[s]
            break
        else:
            opened = opened[:len(opened)-1]
    return opened

def D10_1(L):
    c = 0
    for line in L:
        try:
            c += CheckLine(line)
        except:
            pass
    return c

def D10_2(L):
    Scores = [0]
    Scores_dict = {'(' : 1, '[' : 2, '{' : 3, '<' : 4}
    for line in L:
        s_added = ""
        opened = CheckLine(line)
        if Scores[-1] != 0:
            Scores.append(0)
        if type(opened) == int :
            # print(opened)
            pass
        else :
            s_added += "".join(opened)
            #print(s_added)
            for s in s_added[::-1]:
                Scores[-1] *= 5
                Scores[-1] += Scores_dict[s]
    if 0 in Scores:
        Scores.remove(0)
    return sorted(Scores)[len(Scores) // 2]

print(D10_2(L10))
