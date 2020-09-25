#reads fcl file into a dictionary
# Written by Josh LaBounty
def fclReader(fileName):
    fclDict = {}
    tempDicts = []
    tempKeys = []
    dictCounter = 0
    with open(str(fileName)) as file:
        for line in file:
            #print(line)
            #if((":" in line or "}" in line or "{" in line) and ("//" not in line)):
            if(("//" not in line) and ("PROLOG" not in line)):
                #print(line)
                if("{" in line):
                    #print("start of new nested dict")
                    newDict = {}
                    tempDicts.append(newDict)
                    key = line.split(":")[0].rstrip().lstrip()
                    #print(key)
                    tempKeys.append(key)
                    dictCounter += 1
                elif("}" in line):
                    #print("end of nested dict")
                    dictCounter -= 1
                    if(len(tempDicts) > 1):
                        thisKey = tempKeys[dictCounter]
                        previousKey = tempKeys[dictCounter - 1]
                        #print(thisKey, previousKey)
                        tempDicts[dictCounter - 1][thisKey] = tempDicts[dictCounter]
                        tempDicts = tempDicts[:-1]
                        tempKeys = tempKeys[:-1]
                    else:
                        thisKey = tempKeys[dictCounter]
                        fclDict[thisKey] = tempDicts[dictCounter]
                        tempDicts = tempDicts[:-1]
                        tempKeys = tempKeys[:-1]
                        
                elif(":" in line):
                    #print(line.split(":"))
                    key = line.split(":")[0].rstrip().lstrip()
                    value = line.split(":")[1].rstrip().lstrip()
                    try:
                        vi = float(value)
                        value = vi
                    except:
                        value = value
                    #print(key, value)
                    if(dictCounter > 0):
                        tempDicts[dictCounter-1][key] = value
                    else:
                        fclDict[key] = value
    return fclDict
def printDictInFcl(d, file, tabCounter = 0):
    tabString = ""
    for i in range(tabCounter):
        tabString+="\t"
    for di in d:
        di_val = d[di]
        if(isinstance(di_val,dict)):
            #print("nested")
            file.writelines([tabString,di," : {","\n"])
            printDictInFcl(di_val, file, tabCounter+1)
            file.writelines([tabString,"}","\n"])
        else:
            file.writelines([tabString,di," : ", str(di_val),"\n"])
#writes fcl to specified file
def fclWriter(fclDict, fileName, comment = ""):
    with open(fileName,"w") as file:
        file.writelines(["// FCL Generated from python dictionary","\n"])
        file.writelines(["// Comment: "+str(comment),"\n"])
        file.writelines(["BEGIN_PROLOG","\n"])
        printDictInFcl(fclDict, file)   
        file.writelines(["END_PROLOG","\n"])
    file.close()
