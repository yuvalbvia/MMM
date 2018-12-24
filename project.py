import csv
LOCATION = r"C:\Users\noabe\Desktop\Dudu burstein"
f = open(LOCATION+"\hmm_result_for_1M.tsv")
rawDB = f.readlines()
rawDB.pop(0)
rawDB.pop(1)

headers = ["target name", "accession", "query name","accession", "E-value (full sequence)", "score(full sequence)","Bias(full sequence)","E-value(best 1 domain)","score (best 1 domain)","bias (best 1 domain)","exp","reg","clu","ov","env","dom","rep","inc","description of target"]
DB = {}
for line in rawDB:
    parsed = line.split("\t")
    if parsed[0].startswith("ERR"):
        score = parsed[6]
        #eValue = parsed[5]
        protName = parsed[0]
        if DB.get(protName) != None:
            if DB[protName][6] < score:
                DB[protName] = parsed
        else:
            DB[protName] = parsed

with open(LOCATION+'\OUTPUT_bestscore.csv', 'w') as f:
    writer = csv.writer(f,lineterminator='\n')
    writer.writerow([g for g in headers])
    writer.writerows(DB.values())
