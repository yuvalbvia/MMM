import csv
import json

f = open("MMM/data/ICGC-BRCA.json")
data = f.read()  # samples
d = json.loads(data)
f.close()
rows = d.keys()
columns = dict()
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
               '20', '21', '22', 'X', 'Y']

columns = [i for i in range(96)]

mat = list()
mat.append(columns)
names = list()

for item in d:
    names.append(item)
    mut_count = [0] * 96
    for ch in chromosomes:
        vec = d[item][ch]["Sequence"]
        for i in range(len(vec)):
            mut_count[vec[i]] += 1
    mut_count.insert(0, item)
    print(mut_count)
    mat.append(mut_count)

print(names)

with open("ICGC-BRCA.tsv", "w") as w:
    writer = csv.writer(w, delimiter="\t", lineterminator='\n')
    writer.writerows(mat)
w.close()
