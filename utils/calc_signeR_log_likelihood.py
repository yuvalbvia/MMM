import csv
import json

from MMM_many import MMM_many

f = open("/Users/yuvalbarzam/Documents/TAU/Sadna/MMM/data/ehat_signeR.csv")
data = csv.reader(f)
normalized_pi_matrix = []
firstline = True
for row in data:
    if firstline:  # skip first line
        firstline = False
        continue
    row = row[1:]
    norm = sum(map(float, row))
    newList = []
    for x in row:
        newList.append(float(x) / float(norm))
    normalized_pi_matrix.append(newList)
f.close()

f = open("/Users/yuvalbarzam/Documents/TAU/Sadna/MMM/data/phat_signeR.csv")
normalized_E_matrix = []
emat = csv.reader(f)
firstline=True
for row in emat:
    if firstline:
        firstline=False
        continue
    row = row[1:]
    # hacky but works - for some reason the reader reads one extra line
    if row[0] == "":
        break
    row = [float(x) for x in row]
    normalized_E_matrix.append(row)
f.close()

# calc likelihood for signeR results

MMM_instance = MMM_many(12, 560, 96, e_matrix_0=normalized_E_matrix, pi_matrix_0=normalized_pi_matrix, threshold=0.3)
input = json.loads(open("/Users/yuvalbarzam/Documents/TAU/Sadna/MMM/data/ICGC-BRCA.json").read())
mut_count_mat = MMM_instance.get_muation_counts_matrix(input)
likelihood_vector = MMM_instance.get_log_likelihood_vector_per_person(MMM_instance.pi__matrix_0,
                                                                      MMM_instance.e_matrix_0,
                                                                      mut_count_mat)

with open("signeR_log_likelihood_per_person.csv", "w") as w3:
    writer = csv.writer(w3, lineterminator='\n')
    for ll in likelihood_vector:
        writer.writerow([ll])
w3.close()
