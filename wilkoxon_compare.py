import csv
from scipy.stats import wilcoxon

if __name__ == '__main__':
    f = open("/Users/yuvalbarzam/Documents/TAU/Sadna/MMM_log_likelihood_per_person.csv")
    data = csv.reader(f)
    MMM_likelihood = []
    for row in data:
        MMM_likelihood.append(float(row[0]))
    f2 = open("/Users/yuvalbarzam/Documents/TAU/Sadna/signeR_log_likelihood_per_person.csv")
    data2 = csv.reader(f2)
    signeR_likelihood = []
    for row in data2:
        signeR_likelihood.append(float(row[0]))
    f.close()
    f2.close()

    w = wilcoxon(MMM_likelihood, signeR_likelihood)
    print(w)
