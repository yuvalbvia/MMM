import csv
from scipy.stats import wilcoxon


def get_wilcoxon_result(filepath_1, filepath_2):
    f1 = open(filepath_1)
    data1 = csv.reader(f1)
    likelihood1 = []
    for row in data1:
        likelihood1.append(float(row[0]))

    f2 = open(filepath_2)
    data2 = csv.reader(f2)
    likelihood2 = []
    for row in data2:
        likelihood2.append(float(row[0]))

    f1.close()
    f2.close()

    w = wilcoxon(likelihood1, likelihood2)
    return w


if __name__ == '__main__':
    print("Using Wilcoxon to compare between log likelihoods of:")
    MMM_file_path = "./data/results/MMM_random/MMM_log_likelihood_per_person.csv"
    signeR_file_path = "./data/results/signeR/signeR_log_likelihood_per_person.csv"
    convex_file_path = "./data/results/convex/MMM_convex_log_likelihood_per_person.csv"
    cosmic_file_path = "./data/results/MMM_with_cosmic/cosmic_log_likelihood_per_person.csv"

    w = get_wilcoxon_result(MMM_file_path, signeR_file_path)
    print("Regular MMM and signeR: {}".format(w))

    w = get_wilcoxon_result(convex_file_path, signeR_file_path)
    print("Convex MMM and signeR: {}".format(w))

    w = get_wilcoxon_result(cosmic_file_path, MMM_file_path)
    print("MMM with cosmic matrix and MMM with random matrix: {}".format(w))


