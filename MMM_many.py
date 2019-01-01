import json
import numpy as np
from scipy.special import logsumexp


class MMM_many:
    """
    This class will perform the MMM algorithm for an in put of many people, using the previous MMM
    """

    def __init__(self, sig_num, num_of_people, num_of_mut, e_matrix, pi):
        self.sig_num = sig_num
        self.num_of_people = num_of_people
        self.num_of_mut = num_of_mut
        self.e_matrix = e_matrix
        self.pi = pi

    def expectation(self):
        pass

    def maximization(self):
        pass

    def fit(self):
        pass

    def log_likelihood(self):
        pass

    def get_mutation_count_np_array(self, data):
        mutations = np.zeros(self.num_of_mut)
        for i in range(len(data)):
            for j in data[i + 1]:
                mutations[j] += 1
        return mutations

    def get_muation_counts_matrix(self, raw_data):  # removes "Sequence" and sample names from the input
        mut_count_mat = np.zeros((len(raw_data), self.num_of_mut))
        j = 0
        for sample in raw_data.keys():
            chromosomes_dict = raw_data[sample]
            clean_chromosomes_dict = self.get_clean_chromosome_dict(chromosomes_dict)
            mutation_counts = self.get_mutation_count_np_array(clean_chromosomes_dict)
            mut_count_mat[j] = mutation_counts
            j += 1
        return mut_count_mat

    def get_clean_chromosome_dict(self, raw_chromosomes_dict):
        new_chromosones_dict = dict()
        for i in raw_chromosomes_dict.keys():
            if i == 'X':
                new_chromosones_dict[21] = raw_chromosomes_dict[i]["Sequence"]
            elif i == 'Y':
                new_chromosones_dict[22] = raw_chromosomes_dict[i]["Sequence"]
            else:
                new_chromosones_dict[int(i)] = raw_chromosomes_dict[i]["Sequence"]
        return new_chromosones_dict


if __name__ == '__main__':
    input = json.loads(open("data/ICGC-BRCA.json").read())
    MMM_instance = MMM_many(12, len(input), 96, None, None)
    result = MMM_instance.get_muation_counts_matrix(input)
    # result_str = json.dumps(result.tolist())
    # print(type(result_str))
    # file = open('mutation_counts_matrix.txt', 'w')
    # file.write(result_str)
    # file.close()


