import json
import numpy as np
from scipy.special import logsumexp


class MMM:

    def __init__(self, e_matrix, data, pi=None):
        self.e_matrix = e_matrix
        self.data = get_clean_data(data)
        self.pi = pi
        self.mutation_counts = get_mutation_count_dict(data)  # Bj vector

    def get_random_signature_probs(self):
        sigs = 12*[0]
        for i in range(12):
            sigs[i] = np.random.randint(0,100)

        tmp = sum(sigs)

        for i in range(12):
            sigs[i] = sigs[i]/tmp

        return sigs

    def expectation(self, signature_num, mutation_num):
        mutation_count = self.mutation_counts[mutation_num]
        signature_prob = self.get_signature_prob_given_mutation(signature_num, mutation_count)
        return mutation_count * signature_prob

    def maximization(self, i):
        Ai = self.Ai_calculate(self, i)
        Ai_sum =0
        for k in range (0,12):
            Ai_sum += self.Ai_calculate(self, k)
        return Ai/Ai_sum

    def Ai_calculate(self, i):
        Ai=0
        for j in range (0,96):
            Ai += self.expectation(self, i, j)
        return Ai
    
    def calculate_log_prob_X_theta(self, x, theta):
        sum=0
        mul=1
        for t in range (0,x.len()):
            prob_for_mutation(self,x[t])
        return sum

    def fit(self, threshold, max_iterations):
        pass

    def likelihood(self):
        pass

    def get_signature_prob_given_mutation(self, signature_num, mutation_num):
        nominator = self.pi[signature_num] * self.e_matrix[signature_num][mutation_num]
        denominator = self.prob_for_mutation(mutation_num)
        result = nominator/denominator
        return result

    def prob_for_mutation(self, mutation_num):
        ln_pi = np.log(self.pi)
        ln_Ej = np.log(self.e_matrix[:, [mutation_num]])
        sum_vector = ln_pi + ln_Ej
        denominator = logsumexp(sum_vector)
        return denominator


def get_clean_data(data):  # removes "Sequence" from the input
    for i in data.keys():
        for j in data[i]:
            data[i][j] = data[i][j]["Sequence"]
    return data


def get_mutation_count_dict(data):
    mutations = dict()
    for i in data.keys():
        for j in data[i]:
            for k in data[i][j]:
                if k in mutations.keys():
                    mutations[k] += 1
                else:
                    mutations[k] = 1
    return mutations


if __name__ == '__main__':
    sig_mat = np.load("data/BRCA-signatures.npy")  # 12 by 96 matrix. 12 signatures by 96 mutations in each signature
    json_data = open("data/ICGC-BRCA.json").read()
    data = json.loads(json_data)  # dictionary with data as: sample: chromosome#: "sequence": list of mutation#s

    a = np.arange(2)
    print(a)
    print(logsumexp(a))
    #
    # MMM_instance = MMM(sig_mat, data)
    # MMM_instance.expectation()
    # MMM_instance.maximization()


