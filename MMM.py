import json
import numpy as np
from scipy.special import logsumexp


class MMM:

    def __init__(self, e_matrix, data, pi_0=None):
        self.e_matrix = e_matrix
        self.data = get_clean_data(data)
        self.pi_0 = pi_0
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
        mutation_count = self.mutation_counts[mutation_num]  # Bj
        signature_prob = self.get_signature_prob_given_mutation(signature_num, mutation_count)  # P[Y=i|X=j]
        Eij = mutation_count * signature_prob

        return Eij

    def maximization(self, i):
        Ai = self.Ai_calculate(i)
        Ai_sum = 0
        for k in range(0, 12):
            Ai_sum += self.Ai_calculate(k)
        return Ai/Ai_sum  # pi_i

    def Ai_calculate(self, i):
        Ai = 0
        for j in range(0, 96):
            Ai += self.expectation(i, j)
        return Ai

    def fit(self, threshold, max_iterations, input):
        Eij_mat = []
        k = 0
        while k < max_iterations or :
            for i in range(0, 12):
                for j in range(0, 96):
                    Eij_mat[i][j] = self.expectation(i, j)


    def check_convergence(self, pi_1, pi_0):
        convergence = self.log_likelihood(len(self.data)) -


    def log_likelihood(self, x, pi):   # x is the data["input"] - a vector of mutation numbers
        sum = 0
        for t in range(0, x.len()):
            self.prob_for_mutation(x[t], pi)
        return sum

    def get_signature_prob_given_mutation(self, signature_num, mutation_num): # P[Yi|Xj]
        nominator = self.pi_0[signature_num] * self.e_matrix[signature_num][mutation_num]
        denominator = self.prob_for_mutation(mutation_num)
        result = nominator/denominator
        return result

    def prob_for_mutation(self, mutation_num, pi):
        ln_pi = np.log(pi)   #pi_1
        ln_Ej = np.log(self.e_matrix[:, [mutation_num]])   #e_1
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
    sig_mat = np.load("data/BRCA-signatures.npy")  # 12 signatures by 96 mutations numpy matrix
    json_data = open("data/ICGC-BRCA.json").read()
    data = json.loads(json_data)  # dictionary with data as: sample: chromosome#: "sequence": list of mutation#s
    input = json.loads(open("data/example.json").read())['input']
    threshold = 0
    max_iterations = 1000

    MMM_instance = MMM(sig_mat, data)
    MMM_instance.fit(threshold, max_iterations, input)

