import json
import numpy as np
from numpy import ma
from scipy.special import logsumexp


class MMM:

    def __init__(self, e_matrix, data, pi_0, input):
        self.e_matrix = e_matrix
        self.data = get_clean_data(data)
        self.pi_0 = pi_0
        self.mutation_counts = get_mutation_count_dict(data, input)  # Bj vector
        self.input = input

    def expectation(self, signature_num, mutation_num):
        mutation_count = self.mutation_counts[mutation_num]  # Bj
        signature_prob = self.get_signature_prob_given_mutation(signature_num, mutation_num)  # P[Y=i|X=j]
        Eij = mutation_count * signature_prob

        return Eij

    def maximization(self, i, Eij_mat):
        Ai = self.Ai_calculate(i, Eij_mat)
        Ai_sum = 0
        for k in range(0, 12):
            Ai_sum += self.Ai_calculate(k, Eij_mat)
        return Ai/Ai_sum  # pi_i

    def Ai_calculate(self, i, Eij_mat): # sums the i'th row in Eij_mat
        row = Eij_mat[i]
        return sum(row)

    def fit(self, threshold, max_iterations):
        Eij_mat = [[0] * 96 for i in range(12)]
        k = 0
        pi_vector = [0 for i in range(12)]
        pi_0 = self.pi_0
        while k < max_iterations and not self.is_converged(pi_vector, pi_0, threshold):
            for i in range(0, 12):
                for j in range(0, 96):
                    Eij_mat[i][j] = self.expectation(i, j)
            for i in range(0, 12):
                pi_vector[i] = self.maximization(i, Eij_mat)
            k += 1
            pi_0 = self.pi_0
            print(k)
            self.pi_0 = pi_vector

    def is_empty(self, pi):
        for num in pi:
            if num != 0:
                return False
        return True

    def is_converged(self, pi_1, pi_0, threshold):
        if self.is_empty(pi_1):
            return False
        convergence = self.log_likelihood(self.input, pi_1) - self.log_likelihood(self.input, pi_0)
        print("convergence")
        print(pi_0, pi_1)
        print(convergence)
        return convergence < threshold

    def log_likelihood(self, x, pi):   # x is the data["input"] - a vector of mutation numbers
        sum = 0
        for t in range(0, len(x)):
            sum += self.prob_for_mutation(x[t], pi)
        return sum

    def get_signature_prob_given_mutation(self, signature_num, mutation_num):  # P[Yi|Xj]
        nominator = self.pi_0[signature_num] * self.e_matrix[signature_num][mutation_num]
        denominator = self.prob_for_mutation(mutation_num, self.pi_0)
        result = nominator/denominator
        return result

    def prob_for_mutation(self, mutation_num, pi):

        ln_pi = np.ma.filled(np.log(ma.masked_equal(pi, 0)), 0)
        e = self.e_matrix[:, [mutation_num]]
        ln_Ej = np.ma.filled(np.log(ma.masked_equal(e, 0)), 0)
        sum_vector = ln_pi + ln_Ej
        prob = logsumexp(sum_vector)
        return prob


def get_clean_data(data):  # removes "Sequence" from the input
    for i in data.keys():
        for j in data[i]:
            data[i][j] = data[i][j]["Sequence"]
    return data


def get_mutation_count_dict(data, input):
    mutations = dict()
    for i in data.keys():
        for j in data[i]:
            for k in data[i][j]:
                if k in mutations.keys():
                    mutations[k] += 1
                else:
                    mutations[k] = 1

    for l in mutations:
        if l not in input:
            mutations.pop(l)
    return mutations


def get_random_signature_probs():
    sigs = 12*[0]
    for i in range(12):
        sigs[i] = np.random.randint(0,100)

    tmp = sum(sigs)

    for i in range(12):
        sigs[i] = sigs[i]/tmp

    return sigs


if __name__ == '__main__':
    sig_mat = np.load("data/BRCA-signatures.npy")  # 12 signatures by 96 mutations numpy matrix
    json_data = open("data/ICGC-BRCA.json").read()
    data = json.loads(json_data)  # dictionary with data as: sample: chromosome#: "sequence": list of mutation#s
    input = json.loads(open("data/example.json").read())['input']
    initial_pi = json.loads(open("data/example.json").read())['initial_pi']
    threshold = 0.001
    max_iterations = 1000

    MMM_instance = MMM(sig_mat, data, initial_pi, input)
    MMM_instance.fit(threshold, max_iterations)
    print(MMM_instance.pi_0) # print the new pi vector

