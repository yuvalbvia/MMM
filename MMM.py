import json
import numpy as np
from numpy import ma
from scipy.special import logsumexp


class MMM:

    def __init__(self, e_matrix, data, pi_0):
        self.e_matrix = np.log(e_matrix)
        self.pi_0 = np.log(pi_0)
        self.data = get_clean_data(data)
        self.mutation_counts = get_mutation_count_dict(self.data)  # Bj vector
        self.threshold = threshold

    def expectation(self):
        log_mut = np.log(self.mutation_counts)
        p_xy = self.pi_0 + self.e_matrix.T
        denominator = logsumexp(p_xy, axis=1)
        log_Emat = log_mut + p_xy.T - denominator
        return log_Emat

    def maximization(self, log_Emat):
        logA = self.get_logA(log_Emat)
        pi_1 = logA - logsumexp(logA)
        return pi_1

    def get_logA(self, log_Emat):
        logA = logsumexp(log_Emat, axis=1)
        return logA

    def fit(self, threshold, max_iterations):
        k = 0
        pi_0 = self.pi_0
        pi_1 = pi_0
        while k == 0 or (k < max_iterations and not self.is_converged(pi_1, pi_0, threshold)):
            self.pi_0 = pi_1
            log_Emat = self.expectation()
            pi_1 = self.maximization(log_Emat)
            k += 1
        self.pi_0 = pi_1

    def is_converged(self, pi_1, pi_0, threshold):
        convergence = self.log_likelihood(pi_1) - self.log_likelihood(pi_0)
        return convergence < threshold

    def log_likelihood(self, pi):   # x is the data["input"] - a vector of mutation numbers
        p_xy = pi + self.e_matrix.T
        result = np.sum(logsumexp(p_xy, axis=1))
        return result


def get_clean_data(data):  # removes "Sequence" from the input
    for i in data.keys():
        data[i] = data[i]["Sequence"]
    return data


def get_mutation_count_dict(data):
    mutations = np.zeros(96)
    for i in range(len(data)):
        for j in data[i + 1]:
            mutations[j] += 1
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
    input = json.loads(open("data/example.json").read())['input']   # input vector of mutations appearances
    fixed_input = {1: {'Sequence': input}}
    initial_pi = json.loads(open("data/example.json").read())['initial_pi']
    threshold = 0.001
    max_iterations = 1000

    if not initial_pi:
        initial_pi = get_random_signature_probs()

    MMM_instance = MMM(sig_mat, fixed_input, initial_pi)
    MMM_instance.fit(threshold, max_iterations)
    print(np.exp(MMM_instance.pi_0))  # print the new pi vector

