import json
import numpy as np
from scipy.special import logsumexp


class MMM:

    def __init__(self, e_matrix, pi_0, threshold):  # class has pi,e,threshold
        self.e_matrix = np.log(e_matrix)
        self.pi_0 = np.log(pi_0)
        self.threshold = threshold

    def expectation(self, mutation_counts):
        log_mut = np.log(mutation_counts)
        p_xy = self.pi_0 + self.e_matrix.T
        denominator = logsumexp(p_xy, axis=1)
        log_Emat = log_mut + p_xy.T - denominator
        logA = logsumexp(log_Emat, axis=1)
        expectation = (log_Emat, logA)
        return expectation

    def maximization(self, expectation_res):
        pi_1 = expectation_res[1] - logsumexp(expectation_res[1])
        return pi_1

    def fit(self, max_iterations, data):
        mutation_counts = self.get_mutation_count_np_array(data)
        k = 0
        expectation_res = self.expectation(mutation_counts)
        pi_1 = self.maximization(expectation_res)
        convergence = self.log_likelihood(pi_1, mutation_counts) - self.log_likelihood(self.pi_0, mutation_counts)

        while k < max_iterations and convergence >= self.threshold:
            self.pi_0 = pi_1
            expectation_res = self.expectation(mutation_counts)
            pi_1 = self.maximization(expectation_res)
            k += 1
            convergence = self.log_likelihood(pi_1, mutation_counts) - self.log_likelihood(self.pi_0, mutation_counts)
        self.pi_0 = pi_1

    def log_likelihood(self, pi, mutation_counts):   # x is the data["input"] - a vector of mutation numbers
        p_xy = pi + self.e_matrix.T
        result = np.sum(logsumexp(p_xy, axis=1) * mutation_counts)
        return result

    def get_mutation_count_np_array(self, data):
        mutations = np.zeros(96)
        for i in range(len(data)):
            for j in data[i + 1]:
                mutations[j] += 1
        return mutations


def get_clean_data(raw_data):  # removes "Sequence" from the input
    for i in raw_data.keys():
        raw_data[i] = raw_data[i]["Sequence"]
    return raw_data


def get_random_signature_probs():
    sigs = 12*[0]
    for i in range(12):
        sigs[i] = np.random.randint(1, 100)

    tmp = sum(sigs)

    for i in range(12):
        sigs[i] = sigs[i]/tmp

    return sigs


if __name__ == '__main__':
    sig_mat = np.load("data/BRCA-signatures.npy")  # 12 signatures by 96 mutations numpy matrix
    input = json.loads(open("data/example.json").read())['input']   # input vector of mutations appearances
    fixed_input = {1: {'Sequence': input}}
    data = get_clean_data(fixed_input)
    initial_pi = json.loads(open("data/example.json").read())['initial_pi']
    threshold = 0.001
    max_iterations = 1000

    if not initial_pi:
        initial_pi = get_random_signature_probs()

    MMM_instance = MMM(sig_mat, initial_pi, threshold)
    MMM_instance.fit(max_iterations, data)
    print(np.exp(MMM_instance.pi_0))  # print the new pi vector

