import csv
import json

import numpy as np
from scipy.special import logsumexp

from MMM import MMM


class MMM_many:
    """
    This class will perform the MMM algorithm for an in put of many people, using the previous MMM
    """

    def __init__(self, sig_num, num_of_ppl, num_of_mut, e_matrix_0, pi_matrix_0, threshold):
        self.sig_num = sig_num
        self.num_of_ppl = num_of_ppl
        self.num_of_mut = num_of_mut
        self.e_matrix_0 = np.log(e_matrix_0)
        self.pi__matrix_0 = np.log(pi_matrix_0)
        self.threshold = threshold

    def expectation_many(self, mut_count_matrix):
        # return tuple: (logE0 = sum of all ei, and logA matrix of logAi for each person)
        log_Emat_3d = np.zeros((self.num_of_ppl, self.sig_num, self.num_of_mut))
        log_A_mat = np.zeros((self.num_of_ppl, self.sig_num))
        for i in range(self.num_of_ppl):
            MMM_i = MMM(np.exp(self.e_matrix_0), np.exp(self.pi__matrix_0[i]), self.threshold)
            exp_res = MMM_i.expectation(mut_count_matrix[i])
            log_Emat_3d[i] = exp_res[0]
            log_A_mat[i] = exp_res[1]

        log_Emat = logsumexp(log_Emat_3d, axis=0)
        return log_Emat, log_A_mat

    def maximization_many(self, expectation_res, is_normalized=False):
        # returns tuple: (logE1 = new E, and logPi_1 matrix of pi_1 for each person)
        log_Emat_1 = expectation_res[0].T - logsumexp(expectation_res[0], axis=1)
        log_Emat_1 = log_Emat_1.T

        if is_normalized:
            log_pi_1_mat = expectation_res[1]
        else:
            log_pi_1_mat = expectation_res[1] - logsumexp(expectation_res[1])

        return log_Emat_1, log_pi_1_mat

    def fit_many(self, max_iter, data, normalize=False):
        mut_count_mat = self.get_muation_counts_matrix(data)
        mut_count_mat_copy = mut_count_mat.copy()
        if normalize:
            mut_count_mat = mut_count_mat.T / mut_count_mat.sum(axis=1).T
            mut_count_mat = mut_count_mat.T
        k = 0
        expectation_res = self.expectation_many(mut_count_mat)
        maximization_res = self.maximization_many(expectation_res, is_normalized=normalize)
        e_matrix_1 = maximization_res[0]
        pi_matrix_1 = maximization_res[1]

        print("Starting fitting process:")

        while k < max_iter:
            print("Iteration number:{!s}".format(k))
            if k == 100:
                ll_0 = self.log_likelihood(self.pi__matrix_0, self.e_matrix_0, mut_count_mat_copy)
            elif k >= 100:
                ll_1 = self.log_likelihood(pi_matrix_1, e_matrix_1, mut_count_mat_copy)
                convergence = ll_1 - ll_0
                print("Convergence:{!s}".format(convergence))
                ll_0 = ll_1
                if convergence < self.threshold:
                    print("The final likelihood is:{!s}".format(ll_1))
                    break
            self.pi__matrix_0 = pi_matrix_1
            self.e_matrix_0 = e_matrix_1
            expectation_res = self.expectation_many(mut_count_mat)
            maximization_res = self.maximization_many(expectation_res, is_normalized=normalize)
            e_matrix_1 = maximization_res[0]
            pi_matrix_1 = maximization_res[1]
            k += 1
        self.pi__matrix_0 = pi_matrix_1
        self.e_matrix_0 = e_matrix_1

    def log_likelihood(self, pi_matrix, e_matrix, mut_count_matrix):
        sum_ll = 0
        for i in range(self.num_of_ppl):  # the length of pi_matrix
            MMM_i = MMM(np.exp(e_matrix), np.exp(pi_matrix[i]), self.threshold)
            sum_ll += MMM_i.log_likelihood(pi_matrix[i], mut_count_matrix[i])
        # to do - take a look if the following two lines are better than what we did
        # tmp = np.dot(np.exp(pi_matrix), np.exp(e_matrix)) //the probability that the MMM of the ith person shows the jth mutation
        # sum_ll = np.sum(np.log(tmp) * mut_count_matrix) // log(temp) = log likelihood
        return sum_ll

    def get_mutation_count_np_array(self, data):
        mutations = np.zeros(self.num_of_mut)
        for mut in data:
            mutations[mut] += 1
        return mutations

    def get_muation_counts_matrix(self, raw_data):
        # returns numpy matrix in which lines are people and columns are mutations
        mut_count_mat = np.zeros((len(raw_data), self.num_of_mut))
        j = 0
        for sample in raw_data.keys():
            chromosomes_list = raw_data[sample]
            mutation_counts = self.get_mutation_count_np_array(chromosomes_list)
            mut_count_mat[j] = mutation_counts
            j += 1
        return mut_count_mat

    def get_log_likelihood_vector_per_person(self, pi_matrix, e_matrix, mut_count_matrix):
        likelihood_vector = []
        for i in range(self.num_of_ppl):  # the length of pi_matrix
            MMM_i = MMM(np.exp(e_matrix), np.exp(pi_matrix[i]), self.threshold)
            log_likelihood = MMM_i.log_likelihood(pi_matrix[i], mut_count_matrix[i])
            likelihood_vector.append(log_likelihood)
        return likelihood_vector


def get_random_probs_mat(num_of_rows, num_of_cols):
    matrix = [num_of_cols * [0] for i in range(num_of_rows)]
    for i in range(num_of_rows):
        for j in range(num_of_cols):
            matrix[i][j] = np.random.randint(1, 100)
        tmp = sum(matrix[i])
        for j in range(num_of_cols):
            matrix[i][j] = matrix[i][j] / tmp
    return matrix


def run_MMM_with_cosmic(train, sigs_probs, max_iterations=10000, threshold=0.3):
    Emat_probs = np.load("data/BRCA-signatures.npy")

    MMM_instance = MMM_many(12, len(train), 96, Emat_probs, sigs_probs, threshold)
    MMM_instance.fit_many(max_iterations, train, normalize=True)

    handle_MMM_result(MMM_instance, "cosmic_pi_output", "cosmic_E_output", "cosmic_log_likelihood_per_person")


def run_MMM(train, sigs_probs, max_iterations=10000, threshold=0.3):
    Emat_probs = get_random_probs_mat(12, 96)

    MMM_instance = MMM_many(12, len(train), 96, Emat_probs, sigs_probs, threshold)
    MMM_instance.fit_many(max_iterations, train, normalize=True)

    handle_MMM_result(MMM_instance, "MMM_pi_output", "MMM_E_output", "MMM_log_likelihood_per_person")


def handle_MMM_result(MMM_instance, pi_output_filename, E_output_filename, log_likelihood_filename):
    result_pi = np.exp(MMM_instance.pi__matrix_0)
    with open("{!s}.csv".format(pi_output_filename), "w") as w:
        writer = csv.writer(w, lineterminator='\n')
        writer.writerows(result_pi)
    w.close()

    result_e = np.exp(MMM_instance.e_matrix_0)
    with open("{!s}.csv".format(E_output_filename), "w") as w2:
        writer = csv.writer(w2, lineterminator='\n')
        writer.writerows(result_e)
    w2.close()

    test_data = json.loads(open("data/test_data.json").read())
    mut_count_matrix = MMM_instance.get_muation_counts_matrix(test_data)
    likelihood_vector = MMM_instance.get_log_likelihood_vector_per_person(MMM_instance.pi__matrix_0,
                                                                          MMM_instance.e_matrix_0,
                                                                          mut_count_matrix)

    with open("{!s}.csv".format(log_likelihood_filename), "w") as w3:
        writer = csv.writer(w3, lineterminator='\n')
        for ll in likelihood_vector:
            writer.writerow([ll])
    w3.close()


if __name__ == '__main__':
    train = json.loads(open("data/train_data.json").read())
    sigs_probs = get_random_probs_mat(len(train), 12)
    run_MMM(train, sigs_probs, 10000, 0.3)
    run_MMM_with_cosmic(train, sigs_probs, 10000, 0.3)
