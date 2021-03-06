import csv
import json

import numpy as np
from scipy.special import logsumexp

from MMM import MMM


class MMM_convex:
    """
    This class will perform the MMM algorithm for an in put of many people, using the gradient descent technique:
    Optimizing the E and Pi matrices separately, consequently
    """

    def __init__(self, sig_num, num_of_ppl, num_of_mut, e_matrix_0, pi_matrix_0, thresh, thresh_e, thresh_pi):
        self.sig_num = sig_num
        self.num_of_ppl = num_of_ppl
        self.num_of_mut = num_of_mut
        self.e_matrix_0 = np.log(e_matrix_0)
        self.pi__matrix_0 = np.log(pi_matrix_0)
        self.thresh = thresh
        self.threshold_e = thresh_e
        self.threshold_pi = thresh_pi

    def expectation(self, mut_count_matrix):
        # return tuple: (logE0 = sum of all ei, and logA matrix of logAi for each person)
        log_Emat_3d = np.zeros((self.num_of_ppl, self.sig_num, self.num_of_mut))
        log_A_mat = np.zeros((self.num_of_ppl, self.sig_num))
        for i in range(self.num_of_ppl):
            MMM_i = MMM(np.exp(self.e_matrix_0), np.exp(self.pi__matrix_0[i]), self.thresh)
            exp_res = MMM_i.expectation(mut_count_matrix[i])
            log_Emat_3d[i] = exp_res[0]
            log_A_mat[i] = exp_res[1]

        log_Emat = logsumexp(log_Emat_3d, axis=0)
        return log_Emat, log_A_mat

    def maximize_e(self, log_Emat):
        log_Emat_1 = log_Emat.T - logsumexp(log_Emat, axis=1)
        log_Emat_1 = log_Emat_1.T

        return log_Emat_1

    def maximize_pi(self, log_A_mat, is_normalized=False):
        if is_normalized:
            log_pi_1_mat = log_A_mat
        else:
            log_pi_1_mat = log_A_mat - logsumexp(log_A_mat)

        return log_pi_1_mat

    def fit_convex(self, data, normalize=False):
        print("Starting the fitting process")
        mut_count_mat = self.get_muation_counts_matrix(data)
        mut_count_mat_copy = mut_count_mat.copy()
        if normalize:
            mut_count_mat = mut_count_mat.T / mut_count_mat.sum(axis=1).T
            mut_count_mat = mut_count_mat.T
        final_score = -1 * np.inf
        num_iter = 0
        while True:
            num_iter += 1
            print("Iteration number: {!s}".format(num_iter))
            score = -1 * np.inf

            print("Starting to maximize E matrix")
            while True:
                log_Emat, _ = self.expectation(mut_count_mat)
                e_matrix_1 = self.maximize_e(log_Emat)
                new_score = self.log_likelihood(self.pi__matrix_0, e_matrix_1, mut_count_mat_copy)
                convergence = new_score - score
                print("Covergence:{!s}".format(convergence))
                score = new_score
                self.e_matrix_0 = e_matrix_1
                if convergence < self.threshold_e:
                    break

            print("Current log likelihood: {!s}".format(score))
            print("Starting to maximize pi matrix")
            while True:
                _, log_A_mat = self.expectation(mut_count_mat)
                pi_matrix_1 = self.maximize_pi(log_A_mat, is_normalized=normalize)
                new_score = self.log_likelihood(pi_matrix_1, self.e_matrix_0, mut_count_mat_copy)
                convergence = new_score - score
                print("Covergence:{!s}".format(convergence))
                score = new_score
                self.pi__matrix_0 = pi_matrix_1
                if convergence < self.threshold_pi:
                    break

            print("Current log likelihood: {!s}".format(score))
            print("Finished maximizing both E and pi")
            new_final_score = self.log_likelihood(self.pi__matrix_0, self.e_matrix_0, mut_count_mat_copy)
            convergence = new_final_score - final_score
            final_score = new_final_score
            print("Current log likelihood for iteration {!s} is: {!s}".format(num_iter, final_score))
            if convergence < self.thresh:
                break

    def log_likelihood(self, pi_matrix, e_matrix, mut_count_matrix):
        sum_ll = 0
        for i in range(self.num_of_ppl):  # the length of pi_matrix
            MMM_i = MMM(np.exp(e_matrix), np.exp(pi_matrix[i]), self.thresh)
            sum_ll += MMM_i.log_likelihood(pi_matrix[i], mut_count_matrix[i])
        # to do - take a look if the following two lines are better than what we did
        # tmp = np.dot(np.exp(pi_matrix), np.exp(e_matrix))
        # sum_ll = np.sum(np.log(tmp) * mut_count_matrix)
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
            MMM_i = MMM(np.exp(e_matrix), np.exp(pi_matrix[i]), self.thresh)
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


def run_MMM_convex(train_data, thresh_e=0.001, thresh_pi=0.001, thresh=0.3):
    sigs_probs = get_random_probs_mat(len(train_data), 12)
    Emat_probs = get_random_probs_mat(12, 96)

    MMM_instance = MMM_convex(12, len(train_data), 96, Emat_probs, sigs_probs, thresh, thresh_e, thresh_pi)
    MMM_instance.fit_convex(train_data, normalize=True)

    handle_MMM_result(MMM_instance, "MMM_convex_pi_output", "MMM_convex_E_output",
                      "MMM_convex_log_likelihood_per_person")


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
    train_data = json.loads(open("data/train_data.json").read())
    run_MMM_convex(train_data)

