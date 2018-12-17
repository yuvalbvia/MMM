import json
import numpy as np
from scipy.special import logsumexp


class MMM:

    def __init__(self, e_matrix, data, pi=None):
        self.e_matrix = e_matrix
        self.data = clean_data(data)
        self.pi = pi
        self.mutation_counts = get_mutation_counts(data) #Bj

    def get_random_signature_probs(self):
        sigs = 12*[0]
        for i in range(12):
            sigs[i] = np.random.randint(0,100)

        tmp = sum(sigs)

        for i in range(12):
            sigs[i] = sigs[i]/tmp

        return sigs


    def get_mutation_count_dict(self):
        mutations = dict()
        for i in self.data.keys():
            for j in self.data[i]:
                for k in self.data[i][j]:
                    if k in mutations.keys():
                        mutations[k]+=1
                    else:
                        mutations[k]=1
        return mutations

    def clean_data(self): #removes "Sequence" from the input
        for i in self.data.keys():
            for j in self.data[i]:
                self.data[i][j]= self.data[i][j]["Sequence"]

    def expectation(self):
        pass

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

    def get_signature_prob_given_mutation(self, signature, mutation):
        nominator = self.pi[signature] * self.e_matrix[signature][mutation]


def clean_data(data):
    pass


def get_mutation_counts(data):
    pass


if __name__ == '__main__':
    sig_mat = np.load("data/BRCA-signatures.npy")
    json_data = open("data/ICGC-BRCA.json").read()
    data = json.loads(json_data)  # dictionary with data as: sample: chromosome#: "sequence": list of mutation#s

    MMM_instance = MMM(sig_mat, data)
    MMM_instance.expectation()
    MMM_instance.maximization()


