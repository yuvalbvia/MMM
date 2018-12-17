import json
import numpy as np
from scipy.special import logsumexp


class MMM:

    sig_mat = np.load("/Users/sharonz/Downloads/data/BRCA-signatures.npy")
    json_data = open("/Users/sharonz/Downloads/data/ICGC-BRCA.json").read()
    data = json.loads(json_data)  # dictionary with data as: sample: chromosome#: "sequence": list of mutation#s


    def __init__(self, eMat, pi=None):
        pass

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

    def convergence(self,):
        theta0=


    def fit(data, threshold, max_iterations):
        pass

    def likelihood(data):
        pass

