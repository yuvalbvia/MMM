import json
import numpy as np
from scipy.special import logsumexp



class MMM:

    sig_mat = np.load("data/BRCA-signatures.npy")
    json_data = open("data/ICGC-BRCA.json").read()
    data = json.loads(json_data)  # dictionary with data as: sample: chromosome#: "sequence": list of mutation#s

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

    def fit(data, threshold, max_iterations):
        pass

    def likelihood(data):
        pass



