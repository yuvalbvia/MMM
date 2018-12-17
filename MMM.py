import json
import numpy as np



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

    def maximization(self):
        pass

    def fit(data, threshold, max_iterations):
        pass

    def likelihood(data):
        pass



