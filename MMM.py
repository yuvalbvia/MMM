import json
import numpy as np


class MMM:

    print("hi")
    sig_mat = np.load("/Users/sharonz/Downloads/data/BRCA-signatures.npy")
    json_data = open("/Users/sharonz/Downloads/data/ICGC-BRCA.json").read()
    data = json.loads(json_data)  # dictionary with data as: sample: chromosome#: "sequence": list of mutation#s


    def __init__(self, eMat, pi=None):
        pass

    def expectation(self):
        pass

    def maximization(self):
        pass

    def fit(data, threshold, max_iterations):
        pass

    def likelihood(data):
        pass

