import json
import numpy as np
from scipy.special import logsumexp


class MMM_many:
    """
    This class will perform the MMM algorithm for an in put of many people, using the previous MMM
    """

    def __init__(self, sig_num, num_of_people, e_matrix, pi):
        self.sig_num = sig_num
        self.num_of_people = num_of_people
        self.e_matrix = e_matrix
        self.pi = pi


if __name__ == '__main__':
    input = json.loads(open("data/ICGC-BRCA.json").read())


