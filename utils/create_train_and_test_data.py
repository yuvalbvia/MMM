import json
import numpy as np


def get_train_test_data(raw_data, write_to_file=False):
    # returns tuple: train, test.
    # train is 90% of data, test is 10%
    # both are dictionaries in which: keys = sample names, values = array of mutations per sample
    train = dict()
    test = dict()
    samples_mut_matrix = dict()
    for sample in raw_data.keys():
        chrom_dict = raw_data[sample]
        sample_mut = []
        for chrom in chrom_dict.keys():
            mutations = chrom_dict[chrom]["Sequence"]
            sample_mut.extend(mutations)
        samples_mut_matrix[sample] = sample_mut

    for sample in samples_mut_matrix:
        a = samples_mut_matrix[sample]
        ninety_per = int(0.9 * len(a))
        np.random.shuffle(a)
        first = a[:ninety_per]
        last = a[ninety_per:]
        train[sample] = first
        test[sample] = last

    if write_to_file:
        json.dump(train, open("train_data.json", 'w'))
        json.dump(test, open("test_data.json", "w"))

    return train, test


if __name__ == '__main__':
    input = json.loads(open("data/ICGC-BRCA.json").read())
    train, test = get_train_test_data(input)
