import re
import os
import matplotlib.pyplot as plt
import numpy as np


def estimate_marginal_importance(filename):
    likelihoods, priors, jacobians, variational_weights = read_vbis(filename)
    n = len(likelihoods)

    log_estimated_marginal_sum = 0.0
    log_samples = likelihoods + priors + jacobians - variational_weights
    log_estimated_marginal_sum = np.logaddexp.reduce(log_samples)

    log_estimated_marginal = log_estimated_marginal_sum - np.log(n)
    return log_estimated_marginal

def count_lines(filename):
    with open(filename, 'r', encoding="UTF-8") as file:
        count = sum(1 for _ in file)
    return count

def read_vbis(filename):
    with open(filename, 'r', encoding="UTF-8") as file:
        line_count = count_lines(filename)
        if line_count == 0:
            return None

        likelihood = np.zeros((line_count))
        prior = np.zeros((line_count))
        variational_weight = np.zeros((line_count))
        jacobian = np.zeros((line_count))

        for i, line in enumerate(file):
            match = re.search(r'&lnL=(-?\d+\.\d+)', line)
            if match:
                likelihood[i] = float(match.group(1))

            match = re.search(r'&lnPr=(-?\d+\.\d+)', line)
            if match:
                prior[i] = float(match.group(1))

            match = re.search(r'&lnQ=(-?\d+\.\d+)', line)
            if match:
                variational_weight[i] = float(match.group(1))

            match = re.search(r'&lnJac=(-?\d+\.\d+)', line)
            if match:
                jacobian[i] = float(match.group(1))

    nonzero_indices = np.nonzero(likelihood)

    return likelihood[nonzero_indices], prior[nonzero_indices], jacobian[nonzero_indices], variational_weight[nonzero_indices]


def compute_marginal_for_files(filenames):
    marginal_values = []

    for filename in filenames:
        if os.path.isfile(filename):
            marginal = estimate_marginal_importance(filename)
            if marginal is not None:
                marginal_values.append(marginal)
        else:
            marginal_values.append(None)
            # raise FileExistsError(filename)

    return marginal_values


def plot_over_booots():
    # plot the estimated marginal for different number of boosts
    directory = os.path.join("..", "analysis", "ds1", "vi", "up_nj", "")
    max_n = 10
    filenames = [os.path.join(directory, f"d3_lr1_i3_b{i}_jc69_exp_con/samples.t") for i in range(1,max_n+1)]
    boost_numbers = list(range(1, max_n+1))
    marginal_values = compute_marginal_for_files(filenames)

    for (i, marginal) in enumerate(marginal_values):
        print(f"Boosts{i+1}: -> marginal: {marginal}")

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(4, 4)
    ax.set_position([0.20, 0.11, 0.79, 0.88])
    plt.rcParams.update({'font.size': 11})
    plt.plot(boost_numbers, marginal_values, marker='o', color='k')
    plt.xlabel('Boosts (K)')
    plt.ylabel('Estimated Marginal Probability')
    plt.show()

    # path_save = os.path.join(".", "out", "boosts_marginal_exp_jc69_n500.pdf")
    # plt.savefig(path_save, format="pdf")


def print_over_ds():
    # print the estimated marginal for different data sets
    directory = os.path.join("..", "analysis")
    filenames = [os.path.join(directory, f"ds{i}", "vi", "up_nj", "d3_lr1_i1_b1_jc69_exp_ds8", "samples.t") for i in range(1, 9)]
    marginal_values = compute_marginal_for_files(filenames)
    for (i, marginal) in enumerate(marginal_values):
        print(f"DS{i+1}: -> marginal: {marginal}")

if __name__ == "__main__":
    print_over_ds()
