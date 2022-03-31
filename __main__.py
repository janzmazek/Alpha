from metabolism.metabolism import Alpha
from secretion.alpha.secretion import run_montefusco_parallel
import numpy as np
import matplotlib.pyplot as plt
import pathlib


G = np.arange(0, 10.5, 0.5)
AA =[0, 0.5, 1]

def plot_results(amino_acids):
    fig, ax = plt.subplots()
    for AA in amino_acids:
        results_path = pathlib.Path(__file__).parent / "results" / f"{int(AA*100):03}%AA.txt"
        results = np.loadtxt(results_path).T
        ax.plot(results[0], results[2])
    plt.show()

def main():
    alpha = Alpha()

    try:
        plot_results(AA)
    except IOError as e:
        print(e)
        for a, res in run_montefusco_parallel(alpha, G, AA):
            results_path = pathlib.Path(__file__).parent / "results" / f"{int(a*100):03}%AA.txt"
            np.savetxt(results_path, res)
        plot_results(AA)


if __name__ == "__main__":
    main()