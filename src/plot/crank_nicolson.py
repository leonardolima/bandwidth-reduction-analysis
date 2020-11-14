#!/usr/bin python

import matplotlib.pyplot as plt
import pandas as pd
import itertools

def plot_runtime():
    N = [150, 222, 315, 416, 540, 670, 825, 984, 1170, 1358, 1575]
    before = [15, 42, 106, 231, 491, 919, 1690, 2839, 4750, 7406, 11338]
    after = [12, 29, 56, 103, 170, 272, 409, 604, 858, 1185, 1602]

    df = pd.DataFrame({'N': N, 'before': before, 'after': after})

    ax = plt.gca()

    df.plot(kind='line', x='N', y='before', label='Execution time', ax=ax, color="black", style='-')
    df.plot(kind='line', x='N', y='after', label='Execution time after bandwidth reduction', ax=ax, color="black", style='-.')
    #plt.title('Crank-Nicolson scheme\'s runtime for solving the Heat Equation in 2D (nsteps=1000)')

    # for (i,j) in zip(N[6:],before[6:]):
    #     plt.annotate(' '+str(j)+'s', (i-175,j), fontsize=9)

    # for (i,j) in zip(N[6:],after[6:]):
    #     plt.annotate(' '+str(j)+'s', (i-100,j+250), fontsize=9)

    plt.xticks(N, fontsize=8)
    plt.yticks(range(0, 12000, 1000))

    plt.xlabel('Matrix dimension')
    plt.ylabel('seconds')
    plt.grid(True)
    plt.savefig('cnexectime.png', bbox_inches='tight')
    #plt.show()

def plot_bandwidth_reduction():
    N = [150, 222, 315, 416, 540, 670, 825, 984, 1170, 1358, 1575]
    before = [26, 55, 63, 91, 100, 153, 165, 209, 222, 299, 315]
    after = [11, 14, 16, 19, 21, 24, 26, 29, 31, 34, 36]

    df = pd.DataFrame({'N': N, 'before': before, 'after': after})

    ax = plt.gca()

    df.plot(kind='line', x='N', y='before', label='Bandwidth', ax=ax, color="black", style='-')
    df.plot(kind='line', x='N', y='after', label='Bandwidth after reduction', ax=ax, color="black", style='-.')
    #plt.title('Bandwidth reduction of the Cuthill-McKee algorithm')

    # for (i,j) in zip(N[1:-1],before[1:-1]):
    #     plt.annotate(' '+str(j), (i+10,j-10), fontsize=9)

    # for (i,j) in zip(N[1:],after[1:]):
    #     plt.annotate(' '+str(j), (i,j+7), fontsize=9)

    # plt.annotate('315', (1555,300), fontsize=9)

    plt.xticks(N, fontsize=8)
    plt.yticks(range(0, 325, 25))

    plt.xlabel('Matrix dimension')
    plt.ylabel('Matrix bandwidth')
    plt.grid(True)
    plt.savefig('cnbwreduc.png', bbox_inches='tight')
    #plt.show()

def main():
    #plot_runtime()
    plot_bandwidth_reduction()

if __name__ == "__main__":
    main()
