#!/usr/bin python

import matplotlib.pyplot as plt
import pandas as pd

def plot_time_analysis():
    N = [150, 222, 315, 416, 540, 670]
    before = [19, 52, 127, 276, 582, 1093]
    after = [18, 45, 86, 160, 263, 416]

    df = pd.DataFrame({'N': N, 'before': before, 'after': after})

    ax = plt.gca()

    df.plot(kind='line', x='N', y='before', label='Before applying Cuthill-McKee algorithm', ax=ax)
    df.plot(kind='line', x='N', y='after', label='After applying Cuthill-McKee algorithm', ax=ax)
    plt.title('Time complexity analysis of the Crank-Nicolson scheme for solving 2D Heat Equation (nsteps=1000)')

    plt.xticks(N)
    plt.yticks(before+after)

    plt.xlabel('N')
    plt.ylabel('seconds')
    plt.show()

def plot_bandwidth_reduction():
    N = [150, 222, 315, 416, 540, 670]
    before = [26, 55, 63, 91, 100, 153]
    after = [11, 14, 16, 19, 21, 24]

    df = pd.DataFrame({'N': N, 'before': before, 'after': after})

    ax = plt.gca()

    df.plot(kind='line', x='N', y='before', label='Before applying Cuthill-McKee algorithm', ax=ax)
    df.plot(kind='line', x='N', y='after', label='After applying Cuthill-McKee algorithm', ax=ax)
    plt.title('Bandwidth of Matrices')

    plt.xticks(N)
    plt.yticks(before+after)

    plt.xlabel('N')
    plt.ylabel('Bandwidth')
    plt.show()

def main():
    #plot_time_analysis()
    plot_bandwidth_reduction()

if __name__ == "__main__":
    main()
