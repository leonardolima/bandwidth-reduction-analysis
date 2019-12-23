#!/usr/bin python

import matplotlib.pyplot as plt
import pandas as pd

def plot_time_analysis():
    N = [150, 222, 315, 416]
    before = [170, 533, 1512, 3552]
    after = [163, 524, 1461, 3468]

    df = pd.DataFrame({'N': N, 'before': before, 'after': after})

    ax = plt.gca()

    df.plot(kind='line', x='N', y='before', label='Before applying Cuthill-McKee algorithm', ax=ax)
    df.plot(kind='line', x='N', y='after', label='After applying Cuthill-McKee algorithm', ax=ax)
    plt.title('Time complexity analysis of the Crank-Nicolson scheme for solving 2D Heat Equation')

    plt.xticks(N)
    plt.yticks(before+after[-2:])

    plt.xlabel('N')
    plt.ylabel('seconds')
    plt.show()

def plot_bandwidth_reduction():
    N = [150, 222, 315, 416]
    before = [26, 55, 63, 91]
    after = [11, 14, 16, 19]

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
    #plot_bandwidth_reduction()

if __name__ == "__main__":
    main()
