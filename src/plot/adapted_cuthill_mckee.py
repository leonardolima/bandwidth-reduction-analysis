#!/usr/bin python

import matplotlib.pyplot as plt
import pandas as pd
import itertools

def plot_bandwidth_comparison():

    files = ['example.txt', 'DateTimeNow_yxmc.txt', 'CountRecords_yxmc.txt', 'WeightedAvg_yxmc.txt',
             'parallel.txt', 'FooterMacro_yxmc.txt', 'Legend_Splitter_yxmc.txt',
             'Base64_Encoder_yxmc.txt', 'SelectRecords_yxmc.txt', 'Date_Filter_yxmc.txt',
             'RandomRecords_yxmc.txt', 'PearsonCorrCoeff_yxmc.txt', 'Frequency_yxmc.txt',
             'HeaderMacro_yxmc.txt', 'Cleanse_yxmc.txt', 'SpearmanCorrCoeff_yxmc.txt',
             'PieWedgeTradeArea_yxmc.txt', 'Imputation_yxmc.txt', 'Imputation_v2_yxmc.txt',
             'MultiFieldBinning_yxmc.txt', 'HeatMap_yxmc.txt', 'Google_Analytics_yxmc.txt',
             'Google_Analytics_v5_yxmc.txt', 'Google_Analytics_v6_yxmc.txt']

    N_unfiltered = [6, 7, 9, 19, 20, 21, 23, 23, 24, 25, 29, 34, 34, 34, 35, 48, 58, 61, 81, 81, 95, 108, 123, 123]
    bw_before_unfiltered = [2, 4, 3, 15, 16, 13, 15, 15, 20, 24, 23, 27, 25, 22, 24, 22, 44, 55, 71, 78, 67, 83, 74, 74]
    bw_after_unfiltered = [2, 3, 3, 14, 16, 13, 18, 8, 14, 20, 18, 17, 27, 12, 24, 12, 39, 34, 71, 75, 66, 64, 102, 102]

    N, bw_before, bw_after = N_unfiltered, bw_before_unfiltered, bw_after_unfiltered
    # N, bw_before, bw_after = [], [], []

    # for i in range(0,len(N_unfiltered)):
    #     if i > 0:
    #         if N_unfiltered[i] != N_unfiltered[i-1]:
    #             N.append(N_unfiltered[i])
    #             bw_before.append(bw_before_unfiltered[i])
    #             bw_after.append(bw_after_unfiltered[i])
    #     if bw_before_unfiltered[i] == bw_after_unfiltered[i]:
    #         N.append(N_unfiltered[i])
    #         bw_before.append(bw_before_unfiltered[i])
    #         bw_after.append(bw_after_unfiltered[i])

    df = pd.DataFrame({'N': N, 'bw_before': bw_before, 'bw_after': bw_after})

    ax = plt.gca()

    df.plot(kind='line', x='N', y='bw_before', label='Bandwidth', ax=ax, color="black", style='-')
    df.plot(kind='line', x='N', y='bw_after', label='Bandwidth after applying ACM', ax=ax, color="black", style='-.')
    #plt.title('Crank-Nicolson scheme\'s runtime for solving the Heat Equation in 2D (nsteps=1000)')

    # for (i,j) in zip(N[6:],before[6:]):
    #     plt.annotate(' '+str(j)+'s', (i-175,j), fontsize=9)

    # for (i,j) in zip(N[6:],after[6:]):
    #     plt.annotate(' '+str(j)+'s', (i-100,j+250), fontsize=9)

    plt.xticks(range(0, 140, 10), fontsize=8)
    plt.yticks(range(0, 120, 10), fontsize=8)

    plt.xlabel('Matrix dimension')
    plt.ylabel('Bandwidth')
    plt.grid(True)
    plt.savefig('acm_bw.png', bbox_inches='tight')
    plt.show()

def plot_peak_comparison():
    files = ['example.txt', 'DateTimeNow_yxmc.txt', 'CountRecords_yxmc.txt', 'WeightedAvg_yxmc.txt',
             'parallel.txt', 'FooterMacro_yxmc.txt', 'Legend_Splitter_yxmc.txt',
             'Base64_Encoder_yxmc.txt', 'SelectRecords_yxmc.txt', 'Date_Filter_yxmc.txt',
             'RandomRecords_yxmc.txt', 'PearsonCorrCoeff_yxmc.txt', 'Frequency_yxmc.txt',
             'HeaderMacro_yxmc.txt', 'Cleanse_yxmc.txt', 'SpearmanCorrCoeff_yxmc.txt',
             'PieWedgeTradeArea_yxmc.txt', 'Imputation_yxmc.txt', 'Imputation_v2_yxmc.txt',
             'MultiFieldBinning_yxmc.txt', 'HeatMap_yxmc.txt', 'Google_Analytics_yxmc.txt',
             'Google_Analytics_v5_yxmc.txt', 'Google_Analytics_v6_yxmc.txt']

    N_unfiltered = [6, 7, 9, 19, 20, 21, 23, 23, 24, 25, 29, 34, 34, 34, 35, 48, 58, 61, 81, 81, 95, 108, 123, 123]
    peak_random_unfiltered = [15, 147, 166, 482, 331, 588, 654, 662, 495, 303, 717, 1441, 892, 778,  680, 1565, 1051, 3218, 3259, 2114, 1784, 7789, 9122, 9120]
    peak_acm_unfiltered = [15, 147, 166, 503, 331, 588, 781, 662, 656, 303, 717, 1431, 864, 582, 1178, 1687,1202, 3293, 3403, 2114, 2145, 10260, 10781, 10780]

    N, peak_random, peak_acm = N_unfiltered, peak_random_unfiltered, peak_acm_unfiltered
    # Filtering
    # N, peak_random, peak_acm = [], [], []

    # for i in range(0,len(N_unfiltered)):
    #     if i > 0:
    #         if N_unfiltered[i] != N_unfiltered[i-1]:
    #             N.append(N_unfiltered[i])
    #             peak_random.append(peak_random_unfiltered[i])
    #             peak_acm.append(peak_acm_unfiltered[i])
    #     if peak_random_unfiltered[i] == peak_acm_unfiltered[i]:
    #         N.append(N_unfiltered[i])
    #         peak_random.append(peak_random_unfiltered[i])
    #         peak_acm.append(peak_acm_unfiltered[i])

    df = pd.DataFrame({'N': N, 'peak_random': peak_random, 'peak_acm': peak_acm})

    ax = plt.gca()

    df.plot(kind='line', x='N', y='peak_random', label='PMU', ax=ax, color="black", style='-')
    df.plot(kind='line', x='N', y='peak_acm', label='PMU applying ACM', ax=ax, color="black", style='-.')
    #plt.title('Crank-Nicolson scheme\'s runtime for solving the Heat Equation in 2D (nsteps=1000)')

    # for (i,j) in zip(N[6:],before[6:]):
    #     plt.annotate(' '+str(j)+'s', (i-175,j), fontsize=9)

    # for (i,j) in zip(N[6:],after[6:]):
    #     plt.annotate(' '+str(j)+'s', (i-100,j+250), fontsize=9)

    plt.xticks(range(0,130,10), fontsize=8)
    # plt.xticks(N, fontsize=8)
    plt.yticks(range(0, 13000, 1000), fontsize=8)

    plt.xlabel('Matrix dimension')
    plt.ylabel('Peak memory usage')
    plt.grid(True)
    plt.savefig('acm_peak.png', bbox_inches='tight')
    plt.show()

def main():
    plot_peak_comparison()
    plot_bandwidth_comparison()

if __name__ == "__main__":
    main()
