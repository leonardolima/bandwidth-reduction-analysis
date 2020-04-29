#!/usr/bin python

import matplotlib.pyplot as plt
import pandas as pd

def remove_too_close(l):
    aux_l = [l[0]]

    for i in range(1, len(l)):
        if l[i]-l[i-1] > 2:
            aux_l.append(l[i])

    return aux_l


def plot_bandwidth_comparison():

    files = ['DateTimeNow_yxmc.txt', 'CountRecords_yxmc.txt', 'WeightedAvg_yxmc.txt',
             'parallel.txt', 'FooterMacro_yxmc.txt', 'Legend_Splitter_yxmc.txt',
             'Base64_Encoder_yxmc.txt', 'SelectRecords_yxmc.txt', 'Date_Filter_yxmc.txt',
             'RandomRecords_yxmc.txt', 'PearsonCorrCoeff_yxmc.txt', 'Frequency_yxmc.txt',
             'HeaderMacro_yxmc.txt', 'Cleanse_yxmc.txt', 'SpearmanCorrCoeff_yxmc.txt',
             'PieWedgeTradeArea_yxmc.txt', 'Imputation_yxmc.txt', 'Imputation_v2_yxmc.txt',
             'MultiFieldBinning_yxmc.txt', 'HeatMap_yxmc.txt', 'Google_Analytics_yxmc.txt',
             'Google_Analytics_v5_yxmc.txt', 'Google_Analytics_v6_yxmc.txt']

    N = [7, 9, 19, 20, 21, 23, 23, 24, 25, 29, 34, 34, 34,
         35, 48, 58, 61, 81, 81, 95, 108, 123, 123]

    before = [4, 3, 15, 16, 13, 15, 15, 20, 24, 23, 22, 27,
              25, 24, 22, 44, 55, 71, 78, 67, 83, 74, 74]
    after =  [3, 3, 8, 11, 13, 12, 19, 9, 12, 9, 10, 15,
              26, 27, 7, 22, 44, 71, 74, 60, 70, 98, 98]

    df = pd.DataFrame({'files': files, 'N': N, 'before': before, 'after': after})

    ax = plt.gca()

    df.plot(kind='line', x='N', y='before', label='Before', ax=ax, style='.-')
    df.plot(kind='line', x='N', y='after', label='After', ax=ax, style='.-')
    # df.plot(kind='scatter', x='N', y='before', label='Before', ax=ax)
    # df.plot(kind='scatter', x='N', y='after', label='After', ax=ax)

    plt.title('Before and after applying the adapted Cuthill-McKee algorithm')

    # yticks = before+after
    # yticks.sort()

    # plt.xticks(remove_too_close(N))
    # plt.yticks(remove_too_close(yticks))

    plt.grid(True, linestyle='--')

    plt.xlabel('Dimension')
    plt.ylabel('Bandwidth')
    # plt.savefig('adapted_results.png')
    plt.show()

def main():
    plot_bandwidth_comparison()

if __name__ == "__main__":
    main()
