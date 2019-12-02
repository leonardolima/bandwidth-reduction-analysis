#!/usr/bin python

import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("out.csv", sep=',', names=['N', 'lubksb', 'tridag'])

ax = plt.gca()

df.plot(kind='line', x='N', y='lubksb', ax=ax)
df.plot(kind='line', x='N', y='tridag', ax=ax)

plt.title('Comparison of lubksb and tridag methods')
plt.xlabel('N')
plt.ylabel('milliseconds')
plt.show()
