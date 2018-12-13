import pandas as pd
import os

from frbpoppy import paths

path = os.path.join(paths.surveys(), 'surveys.csv')

df = pd.read_csv(path)
df = df.set_index('survey')
survey = df[df.index == 'apertif'].squeeze()

print(df.index)
