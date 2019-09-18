# -*- coding: future_fstrings -*-
"""Grab survey parameters, and output to a latex table."""
import pandas as pd
from frbpoppy import paths

df_path = paths.surveys() + 'surveys.csv'
df = pd.read_csv(df_path, dtype=object)

# Only included surveys relevant to paper
surveys = ['apertif',
           'askap-fly',
           'askap-incoh',
           'htru',
           'palfa',
           'parkes',
           'guppi',
           'perfect']
df = df[df.survey.isin(surveys)]

# Merge columns
coords = ['RA', 'DEC', 'Galactic longitude', 'Galactic latitude']
for c in coords:
    mi = df[f'minimum {c} (deg)']
    ma = df[f'maximum {c} (deg)']
    df[f'{c} (deg)'] = mi + ' -- ' + ma

# Drop unneccessary columns
mis = [c for c in df.columns if 'minimum' in c]
mas = [c for c in df.columns if 'maximum' in c]
mas = [e for e in mas if e != 'maximum pulse width (ms)']
delete = mis + mas + ['fractional uptime [0-1]', 'reference']
df.drop(columns=delete, inplace=True)

# Conversion table
convert = {'survey': r'Survey',
           'survey degradation factor': r'$\beta$',
           'antenna gain (K/Jy)': r'G (K/Jy)',
           'integration time (s)': r'$t_{\textrm{int}}$ (s)',
           'sampling time (ms)': r'$t_{\textrm{samp}}$ (ms)',
           'receiver temperature (K)': r'$T_{\textrm{rec}}$ (K)',
           'centre frequency (MHz)': r'$\nu_{\textrm{c}}$ (MHz)',
           'bandwidth (MHz)': r'BW (MHz)',
           'channel bandwidth (MHz)': r'BW$_{\textrm{channel}}$ (MHz)',
           'number of polarizations': r'$n_{\textrm{pol}}$',
           'beam size (deg^2)': r'FoV (\degr$^2$)',
           'signal-to-noise ratio [0-1]': r'S/N',
           'maximum pulse width (ms)': r'$w_{\textrm{eff, max}}$',
           'RA (deg)': r'$\alpha$ (\degr)',
           'DEC (deg)': r'$\delta$ (\degr)',
           'Galactic longitude (deg)': r'$l$ (\degr)',
           'Galactic latitude (deg)': r'$b$ (\degr)',
           'reference': r'Ref'}

df.rename(columns=convert, inplace=True)

# Make a vertical table instead
df.set_index('Survey', inplace=True)
df = df.T

def convert_to_latex(df):
    n_cols = df.shape[1]

    s = ['\\begin{tabular}{' + 'c'*(n_cols+2) + '}',
         '\hline',
         '\hline']

    c = ' & '.join(df.columns) + "\\\\"
    s.extend([f'Parameter & Units & {c}', '\hline'])

    for index, row in df.iterrows():

        # Split units
        entry = index.strip(')').split('(')
        if len(entry) == 2:
            par, unit = entry
        else:
            par = entry[0]
            unit = ''

        line = par + f' & {unit}'

        for e in row:
            v = f' & {e}'

            if isinstance(e, str):
                if e.startswith('http'):
                    v = ' & \\url{' + e + '}'

            line += v
        line += '\\\\'
        s.append(line)

    e = ['\hline',
         '\end{tabular}']

    s.extend(e)

    return '\n'.join(s)


output = paths.code + '/../publication/surv_par.tex'
with open(output, 'w') as tf:
    tf.write(convert_to_latex(df))
