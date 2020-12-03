"""Calculate the actual CHIME repeater fraction."""
import pandas as pd

from frbcat import ChimeRepeaters


def calc_rep_frac():
    df = ChimeRepeaters().df
    df.timestamp = pd.to_datetime(df.timestamp)
    df.sort_values('timestamp', inplace=True)

    starting_time = df['timestamp'].iloc[0]
    srcs_seen_once = []
    srcs_seen_twice = []
    time_rep_seen = []
    for index, row in df.iterrows():
        src = row['name']
        dt = (row['timestamp'] - starting_time).total_seconds() / 86400

        if src not in srcs_seen_once:
            srcs_seen_once.append(src)
            continue
        elif src not in srcs_seen_twice:
            # Timedelta at which repeater detected
            time_rep_seen.append(dt)
            srcs_seen_twice.append(src)

    return time_rep_seen


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from tests.convenience import plot_aa_style, rel_path
    import numpy as np
    from tqdm import tqdm

    # Set up plot style
    plot_aa_style(cols=1)
    f, ax1 = plt.subplots(1, 1)

    # See how the real fraction changes over time
    chime_fracs = []
    dts = calc_rep_frac()
    days = [d for d in range(301)]
    for day in tqdm(days, desc='frbcat'):
        n_rep = sum([dt <= day for dt in dts])
        n_one_offs = 2*day
        try:
            frac = n_rep / (n_rep + n_one_offs)
        except ZeroDivisionError:
            frac = np.nan
        chime_fracs.append(frac)

    ax1.plot(days, chime_fracs, label='chime-frb')

    # Further plot details
    ax1.set_xlabel(r'Time (days)')
    ax1.set_ylabel(r'$f_{\textrm{rep}}$')
    # Equal to $N_{\textrm{repeaters}}/N_{\textrm{detections}}$
    ax1.set_xlim(0, max(days))
    ax1.set_yscale('log')

    # Save figure
    plt.tight_layout()
    plt.savefig(rel_path('plots/obs_rep_frac_chime.pdf'))
    plt.clf()
