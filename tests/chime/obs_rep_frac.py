"""Calculate the actual CHIME repeater fraction."""
import pandas as pd

from frbcat import ChimeRepeaters


def calc_rep_frac():
    df = ChimeRepeaters.df()
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
