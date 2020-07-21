"""Plot CHIME repeater fraction."""
import pandas as pd
import numpy as np
import requests
from collections import defaultdict
from copy import deepcopy


def get_chime_rep():
    # Get CHIME repeaters
    json = requests.post("https://catalog.chime-frb.ca/repeaters", data={}).json()

    d = defaultdict(list)

    for frb_name, src_parameters in json.items():

        ra = src_parameters['ra']['value']
        dec = src_parameters['dec']['value']
        gl = src_parameters['gl']['value']
        gb = src_parameters['gb']['value']
        ymw16 = src_parameters['ymw16']['value']
        ne2001 = src_parameters['ne2001']['value']
        publication = src_parameters['publication']['value']
        localized = src_parameters['localized']['value']

        src_dict = {'ra': ra, 'dec': dec, 'gl': gl, 'gb': gb, 'ymw16': ymw16,
                    'ne2001': ne2001, 'publication': publication,
                    'localized': localized, 'name': frb_name}

        dates = [key for key in src_parameters.keys() if key.isdigit()]

        for date in dates:

            date_dict = deepcopy(src_dict)
            for parameter in json[frb_name][date]:
                value = json[frb_name][date][parameter]['value']

                if type(value) is dict:
                    value = np.nan
                date_dict[parameter] = value

                if 'error_low' in json[frb_name][date][parameter]:
                    err_low = json[frb_name][date][parameter]['error_low']
                    err_high = json[frb_name][date][parameter]['error_high']
                else:
                    err_low = np.nan
                    err_high = np.nan

                date_dict[parameter + '_err_low'] = err_low
                date_dict[parameter + '_err_high'] = err_high

            for key, value in date_dict.items():
                d[key].append(value)

    chime_rep = pd.DataFrame(d)
    chime_rep.timestamp = pd.to_datetime(chime_rep.timestamp)
    return chime_rep


def calc_rep_frac():
    df = get_chime_rep()
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
