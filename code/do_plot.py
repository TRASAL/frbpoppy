"""Start up Bokeh server"""
import os
import subprocess
import sys

from log import pprint

def plot(*pops, files=[], frbcat=True, show=True, mute=True):
    """
    Plot populations with bokeh. Has to save populations before plotting.

    Args:
        *pops (Population, optional): Add the populations you would like to
            see plotted
        files (list, optional): List of population files to plot (currently
            only works with csv files - file an issue if you would like more
            options)
        frbcat (bool, optional): Whether to plot frbcat parameters. Defaults to
            True
        show (bool, optional): Whether to display the plot or not. Mainly used
            for debugging purposes. Defaults to True
    """
    folder = os.path.dirname(__file__)

    if len(pops) > 0:

        # Save populations
        for pop in pops:
            pop.save()

            # Save location
            loc = '../data/results/population_' + pop.name.lower() + '.csv'
            out = os.path.join(folder, loc)
            files.append(out)

    # Command for starting up server
    if show:
        command = 'bokeh serve --show'
    else:
        command = 'python3'

    # Add script path
    script = 'plot.py'
    out = os.path.join(folder, script)
    command += ' ' + out

    # For the arguments
    command += ' --args'

    # Add frbcat
    if not frbcat:
        command += ' -nofrbcat'

    # Add in populations
    for f in files:
        command += ' ' + f

    # Let people know what's happening
    pprint('Plotting populations')
    pprint('Press Ctrl+C to quit')

    # Add method to gracefully quit plotting
    try:
        with open(os.devnull, 'w') as f:
            if mute:
                out = f
            else:
                out = None
            subprocess.run(command.split(' '), stderr=out, stdout=out)
    except KeyboardInterrupt:
        print(' ')
        sys.exit()
