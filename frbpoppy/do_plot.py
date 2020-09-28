"""Start up Bokeh server."""
import os
import subprocess
import sys

from frbpoppy.misc import pprint
from frbpoppy.paths import paths


def plot(*pops, files=[], tns=False, show=True,
         mute=True, port=5006, print_command=False):
    """
    Plot populations with bokeh. Has to save populations before plotting.

    Args:
        *pops (Population, optional): Add the populations you would like to
            see plotted
        files (list, optional): List of population files to plot.
        tns (bool, optional): Whether to plot tns parameters. Defaults to
            True
        show (bool, optional): Whether to display the plot or not. Mainly used
            for debugging purposes. Defaults to True.
        mute (bool): Show output from Bokeh or not
        port (int): The port on which to launch Bokeh
        print_command (bool): Whether to show the command do_plot is running

    """
    if len(pops) > 0:

        # Save populations
        for pop in pops:

            if type(pop) == str:
                name = pop
            else:
                # Check whether empty population
                if pop.n_sources() < 1:
                    pprint(f'Skipping {pop.name} population as no sources')
                    continue
                pop.name = pop.name.lower()
                if '_for_plotting' not in pop.name:
                    pop.name += '_for_plotting'
                name = pop.name
                pop.save()

            # Save location
            file_name = name + '.p'
            out = os.path.join(paths.populations(), file_name)
            files.append(out)

    # Command for starting up server
    command = 'nice -n 19'.split(' ')

    if show:
        command.extend('bokeh serve --show'.split(' '))
    else:
        command.append('python3')

    # Command for changing port
    if port != 5006:
        command.append(f'--port={port}')

    # Add script path
    script = 'plot.py'
    out = os.path.join(os.path.dirname(__file__), script)
    command.append(out)

    # For the arguments
    command.append('--args')

    # Add tns
    command.append('-tns')
    if tns is False:
        command.append('False')
    if tns is True:
        command.append('True')
    elif type(tns) == str and len(tns) > 0:
        command.append(f'{tns}')

    # Add in populations
    for f in files:
        command.append(f'"{f}"')

    # Check whether plotting is needed
    if len(files) == 0 and tns is False:
        pprint('No populations to plot')
        return

    # Let people know what's happening
    pprint('Plotting populations')

    if print_command:
        pprint(' '.join(command))

    pprint('Press Ctrl+C to quit')

    # Add method to gracefully quit plotting
    try:
        with open(os.devnull, 'w') as f:
            if mute:
                out = f
            else:
                out = None
            subprocess.run(command, stderr=out, stdout=out)
    except KeyboardInterrupt:
        print(' ')
        sys.exit()
