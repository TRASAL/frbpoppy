"""Calculate log N log S in various fluence points."""
import math
import matplotlib.pyplot as plt
from frbpoppy import unpickle

for s in ['perfect-z001', 'perfect-z01', 'perfect-z4']:
    # User input
    pop = unpickle(s)
    parms = pop.get('fluence')

    def calc_alpha(xs):
        """Calculate slope of log N/log S."""
        n = len(xs)
        f_0 = min(xs)
        alpha = -1/((1/n)*sum([math.log(f/f_0) for f in xs]))
        alpha *= (n-1)/n  # Removing bias in alpha
        alpha_err = n*alpha/((n-1)*(n-2)**0.5)
        norm = n / (f_0**alpha)  # Normalisation at lowest parameter
        return alpha, alpha_err, norm

    min_p = round(math.log10(min(parms)))
    max_p = round(math.log10(max(parms)))
    # Points at which to determine alpha (all still in log space)
    parm_points = [i for i in range(min_p+1, max_p)]
    # Range in which to determine alpha
    parm_range = [(e-0.4, e+0.4) for e in parm_points]

    alphas = []
    alpha_errs = []
    norms = []

    for i in range(len(parm_points)):
        point = 10**parm_points[i]
        _range = (10**parm_range[i][0], 10**parm_range[i][1])
        p_range = [p for p in parms if _range[0] <= p <= _range[1]]
        alpha, alpha_err, norm = calc_alpha(p_range)

        alphas.append(alpha)
        alpha_errs.append(alpha_err)
        norms.append(norm)

    # Plot
    xs = [10**p for p in parm_points]
    ys = alphas
    yerrs = alpha_errs

    plt.errorbar(xs, ys, yerr=yerrs, fmt='x', label=s.split('-')[-1])

plt.xscale('log')
plt.xlabel(f'Fluence (Jy ms)')
plt.ylabel(r'$\alpha$')
plt.legend()
plt.tight_layout()
plt.savefig('./plots/alpha_over_fluence.pdf')
