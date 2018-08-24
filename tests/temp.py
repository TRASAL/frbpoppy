from bokeh.models.widgets import Panel, Tabs
from bokeh.io import output_file, show
from bokeh.plotting import figure

output_file("slider.html")

tabs = []
for i in range(20):
    p1 = figure(plot_width=300, plot_height=300)
    p1.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5)
    tabs.append(Panel(child=p1, title="circle"))

show(Tabs(tabs=tabs))
