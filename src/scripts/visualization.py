# %%
import numpy as np
import pandas as pd

from lets_plot.mapping import as_discrete
from lets_plot import *
# from plotnine import ggplot, aes, facet_grid, labs, geom_point
# from plotnine.data import mpg

LetsPlot.setup_html()
LetsPlot.set_theme(theme_light())

# %%
palette = scale_color_manual(
    values=["#394449", "#F7C443"]) + scale_fill_manual(values=["#394449", "#F7C443"])
np.random.seed(0)

cov0 = [[1, -.8],
        [-.8, 1]]
cov1 = [[10, .1],
        [.1, .1]]

x0, y0 = np.random.multivariate_normal(mean=[-2, 0], cov=cov0, size=200).T
x1, y1 = np.random.multivariate_normal(mean=[0, 1], cov=cov1, size=200).T

data = dict(
    x=np.concatenate((x0, x1)),
    y=np.concatenate((y0, y1)),
    c=["A"]*200 + ["B"]*200
)

plot_settings = (ggsize(1200, 800) +
                 theme(plot_background=element_rect(fill="#eaeaea"),
                       legend_background=element_rect(fill="#eaeaea")))

p = ggplot(data, aes("x", "y", color="c", fill="c")) + geom_point() + palette
(p + ggmarginal("rb", size=0.4, layer=geom_violin(trim=False, color="black"))
 + ggmarginal("rb", layer=geom_boxplot(aes(group="c"), fill="white", color="white",
                                                                           alpha=.25, outlier_color="red", width=.2)) + plot_settings)
