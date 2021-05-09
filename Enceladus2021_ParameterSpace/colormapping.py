""" for making colormaps, specifically ones which fade to transparency """

import matplotlib.colors as clr

class cmapper:

    green_to_a = {'red':[(0,0,0), (1,0,0)],
      'green':[(0,0.7,0.7),(1,0.7,0.7)],
      'blue':[(0,0,0), (1,0,0)],
      'alpha':[(0,1,1), (1,0,0)]}

    blue_to_a = {'red':[(0,0,0), (1,0,0)],
      'green':[(0,0,0),(1,0,0)],
      'blue':[(0,1,1), (1,1,1)],
      'alpha':[(0,1,1), (1,0,0)]}

    red_to_a = {'red':[(0,1,1), (1,1,1)],
      'green':[(0,0,0),(1,0,0)],
      'blue':[(0,0,0), (1,0,0)],
      'alpha':[(0,1,1), (1,0,0)]}

    k_to_a = {'red':[(0,0,0), (1,0,0)],
      'green':[(0,0,0),(1,0,0)],
      'blue':[(0,0,0), (1,0,0)],
      'alpha':[(0,1,1), (1,0,0)]}


    def transparent_cmap(rgba, name):
        return clr.LinearSegmentedColormap(name, rgba, N=100).reversed()

    def g2a():
        return cmapper.transparent_cmap(cmapper.green_to_a, 'green_to_a')

    def r2a():
        return cmapper.transparent_cmap(cmapper.red_to_a, 'red_to_a')

    def b2a():
        return cmapper.transparent_cmap(cmapper.blue_to_a, 'blue_to_a')

    def k2a():
        return cmapper.transparent_cmap(cmapper.k_to_a, 'black_to_a')
