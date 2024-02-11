
import matplotlib.pyplot as plt
import numpy as np
import math, os, sys
sys.path.append(os.path.dirname(__file__)+'/../')
from EncBmBs_utils import PlotSetup


class pH_T_LineGenerator:

    def __init__(self, DFF):
        """
        Class for drawing lines on a speciation plot exploring trends between
        pH and temperature.

        DFF is a DataFrameFetcher object specified with the necessary
        environmental parameters.
        """

        self.DFF = DFF

        # setup colomaps for plotting
        self.cmap = plt.get_cmap('Set2', len(self.DFF.pHs))
        self.cmaplist = [self.cmap(i) for i in range(self.cmap.N)]
        if len(self.DFF.pHs) == 6:
            self.cmaplist[4] = 'chocolate'
            self.cmaplist[3] = 'yellowgreen'

        self.lineslist = list(PlotSetup.mpl_linestyles.values())[:len(self.DFF.pHs)]



    def plot_spec_lines(self, ax, Zname, salt_lvl='nom',
      other_df=None, return_df=False):

        """
        Plot the parameter Zname and its variation with temperature on the x axis
        and line colors for pH, at the passed salt_lvl which can be a string or
        a float.

        Pass a pandas DataFrame as other_df if you do not want this pH_T_LineGenerator's
        DFF attribute to be used as the source speciation.
        """

        this_df = other_df
        if type(this_df) == type(None):
            this_df = self.DFF.spec_Z(Zname)
        if type(salt_lvl) != type(''):
            this_df = this_df[np.isclose(this_df['salt_lvl'], salt_lvl)]
        else:
            this_df = this_df[this_df['salt_lvl'] == salt_lvl]
        for _i, _pH in enumerate(self.DFF.pHs):
            this_pH_df = this_df[this_df['pH_bo']==_pH]
            ax.plot(this_pH_df['T'], this_pH_df[Zname], c=self.cmaplist[_i], label='O$\degree$C pH = '+str(_pH), ls=self.lineslist[_i])
            ax.set_xlabel('Temperature [K]')
            ax.set_ylabel(Zname)
        if return_df:
            return ax, this_df
        else:
            return ax

    def plot_ratio_lines(self, ax, Ztop, Zbottom, salt_lvl='nom'):
        """
        Plot the ratio between Ztop and Zbottom and its variation with
        temperature on the x axis and line colors for pH, at the passed
        salt_lvl which can be a string or a float.
        """
        this_df = self.DFF.spec()
        this_df = this_df[this_df['salt_lvl'] == salt_lvl]

        this_df['ratio'] = this_df[Ztop] / this_df[Zbottom]

        ax = self.plot_spec_lines(ax, 'ratio', salt_lvl, this_df)
        return ax
