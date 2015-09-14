import numpy as np
import abcd


class OpticalElement:
    def __init__(self):
        pass


class Cavitymode:
    def __init__(self):
        pass


class BowtieOPO:
    def __init__(self, Lc, L1, L2, R1, R2=inf, folding_angle=5, nc=1.8, fixed_first='L', fixed_second='angle'):
        self.Lc = Lc
        self.L1 = L1
        self.L2 = L2
        self.R1 = R1
        self.R2 = R2
        self.folding_angle = folding_angle
        self.nc = nc
        self.fixed_first = fixed_first
        self.fixed_second = fixed_second

        self.update_geometry()


    @property
    def Lc(self):
        return self._Lc

    @Lc.setter
    def Lc(self, Lc):
        self._Lc = Lc
        self.update_geometry()


    @property
    def width_bowtie(self):
        return self._width_bowtie

    @width_bowtie.setter
    def width_bowtie(self, w):
        self._width_bowtie = w

        if self.fixed_folding_angle:
            self._L2 = 2 * w / np.tan(self.folding_angle * np.pi / 180) - self.L1
        else:
            self._folding_angle = 180 / np.pi * np.arctan(2 * w / (self.L1 + self.L2))

    @property
    def L2(self):
        return self._L2




    @property
    def folding_angle(self):
        return self._folding_angle

    def update_geometry(self, update_parameter):
        if update_parameter in ('L', 'L2', 'bow_width', 'folding_angle'):
            fixed = (update_parameter, self.fixed_first)
        else:
            fixed = (self.fixed_first, self.fixed_second)

        if 'L' in fixed:
            L_ = self.L - (self.nc - 1) * self.Lc

        if 'L2' in fixed:
            d = (self.L1 + self.L2) / 2

        if 'L' in fixed:
            if 'L2' in fixed:
                self.L12 = L_ / 2 - d / 2
                self.folding_angle = np.arccos(d / self.L12)
                self.bow_width = d * np.tan(self.folding_angle)
            elif 'bow_width' in fixed:
                self.L12 = L_ / 4 + self.bow_width**2 / L_
                self.folding_angle = np.arcsin(self.bow_width / self.L12)
                self.L2 = 2 * self.L12 * np.cos(self.folding_angle) - self.L1
            elif 'folding_angle' in fixed:
                self.L12 = L_ / (2 * np.cos(self.folding_angle) + 2)
                self.bow_width = self.L12 * np.sin(self.folding_angle)
                self.L2 = 2 * self.L12 * np.cos(self.folding_angle) - self.L1
        else:
            if 'L2' in fixed and 'bow_width' in fixed:
                self.L12 = np.sqrt(self.bow_width**2 + d**2)
                self.folding_angle = np.arccos(d / self.L12)
            elif 'L2' in fixed and 'folding_angle' in fixed:
                self.L12 = d / np.cos(self.folding_angle)
                self.bow_width = d * np.tan(self.folding_angle)
            elif 'bow_width' in fixed and 'folding_angle' in fixed:
                self.L12 = self.bow_width / np.sin(self.folding_angle)
                self.L2 = 2 * self.L12 * np.cos(self.folding_angle) - self.L1
            self.L = 2 * (d + self.L12) + (self.nc - 1) * self.Lc



        if self.fixed_folding_angle:
            self.width_bowtie =

        if self.length_priority = 'total' and self.folding_priority = 'angle':



        d = (self.L1 + self.L2) / 2
        self._L12 = np.sqrt(self.width_bowtie**2 + d**2)
        self.L =

    def _fix_L_L2(self):


    def set_geometry(self, Lc, L1, L2, w_bowtie, R1, R2=inf):
        L12 = sqrt(w_bowtie**2 + ((L1 + L2) / 2)**2)
        self.L = L1 + L2 + 2 * L12 + (self.nc - 1) * Lc
        self.folding_angle = arcsin(w_bowtie / L12)

        self.M = ( *
                  abcd.Minterface(1, n) *
                  abcd.Mprop((d-l)/2) *
                  abcd.Mmirror(-R) *
                  abcd.Mprop(L12) abcd.qpropagate()
                  abcd.Mmirror(-R2) *
                  abcd.Mprop(L2) *
                  abcd.Mmirror(-R2) *
                  abcd.Mprop(L12) *
                  abcd.Mmirror(-R1) *
                  abcd.Mprop((L1 - Lc) / 2) *
                  abcd.Minterface(n, 1) *
                  abcd.Mprop(Lc / 2))




    def get_mode_at(self, mirror, z):
        """Mode at given position after mirror.

        Parameters
        ----------
        mirror : int
            Reference mirror

        z : float
            Distance from mirror. Direction of beam is towards
            the mirror if z is negative.

        Returns
        -------
        mode : Mode
            Mode at given location.
        """