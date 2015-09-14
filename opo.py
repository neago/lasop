import numpy as np
import abcd


class OpticalElement:
    def __init__(self):
        pass


class Cavitymode:
    def __init__(self):
        pass


class BowtieOPO:
    def __init__(self, Lc, L1, L2, R1, R2=np.inf, folding_angle=5*np.pi/180, nc=1.8,
                 fixed_first='L2', fixed_second='folding_angle'):
        self._Lc = Lc
        self._L1 = L1
        self._L2 = L2
        self._R1 = R1
        self._R2 = R2
        self._folding_angle = folding_angle
        self._nc = nc
        self.fixed_first = fixed_first
        self.fixed_second = fixed_second

        self.update_geometry()


    @property
    def Lc(self):
        return self._Lc

    @Lc.setter
    def Lc(self, Lc):
        self._Lc = Lc
        self.update_geometry('Lc')

    @property
    def L(self):
        return self._L

    @L.setter
    def L(self, L):
        self._L = L
        self.update_geometry('L')

    @property
    def L2(self):
        return self._L2

    @L2.setter
    def L2(self, L2):
        self._L2 = L2
        self.update_geometry('L2')

    @property
    def L1(self):
        return self._L1

    @L1.setter
    def L1(self, L1):
        self._L1 = L1
        self.update_geometry('L1')

    @property
    def folding_angle(self):
        return self._folding_angle

    @folding_angle.setter
    def folding_angle(self, folding_angle):
        self._folding_angle = folding_angle
        self.update_geometry('folding_angle')

    @property
    def bow_width(self):
        return self._bow_width

    @bow_width.setter
    def bow_width(self, bow_width):
        self._bow_width = bow_width
        self.update_geometry('bow_width')

    @property
    def nc(self):
        return self._nc

    @nc.setter
    def nc(self, nc):
        self._nc = nc
        self.update_geometry('nc')


    def update_geometry(self, update_parameter=None):
        if update_parameter in ('L', 'L2', 'bow_width', 'folding_angle'):
            if update_parameter != self.fixed_first:
                fixed = (update_parameter, self.fixed_first)
            else:
                fixed = (update_parameter, self.fixed_second)
        else:
            fixed = (self.fixed_first, self.fixed_second)


        if 'L' in fixed:
            L_ = self.L - (self.nc - 1) * self.Lc
            if 'L2' in fixed:
                d = (self.L1 + self.L2) / 2
                self._L12 = L_ / 2 - d / 2
                self._folding_angle = np.arccos(d / self._L12)
                self._bow_width = d * np.tan(self.folding_angle)
            elif 'bow_width' in fixed:
                self._L12 = L_ / 4 + self.bow_width**2 / L_
                self._folding_angle = np.arcsin(self.bow_width / self._L12)
                self._L2 = 2 * self.L12 * np.cos(self.folding_angle) - self.L1
            elif 'folding_angle' in fixed:
                self._L12 = L_ / (2 * np.cos(self.folding_angle) + 2)
                self._bow_width = self._L12 * np.sin(self.folding_angle)
                self._L2 = 2 * self._L12 * np.cos(self.folding_angle) - self.L1
        else:
            if 'L2' in fixed and 'bow_width' in fixed:
                d = (self.L1 + self.L2) / 2
                self._L12 = np.sqrt(self.bow_width**2 + d**2)
                self._folding_angle = np.arccos(d / self._L12)
            elif 'L2' in fixed and 'folding_angle' in fixed:
                print('ok')
                d = (self.L1 + self.L2) / 2
                self._L12 = d / np.cos(self.folding_angle)
                self._bow_width = d * np.tan(self.folding_angle)
            elif 'bow_width' in fixed and 'folding_angle' in fixed:
                self._L12 = self.bow_width / np.sin(self.folding_angle)
                self._L2 = 2 * self._L12 * np.cos(self.folding_angle) - self.L1
            self._L = self.L1 + self.L2 + 2 * self._L12 + (self.nc - 1) * self.Lc



        #
        #
        # self.M = ( *
        #           abcd.Minterface(1, n) *
        #           abcd.Mprop((d-l)/2) *
        #           abcd.Mmirror(-R) *
        #           abcd.Mprop(L12) abcd.qpropagate()
        #           abcd.Mmirror(-R2) *
        #           abcd.Mprop(L2) *
        #           abcd.Mmirror(-R2) *
        #           abcd.Mprop(L12) *
        #           abcd.Mmirror(-R1) *
        #           abcd.Mprop((L1 - Lc) / 2) *
        #           abcd.Minterface(n, 1) *
        #           abcd.Mprop(Lc / 2))




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