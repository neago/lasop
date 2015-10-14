import numpy as np
import abcd
import matplotlib.pyplot as plt


class OpticalElement:
    def __init__(self):
        pass


class Cavitymode:
    def __init__(self):
        pass


class BowtieOPO:
    """
    Implements bowtie geometry.
    """
    def __init__(self, Lc, L1, L, R1, R2=np.inf, folding_angle=5,
                 nc=1.82, coupling=0.1, loss=.002,
                 fixed_first='L', fixed_second='folding_angle'):
        self._Lc = Lc
        self._L1 = L1
        self._L = L
        self._R1 = R1
        self._R2 = R2
        self._folding_angle = folding_angle * np.pi/180
        self._nc = nc
        self.coupling = coupling
        self.loss = loss
        self.fixed_first = fixed_first
        self.fixed_second = fixed_second

        self.elements, self.M, self.q0 = {}, {}, {}
        self.update_geometry()

    def __str__(self):
        s = 'Mode waist 1: {:} / {:}\n'.format(self.mode_waist('h'), self.mode_waist('v'))
        s += 'Mode waist 2: {:} / {:}\n'.format(self.mode_waist('h', 2), self.mode_waist('v', 2))
        s += 'Eccentricity / coupling to circular mode: {:} / {:}\n'.format(self.eccentricity(2), self.match_to_circular())
        s += 'FSR / Finesse / FWHM: {:} / {:} / {:}'.format(self.FSR(), self.finesse(), self.bandwidth())
        return s


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
    def L1(self):
        return self._L1

    @L1.setter
    def L1(self, L1):
        self._L1 = L1
        self.update_geometry('L1')

    @property
    def L2(self):
        return self._L2

    @L2.setter
    def L2(self, L2):
        self._L2 = L2
        self.update_geometry('L2')

    @property
    def R1(self):
        return self._R1

    @R1.setter
    def R1(self, R1):
        self._R1 = R1
        self.update_geometry('R1')

    @property
    def R2(self):
        return self._R2

    @R2.setter
    def R2(self, R2):
        self._R2 = R2
        self.update_geometry('R2')

    @property
    def folding_angle(self):
        return self._folding_angle * 180 / np.pi

    @folding_angle.setter
    def folding_angle(self, folding_angle):
        self._folding_angle = folding_angle * np.pi / 180
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
                self._L12 = L_ / 2 - d
                self._folding_angle = np.arccos(d / self._L12)
                self._bow_width = d * np.tan(self._folding_angle)
            elif 'bow_width' in fixed:
                self._L12 = L_ / 4 + self.bow_width**2 / L_
                self._folding_angle = np.arcsin(self.bow_width / self._L12)
                self._L2 = 2 * self._L12 * np.cos(self._folding_angle) - self.L1
            elif 'folding_angle' in fixed:
                self._L12 = L_ / (2 * np.cos(self._folding_angle) + 2)
                self._bow_width = self._L12 * np.sin(self._folding_angle)
                self._L2 = 2 * self._L12 * np.cos(self._folding_angle) - self.L1
        else:
            if 'L2' in fixed and 'bow_width' in fixed:
                d = (self.L1 + self.L2) / 2
                self._L12 = np.sqrt(self.bow_width**2 + d**2)
                self._folding_angle = np.arccos(d / self._L12)
            elif 'L2' in fixed and 'folding_angle' in fixed:
                d = (self.L1 + self.L2) / 2
                self._L12 = d / np.cos(self._folding_angle)
                self._bow_width = d * np.tan(self._folding_angle)
            elif 'bow_width' in fixed and 'folding_angle' in fixed:
                self._L12 = self.bow_width / np.sin(self._folding_angle)
                self._L2 = 2 * self._L12 * np.cos(self._folding_angle) - self.L1
            self._L = self.L1 + self.L2 + 2 * self._L12 + (self.nc - 1) * self.Lc

        for dir in ['h', 'v']:
            for d, v in zip([self.elements, self.M, self.q0], self.get_abcd(dir)):
                d[dir] = v

    def get_abcd(self, dir='h'):
        if dir == 'h':
            R1 = self.R1 * np.cos(self._folding_angle)
            R2 = self.R2 * np.cos(self._folding_angle)
        else:
            R1 = self.R1 / np.cos(self._folding_angle)
            R2 = self.R2 / np.cos(self._folding_angle)

        elements_z = [self.Lc / 2,
                      self.L1 / 2,
                      self.L1 / 2 + self._L12,
                      self.L1 / 2 + self._L12 + self.L2,
                      self.L1 / 2 + 2 * self._L12 + self.L2,
                      self.L1 + 2 * self._L12 + self.L2 - self.Lc / 2]
        elements_M = [abcd.Minterface(self.nc, 1),
                      abcd.Mmirror(-R1),
                      abcd.Mmirror(-R2),
                      abcd.Mmirror(-R2),
                      abcd.Mmirror(-R1),
                      abcd.Minterface(1, self.nc)]
        elements = [(ez, eM) for (ez, eM) in zip(elements_z, elements_M)]

        M = np.asmatrix(np.identity(2))
        lastpos = 0
        for e in elements:
            M = e[1] * abcd.Mprop(e[0] - lastpos) * M
            lastpos = e[0]
        M = abcd.Mprop(self.Lc / 2) * M

        q0 = find_stable_mode(M, True)

        return elements, M, q0

    def mode_waist(self, dir='h', waist=1):
        if self.q0[dir]:
            if waist == 1:
                return abcd.q2w(self.q0[dir], self.nc)
            else:
                q = abcd.qpropagate(0, self.q0[dir], self.elements[dir],
                                    self.L1 / 2 + self._L12 + self.L2 / 2)
                return abcd.q2w(q)
        else:
            return np.nan

    def mode_width_at(self, z, dir='h'):
        """
        Get mode radius at arbitrary location.
        """
        if self.q0[dir] == None:
            return np.nan

        q = abcd.qpropagate(0, self.q0[dir], self.elements[dir], z)
        if -self.Lc / 2 <= z < self.Lc / 2:
            return abcd.q2w(q, self.nc)
        else:
            return abcd.q2w(q)

    def eccentricity(self, waist=1):
        a = self.mode_waist('h', waist)
        b = self.mode_waist('v', waist)
        if b > a:
            a, b = b, a
        return np.sqrt(1 - b**2 / a**2)

    def match_to_circular(self):
        """
        Coupling efficiency to a circular mode beam.
        """
        w0h = self.mode_waist('h', 2)
        w0v = self.mode_waist('v', 2)
        w0avg = (w0h + w0v) / 2
        coupling = (4 * w0avg**2 * w0h * w0v /
                    (w0h**2 + w0avg**2) / (w0v**2 + w0avg**2))
        return coupling

    def draw_geometry(self, mirror_diam=8, crystal_width=4):
        fig, ax = plt.subplots(1, 1, figsize=(15,5))
        cv = np.cos(self._folding_angle/2)
        sv = np.sin(self._folding_angle/2)
        cw = crystal_width
        mt, md = 6, mirror_diam  # mirror thickness, diameter
        m1x, m1y = (-self.L1/2, 0)
        m2x, m2y = (self.L1/2, 0)
        m3x, m3y = (-self.L2/2, -self.bow_width)
        m4x, m4y = (self.L2/2, -self.bow_width)
        elements = [plt.Rectangle((-self.Lc / 2, -2), self.Lc, cw),
                    plt.Rectangle((m1x - mt*cv - md/2*sv, m1y + mt*sv - md/2*cv), mt, md, -self.folding_angle/2),
                    plt.Rectangle((m2x + md/2*sv, m2y - md/2*cv), mt, md, self.folding_angle/2),
                    plt.Rectangle((m3x - mt*cv + md/2*sv, m3y - mt*sv - md/2*cv), mt, md, self.folding_angle/2),
                    plt.Rectangle((m4x - md/2*sv, m4y - md/2*cv), mt, md, -self.folding_angle/2)]
                    #plt.Polygon([(m1x, m1y), (m2x, m2y), (m3x, m3y), (m4x, m4y)], edgecolor='k', fill=None, lw=1, alpha=.5)]
        for e in elements:
            ax.add_patch(e)

        x1 = np.linspace(m1x, m2x, 200)
        x2 = np.linspace(m3x, m4x + 20, 200)
        x3 = np.linspace(0, self._L12, 200)
        w01 = self.mode_waist('h', 1)
        w02 = self.mode_waist('h', 2)
        w1 = w01 * np.sqrt(1 + (abcd.lam * x1 / (np.pi * w01**2))**2)
        w2 = w02 * np.sqrt(1 + (abcd.lam * x2 / (np.pi * w02**2))**2)
        w3 = w02 * np.sqrt(1 + (abcd.lam * (x3 + self.L2/2) / (np.pi * w02**2))**2)

        rotmat = lambda v: np.array([[np.cos(v), -np.sin(v)],
                                     [np.sin(v), np.cos(v)]])
        upper1 = np.array([[m3x], [m3y]]) + np.dot(rotmat(self._folding_angle), np.c_[x3, 2*w3].T)
        lower1 = np.array([[m3x], [m3y]]) + np.dot(rotmat(self._folding_angle), np.c_[x3, -2*w3].T)
        upper2 = np.array([[m4x], [m4y]]) + np.dot(rotmat(np.pi - self._folding_angle), np.c_[x3, 2*w3].T)
        lower2 = np.array([[m4x], [m4y]]) + np.dot(rotmat(np.pi - self._folding_angle), np.c_[x3, -2*w3].T)

        ax.add_patch(plt.Polygon(np.c_[upper1, np.fliplr(lower1)].T, lw=0, facecolor='r', alpha=.5))
        ax.add_patch(plt.Polygon(np.c_[upper2, np.fliplr(lower2)].T, lw=0, facecolor='r', alpha=.5))


        ax.fill_between(x1, m1y + 2*w1, m1y - 2*w1, lw=0, facecolor='r', alpha=.5)
        ax.fill_between(x1, m1y + w1, m1y - w1, lw=0, facecolor='r', alpha=.5)
        ax.fill_between(x2, m3y + 2*w2, m3y - 2*w2, lw=0, facecolor='r', alpha=.5)
        ax.fill_between(x2, m3y + w2, m3y - w2, lw=0, facecolor='r', alpha=.5)

        plt.axis('scaled')

        #return fig, ax


    def FSR(self):
        return 2.9979e8 / (.001 * self.L)

    def finesse(self):
        r = np.sqrt((1 - self.coupling) * (1 - self.loss))
        return np.pi / 2 / np.arcsin((1 - r) / (2 * np.sqrt(r)))

    def bandwidth(self):
        """
        FWHM bandwidth.
        """
        return self.FSR() / self.finesse()

    def escape_efficiency(self):
        return self.coupling / (self.coupling + self.loss)


        #self.M_horz =


            #
            #
            # self.M = ( *
            #           abcd.Minterface(1, n) *
            #           abcd.Mprop((d-l)/2) *
            #           abcd.Mmirror(-R) *
            #           abcd.Mprop(L12) abcd.qpropagate()
            #           abcd.Mmirror(-R2) *
            #           abcd.Mprop(L2) *
            #           abcd.Mmrror(-R2) *
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

def find_stable_mode(M, quiet=False):
    """
    Find the stable mode for a cavity with system matrix M.

    Parameters
    ----------
    M : matrix (2x2)
        ABCD matrix for a cavity return trip

    Returns
    -------
    q : complex
        q-parameter of the stable cavity mode.
    """
    A, B, C, D = M[0,0], M[0,1], M[1,0], M[1,1]
    det = np.complex((A-D)**2 + 4*B*C)

    solutions = ((A-D + np.sqrt(det)) / (2*C),
                 (A-D - np.sqrt(det)) / (2*C))

    if solutions[0].imag > 0:
        return solutions[0]
    elif solutions[1].imag > 0:
        return solutions[1]
    else:
        if not quiet:
            print('No solutions')
        #print('No solutions, trying an eps difference...')
        #eps = np.finfo(np.float).eps
        #eps = .1
        #return find_stable_mode(M + np.diag([eps, 0]))
        return None