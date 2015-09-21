import numpy as np
import abcd


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
    def __init__(self, Lc, L1, L, R1, R2=np.inf, folding_angle=5*np.pi/180, nc=1.8,
                 fixed_first='L', fixed_second='folding_angle'):
        self._Lc = Lc
        self._L1 = L1
        self._L = L
        self._R1 = R1
        self._R2 = R2
        self._folding_angle = folding_angle
        self._nc = nc
        self.fixed_first = fixed_first
        self.fixed_second = fixed_second

        self.elements, self.M, self.q0 = {}, {}, {}
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
                self._L12 = L_ / 2 - d
                self._folding_angle = np.arccos(d / self._L12)
                self._bow_width = d * np.tan(self.folding_angle)
            elif 'bow_width' in fixed:
                self._L12 = L_ / 4 + self.bow_width**2 / L_
                self._folding_angle = np.arcsin(self.bow_width / self._L12)
                self._L2 = 2 * self._L12 * np.cos(self.folding_angle) - self.L1
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
                d = (self.L1 + self.L2) / 2
                self._L12 = d / np.cos(self.folding_angle)
                self._bow_width = d * np.tan(self.folding_angle)
            elif 'bow_width' in fixed and 'folding_angle' in fixed:
                self._L12 = self.bow_width / np.sin(self.folding_angle)
                self._L2 = 2 * self._L12 * np.cos(self.folding_angle) - self.L1
            self._L = self.L1 + self.L2 + 2 * self._L12 + (self.nc - 1) * self.Lc

        for dir in ['h', 'v']:
            for d, v in zip([self.elements, self.M, self.q0], self.get_abcd(dir)):
                d[dir] = v

    def update_abcd_old(self):
        elements_z = [self.Lc / 2,
                      self.L1 / 2,
                      self.L1 / 2 + self._L12,
                      self.L1 / 2 + self._L12 + self.L2,
                      self.L1 / 2 + 2 * self._L12 + self.L2,
                      self.L1 / 2 + 2 * self._L12 + self.L2 + self.Lc / 2]
        elements_M_horz = [abcd.Minterface(self.nc, 1),
                           abcd.Mmirror(-self.R1 * np.cos(self.folding_angle)),
                           abcd.Mmirror(-self.R2 * np.cos(self.folding_angle)),
                           abcd.Mmirror(-self.R2 * np.cos(self.folding_angle)),
                           abcd.Mmirror(-self.R1 * np.cos(self.folding_angle)),
                           abcd.Minterface(1, self.nc)]
        elements_M_vert = [abcd.Minterface(self.nc, 1),
                           abcd.Mmirror(-self.R1 / np.cos(self.folding_angle)),
                           abcd.Mmirror(-self.R2 / np.cos(self.folding_angle)),
                           abcd.Mmirror(-self.R2 / np.cos(self.folding_angle)),
                           abcd.Mmirror(-self.R1 / np.cos(self.folding_angle)),
                           abcd.Minterface(1, self.nc)]
        elements_horz = zip(elements_z, elements_M_horz)
        elements_vert = zip(elements_z, elements_M_vert)

        M_horz = np.asmatrix(np.identity(2))
        M_vert = np.asmatrix(np.identity(2))
        lastpos = 0
        for eh, ev in zip(elements_horz, elements_vert):
            M_horz = eh[1] * abcd.Mprop(eh[0] - lastpos) * M_horz
            M_vert = ev[1] * abcd.Mprop(ev[0] - lastpos) * M_vert
            lastpos = eh[0]
        M_horz = abcd.Mprop(self.Lc / 2) * M_horz
        M_vert = abcd.Mprop(self.Lc / 2) * M_vert

        q0_horz = find_stable_mode(M_horz, True)
        q0_vert = find_stable_mode(M_vert, True)

        return ((elements_horz, M_horz, q0_horz), (elements_vert, M_vert, q0_vert))

    def get_abcd(self, dir='h'):
        if dir == 'h':
            R1 = self.R1 * np.cos(self.folding_angle)
            R2 = self.R2 * np.cos(self.folding_angle)
        else:
            R1 = self.R1 / np.cos(self.folding_angle)
            R2 = self.R2 / np.cos(self.folding_angle)

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
        if waist == 1:
            return abcd.q2w(self.q0[dir], self.nc)
        else:
            q = abcd.qpropagate(0, self.q0[dir], self.elements[dir],
                                self.L1 / 2 + self._L12 + self.L2 / 2)
            return abcd.q2w(q)

    def mode_width_at(self, z, dir='h'):
        """
        Get mode radius at arbitrary location.
        """
        q = abcd.qpropagate(0, self.q0[dir], self.elements[dir], z)
        if -self.Lc / 2 <= z < self.Lc / 2:
            return abcd.q2w(q, self.nc)
        else:
            return abcd.q2w(q)

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