"""
Formulas for PPKTP - periodically poled potassium titanyl phosphate.
Dispersion, phase matching, nonlinearity, etc.

Wavelengths (lam) are given in µm.

References:

* Konig & Wong, Appl. Phys. Lett. 84, 1644 (2004)
* Kato & Takaoka, Appl. Optics 41, 5040 (2002)
* Emanueli & Arie, Appl. Optics 42, 6661 (2003)
* Fan et al, Appl. Opt. 26, 2390 (1987)
* Fradkin et al, Appl. Phys. Lett. 74, 914 (1999)

"""

def sellmeier(lam, *coefficients):
    """Sellmeier equation for refractive index at 
       wavelength _lam_ (um) for crystal with given coefficients."""
    c0, c1, c2, c3, c4 = coefficients
    return sqrt(c0 + c1/(lam**2 - c2) + c3/(lam**2 - c4))

def deltan_coeff(lam, *coefficients):
    """3rd order polynomial coefficients of the 2nd order
       temperature correction polynomial."""
    a0, a1, a2, a3 = coefficients
    return a0 + a1/lam + a2/lam**2 + a3/lam**3

def nx(lam, T=25, source='Kato2002'):
    """KTP refractive index in x axis"""
    return sellmeier(lam, 3.29100, 0.04140, 0.03978, 9.35522, 31.45571)

def ny(lam, T=25, source='Konig2004'):
    """KTP refractive index in y axis"""
    if source == 'Kato2002':
        ny_20 = sellmeier(lam, 3.45018, 0.04341, 0.04597, 16.98825, 39.43799)
    elif source == 'Fan1987':
        b0, b1, b2, b3 = (2.19229, 0.83547, 0.04970, 0.01621)
        ny_25 = sqrt(b0 + b1 / (1 - b2 / lam**2) - b3 * lam**2)
    elif source == 'Konig2004':
        b0, b1, b2, b3 = (2.09930, 0.922683, 0.0467695, 0.0138408)
        ny_25 = sqrt(b0 + b1 / (1 - b2 / lam**2) - b3 * lam**2)
        
    ny1 = deltan_coeff(lam, 6.2897e-6, 6.3061e-6, -6.0629e-6, 2.6486e-6)
    ny2 = deltan_coeff(lam, -.14445e-8, 2.2244e-8, -3.5770e-8, 1.3470e-8)
    
    if source == 'Kato2002':
        return ny_20 + ny1 * (T - 20) + ny2 * (T - 20)**2
    else:
        return ny_25 + ny1 * (T - 25) + ny2 * (T - 25)**2

def nz(lam, T=25, source='Fradkin1999'):
    """KTP refractive index in z axis"""
    if source == 'Kato2002':
        nz_20 = sellmeier(lam, 4.59423, 0.06206, 0.04763, 110.80672, 86.12171)
    elif source == 'Fradkin1999':
        b0, b1, b2, b3, b4, b5 = (2.12725, 1.18431, .0514852, 
                                  0.66030, 100.00507, 9.68956e-3)
        nz_25 = sqrt(b0 + b1 / (1 - b2 / lam**2) + 
                     b3 / (1 - b4 / lam**2) - b5 * lam**2)
    elif source == 'Fan1987':
        b0, b1, b2, b3 = (2.25411, 1.06543, 0.05486, 0.02140)
        nz_25 = sqrt(b0 + b1 / (1 - b2 / lam**2) - b3 * lam**2)
        
    nz1 = deltan_coeff(lam, 9.9587e-6, 9.9228e-6, -8.9603e-6, 4.1010e-6)
    nz2 = deltan_coeff(lam, -1.1882e-8, 10.459e-8, -9.8136e-8, 3.1481e-8)
    
    if source == 'Kato2002':
        return nz_20 + nz1 * (T - 20) + nz2 * (T - 20)**2
    else:
        return nz_25 + nz1 * (T - 25) + nz2 * (T - 25)**2

    
### Nonlinear conversion efficiency

def pwpm_typeII(lam, T, poling_period, L,
                source_y='Konig2004', source_z='Fradkin1999'):
    """Plane-wave phase matching condition.
    """
    k1 = 2 * pi * ny(lam, T, source_y) / lam
    k2 = 2 * pi * nz(lam, T, source_z) / lam
    k3 = 2 * pi * ny(lam / 2, T, source_y) / (lam / 2)
    dk = k3 - k2 - k1 + 2 * pi / poling_period
    x = dk * L / 2
    return (sin(x) / x)**2

@vectorize
def BKh(L, zR, mu, dk, alpha_fund=0, alpha_harm=1e-5):
    """
    Boyd-Kleinman h function.
    B-K eqs. (2.23-25) with updated real-valued integrand, from 
    [G. Masada](http://www.tamagawa.jp/en/research/quantum/bulletin/pdf/Tamagawa.Vol.3-4.pdf)
    eq. 10).
    
    L:     crystal length [µm]
    zR:    Rayleigh length [µm]
    mu:    focus position in the crystal, range -1 to 1 (0 is center)
    dk:    wavevector mismatch [1/µm]
    alpha_fund: linear absorption at fundamental wavelength [1/µm]
    alpha_harm: linear absorption at second harmonic wavelength [1/µm]
    """
    xi = L/(2*zR)
    alpha = alpha_fund - alpha_harm/2
    kappa = alpha*zR
    sigma = dk*zR
    
    def integrand(t):
        return (cos(sigma*t) + t*sin(sigma*t)) / (1 + t**2)
    
    integral = quad(integrand, -xi*(1-mu), xi*(1+mu))[0]
    
    return 1/(4*xi) * exp(mu*alpha*L) * integral**2

def BKhT_type0(lam, T, poling_period, L, zR, mu, alpha_fund=0, alpha_harm=0, source='Fradkin1999'):
    k1 = 2 * pi * nz(lam, T, source) / lam
    k2 = k1
    k3 = 2 * pi * nz(lam / 2, T, source) / (lam / 2)
    dk = k3 - k2 - k1 - 2 * pi / poling_period
    return BKh(L, zR, mu, dk, alpha_fund, alpha_harm)

def BKhT_typeII(lam, T, poling_period, L, zR, mu, alpha_fund=0, alpha_harm=0, 
                source_y='Konig2004', source_z='Fradkin1999'):
    k1 = 2 * pi * ny(lam, T, source_y) / lam
    k2 = 2 * pi * nz(lam, T, source_z) / lam
    k3 = 2 * pi * ny(lam / 2, T, source_y) / (lam / 2)
    dk = k3 - k2 - k1 + 2 * pi / poling_period  # note the changed sign of the poling wave vector
    return BKh(L, zR, mu, dk, alpha_fund, alpha_harm)

def Enl(lam, deff, L, zR, mu, dk, alpha_fund=0, alpha_harm=0):
    """
    Nonlinear conversion efficiency $E_{nl}$, 
    from Boyd-Kleinman ($P_{SHG} = E_{nl} P_{fund.}^2$).
    Units are µm.
    """
    h = BKh(L, zR, mu, dk, alpha_fund, alpha_harm)
    atten = exp(-(alpha_fund + alpha_harm/2)*L)
    eps0 = 8.854188e-18
    c = 2.99792458e14
    factor = 16 * pi**2 * deff**2 * L / (eps0 * c * lam**3 * nz(lam) * nz(lam/2))
    
    return factor * h * atten

@vectorize
def Enl_matched(lam, deff, L, zR, mu, alpha_fund=0, alpha_harm=0):
    """Enl at optimal temperature.
    Units are µm.
    """
    dk_opt = optimize.fmin(lambda dk: -Enl(lam, deff, L, zR, mu, dk, alpha_fund, alpha_harm), 0,
                          disp=False, full_output=False)
    return Enl(lam, deff, L, zR, mu, dk_opt, alpha_fund, alpha_harm)