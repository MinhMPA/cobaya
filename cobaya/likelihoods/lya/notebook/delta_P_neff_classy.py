def delta_P_neff(self, kp=0.009, zp=3.0, npoints=300):
    """
    Gives the dimensionless linear matter power spectrum and its logarithmic
    derivative at a pivot wavenumber in s/km and redshift.
    """

    # Compute pivot redshift in h/Mpc
    Omegam = self.Omega_m()
    h      = self.h()
    k_piv  = kp * np.sqrt(Omegam*(1.+zp)**3.+1.-Omegam)/(1.+zp)*100

    # Get linear matter power spectrum
    k_range = np.linspace(0.9*k_piv, 1.1*k_piv, npoints)
    Pk_lin  = np.array([self.pk_lin(k*h, zp)*h**3 for k in k_range])

    # Compute Neff
    ln_k_piv = np.log(k_piv)
    ln_P_ln_k_spline = UnivariateSpline(np.log(k_range),np.log(Pk_lin))
    neff = ln_P_ln_k_spline.derivative(1)(ln_k_piv)[()]

    # Compute delta
    delta = k_piv**3.*np.exp(ln_P_ln_k_spline(ln_k_piv))/(2.*np.pi**2.)

    return delta, neff


"""
This next function is for CAMB. 
"""
def get_delta_neff(kp=0.009,zp=3.0,header=camb_params,npoints=1000):
    ### kp     -- pivot scale in [s/km]
    ### zp     -- pivot redshift
    ### header -- contains cosmology parameters
    ### return:
    ###          Delta -- k^3 P_L(k)/(2pi^2) at (kp,zp)
    ###          neff  -- dlnP_L/dlnk at (kp,zp)
    ### set up the redshifts and parameters
    pars = camb.CAMBparams()
    redshifts=[zp]
    
    H0 = header['h0']*100.
    ombh2 = header['Omega_b']*header['h0']**2.
    omch2 = header['Omega_c']*header['h0']**2
    ns = header['ns']
    As = header['As']
    Omegam = header['Omega_b']+header['Omega_c']
    tau = header['tau']

    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, tau=tau, num_massive_neutrinos=0, mnu=0)
    pars.InitPower.set_params(ns=ns,As=As) ## VI: this could probably be optmized for As
    results = camb.get_results(pars)
    kpiv = kp * results.hubble_parameter(zp)/(1.+zp)/(H0/100.) ## convert to h/Mpc
    #kpiv = kp * np.sqrt(Omegam*(1.+zp)**3.+1.-Omegam)/(1.+zp)*100. ## faster covnersion?
    pars.set_matter_power(redshifts=redshifts, kmax=kpiv*1.5)
    #Linear spectra  
    pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)
    kh, z, pk_lin = results.get_matter_power_spectrum(minkh=0.8*kpiv,
                                                      maxkh=kpiv*1.5, npoints = npoints)
    lnkpiv = np.log(kpiv)
    splin = interp1d(np.log(kh),np.log(pk_lin),kind='linear')
    Delta = kpiv**3.*np.exp(splin(lnkpiv))/(2.*np.pi**2.)
    dlnk = np.log(kh[-1]/kh[0])/npoints
    neff = (splin(lnkpiv+dlnk)-splin(lnkpiv))/dlnk
    return Delta[0],neff[0]