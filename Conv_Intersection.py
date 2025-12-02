# -*- coding: utf-8 -*-

# Base packages
import numpy  as np
from mpmath import polylog

def to_doubled_array(a):
    """
    - If `a` is an np.ndarray: return it unchanged.
    - If `a` is a number: return a 1-D float64 array of length 1.
    - If `a` is a list (or tuple): return it as a NumPy array.
    """
    if isinstance(a, np.ndarray):
        return a
    if np.isscalar(a):
        return np.array([a], dtype=np.float64)
    if isinstance(a, (list, tuple)):
        return np.asarray(a).astype(np.float64)
    raise TypeError(f"Unsupported type: {type(a)!r}")
    
# Compute angles
def s_np(a, θ):
    þ = 1e-13    # Threshold for equality
    
    # Stable computation of s(θ)
    sinθ, cosθ = np.sin(θ), np.cos(θ)

    if θ == 0.0:
        sinθ, cosθ = 0.0, 1.0
    elif 2.0 * θ == np.pi:
        sinθ, cosθ = 1.0, 0.0
    elif 2.0 * θ == 3.0 * np.pi:
        sinθ, cosθ = -1.0, 0.0
    elif θ == np.pi:
        sinθ, cosθ = 0.0, -1.0
    elif θ == 2.0 * np.pi:
        sinθ, cosθ = 0.0, 1.0
    elif (np.sin(θ) * a == 1.0) and (a > 1.0):
        if cosθ < 0.0:
            sinθ, cosθ = 1.0 / a, -np.sqrt(a * a - 1.0) / a      # π - α
        else:
            sinθ, cosθ = 1.0 / a,  np.sqrt(a * a - 1.0) / a      # α
    elif (np.sin(θ) * a == -1.0) and (a > 1.0):
        sinθ, cosθ = -1.0 / a, np.sqrt(a * a - 1.0) / a         # 2π - α

    R = 1.0 - (a * a) * (sinθ * sinθ)
    if (abs(np.sin(θ)) * a == 1.0) and (a > 1.0):
        R = 0.0

    # Trim under precision
    if R < -þ:
        return 0.0

    # Clip and process
    if abs(R) <= þ:
        s_out = -a * cosθ
    else:
        s_out = -a * cosθ + np.sqrt(R)

    if a == 1.0:
        if cosθ >= 0.0:
            s_out = 0.0
        else:
            s_out = -2.0 * cosθ

    return s_out
def f_np(a, θ):
    s_aθ = s_np(a, θ)
    if s_aθ == 0.0:
        return 0.0
    return (s_aθ * s_aθ) * (2.0 * np.log(abs(s_aθ)) - 1.0)
def φ_from_np(a, ε):
    þ = 1e-13    # Threshold for equality

    cos  = 1.0 - a * a - ε * ε
    cos /= 2.0 * a * ε

    # Special cases
    if a == 1.0 - ε:              # cos φ = 1
        return 0.0
    elif a == 1.0 + ε:            # cos φ = -1
        return np.pi
    elif np.isclose(a * a, 1.0 + ε * ε, atol=þ, rtol=þ):
        # cos φ = -ε/sqrt(1 + ε^2)
        return np.arctan(ε) + 0.5 * np.pi
    elif np.isclose(a * a, 1.0 - ε * ε, atol=þ, rtol=þ):
        # cos φ = 0
        return 0.5 * np.pi
    elif a == 1.0:                # cos φ = -ε/2
        return np.arccos(-0.5 * ε)

    return np.arccos(cos)
def Φ_np(a, θ):
    '''
        Variable transformation
    '''
    sθ = s_np(a, θ)
    Lθ = (sθ * sθ - 1.0 - a * a) / (2.0 * a)      # Range [-1,1]
    if a == 1.0:
        Lθ = (sθ * sθ * 0.5) - 1.0
    ϕ = 0.5 * np.arccos(Lθ)                       # Range [0, π/2]
    return ϕ

# Define Li_2(z) (extends beyond 1 from spence)
_Li_2 = np.frompyfunc(lambda z: polylog(2, z), 1, 1)
Li_2  = lambda z: np.asarray(_Li_2(z), dtype=np.complex128)

# Additional integrals
def F_np(a, ε):
    þ = 1e-13    # Threshold for equality

    '''
        Variable transformation
    '''
    θ = φ_from_np(a, ε)
    ϕ = Φ_np(a, θ)
    
    if a > 1.0:
        # α activates
        α = np.arcsin(1.0 / a)
        
        # Exact evaluation at breaking point
        if np.isclose(a * a, 1.0 + ε * ε, atol=þ, rtol=þ):
            ϕ = 0.5 * (np.pi - np.arctan(ε))
            α = (0.5 * np.pi) - np.arctan(ε)
        
        # Auxiliary angle = Φ(α)
        ω = (2.0 * α + np.pi) * 0.25
    
    '''
        Functional evaluation
    '''
    if np.isclose(a, 1.0 - ε, atol=þ, rtol=þ):
        # At interval extremes, the mass is zero
        # π ε (2 - ε) cancels out the quadratic term
        T = 0.0
    elif (a > 1.0 - ε) and (a < 1.0):
        # One continuous interval of existence
        T = G_np(a, ϕ) - (1.0 - a * a) * np.pi
    elif np.isclose(a, 1.0, atol=þ, rtol=þ):
        # One functional evaluation
        T = G_np(1.0, Φ_np(1.0, φ_from_np(1.0, ε)))
    elif (a > 1.0) and (a * a < 1.0 + ε * ε):
        # Two intervals of existence and correct branch selection
        T = G_np(a, ϕ) - (2.0 * G_np(a, ω)) - (2.0 * np.pi * np.log(a))
    elif (a * a >= 1.0 + ε * ε):
        # One interval of continuity, no α involved
        ϕ = Φ_np(a, θ + np.pi)
        if np.isclose(a * a, 1.0 + ε * ε, atol=þ, rtol=þ):
            ϕ = (np.pi - np.arctan(ε)) * 0.5
        # Principal boundary value (approached from the upper half-plane):  
        # 2 im Li_2 (a + 0i) = -2π a
        T = -2.0 * np.pi * np.log(a) - G_np(a, ϕ)
        if np.isclose(a, 1.0 + ε, atol=þ, rtol=þ):
            # Φ(φ+π) = π/2 and integral cancels out in the right branch
            T = 0.0
    
    # Scale
    T /= np.pi * 8.0
    
    return T
def G_np(a, ϕ):
    þ = 1e-13    # Threshold for equality

    '''
        Angular quantities
    '''
    cosϕ, cos2ϕ, sin2ϕ  = np.cos(ϕ), np.cos(2.0 * ϕ), np.sin(2.0 * ϕ)
    if np.isclose(2.0 * ϕ, np.pi, atol=þ, rtol=þ):
        cosϕ, cos2ϕ, sin2ϕ = 0.0, -1.0, 0.0
    if np.isclose(ϕ, 0.0, atol=þ, rtol=þ):
        cosϕ, cos2ϕ, sin2ϕ = 1.0, 1.0, 0.0

    '''
        Special cases
    '''
    # G(1; ϕ)
    if np.isclose(a, 1.0, atol=þ, rtol=þ):
        if np.isclose(cosϕ, 0.0, atol=þ, rtol=þ):
            return 0.0
        else:
            z  = -np.exp(2j * ϕ)
            G  = 2.0 * Li_2(z).imag
            G += 2.0 * (1.0 - np.log(2.0 * abs(cosϕ))) * sin2ϕ
            return G

    # G(a; 0)
    if np.isclose(ϕ, 0.0, atol=þ, rtol=þ):
        return 0.0

    # G(a; (2α+π)/2 )
    if a > 1.0:
        α = np.arcsin(1.0 / a)
        # Exact case: 4ϕ == 2α + π
        if np.isclose(ϕ * 4.0, 2.0 * α + np.pi, atol=þ, rtol=þ):
            # Dilog term
            z   = -a * np.exp(2j * ϕ)
            T_1 = 2.0 * Li_2(z).imag
            # Angular term
            T_2 = (1.0 - a * a) * α
            # Log term
            T_3 = (2.0 - np.log(a * a - 1.0)) * np.sqrt(a * a - 1.0)
            return T_1 + T_2 + T_3

    '''
        General case
    '''
    # Dilog term
    z   = -a * np.exp(2j * ϕ)
    T_1 = 2.0 * Li_2(z).imag

    # Angular term
    ang = np.angle(1.0 - z)
    T_2 = (1.0 - a * a) * (2.0 * ϕ - ang)

    # Log term (stable version of 1 + a^2 + 2 a cos 2ϕ)
    l_a = (a - 1.0) * (a - 1.0) + 4.0 * a * (cosϕ * cosϕ)
    T_3 = a * (2.0 - np.log(l_a)) * sin2ϕ
    if np.isclose(sin2ϕ, 0.0, atol=þ, rtol=þ):
        T_3 = 0.0

    # Collect terms:
    G = T_1 + T_2 + T_3
    
    return G

# Stable version of E
def E_np(a,ε, asymp = False):
    """
    Stable E(a; ε) for small ε, with a in [0,+∞).

    Parameters
    ----------
    a : float or np.ndarray
        Radial variable, scalar or array.
    ε : float
        Small parameter, typically ε < 0.5.
    asymp : bool, optional
        If True, use asymptotic expansions (recommended for ε < 1e-6).
    """
    # Preprocess a, ε
    a = to_doubled_array(a)
    ε = to_doubled_array(ε).item()
    
    '''
        Determine evaluation cases
    '''
    masks = [(a <= 1.0 - ε), (a >= 1.0 + ε)]
    masks.append(~masks[0] & ~masks[1])
    
    E = np.zeros_like(a)    # Already covers a ≥ 1 + ε
    
    if  masks[0].any():    # a ≤ 1 - ε:
        E[masks[0]] = np.full(masks[0].sum(), ε * ε * (np.log(ε*ε) - 1))
        
    if masks[2].any():
        # Use asymptotic expansion instead of exact formulae
        if asymp:
            # Map a to [0,1]
            λ = (a[masks[2]] - (1.0 - ε)) / (2.0 * ε)
            
            # Approximate angle
            φ = np.acos(1.0 - 2.0 * λ) + np.sqrt(λ * abs(1.0 - λ))*(ε + 0.75 * (1.0 - 2.0 * λ)*ε*ε)
            #φ = np.asarray([φ_from_np(a_it, ε) for a_it in a[masks[2]]])
            E[masks[2]] = (np.pi - φ) * ε * ε * (np.log(ε*ε) - 1.0) / np.pi
            
            # Branch points
            λ_masks = [(λ < 0.5), (λ > 0.5)]
            # *** a < 1 ***
            λ_L = λ[λ_masks[0]]
            h_1 = 0.25 * (1.0 - 2.0 * λ_L) * np.sqrt(λ_L * (1.0 - λ_L)) / np.pi
            h_2 = 0.25 * (1.0 - 2.0 * λ_L) * ( (1.0 - 2.0 * λ_L) * np.acos(1.0 - 2.0 * λ_L) - 3.0 * np.sqrt(λ_L * (1.0 - λ_L))) / np.pi
            
            F = h_1 * (2.0 * ε*ε*np.log(ε)) + h_2 * (ε * ε)
            up_mask = np.zeros_like(masks[2], dtype = bool)
            up_mask[masks[2]] = λ_masks[0]                        # True only where masks[2] is True AND λ < 0.5
            
            E[up_mask] += 8.0 * F
            
            
            # *** a = 1 ***
            þ = 1e-13    # Threshold for equality
            λ_1i = np.isclose(λ, 0.5, atol=þ, rtol=þ)
            if λ_1i.any():
                F_1 = ( 6.0 * ε * ε * ε * np.log(ε) - 5.0 * ε * ε * ε ) / (18.0 * np.pi)                # 144 = 18 * 8
                E[masks[2]][λ_1i] = np.full( λ_1i.sum(), F_1 )
                
            # *** a > 1 ***
            # Use approximation from symmetry
            λ_U = 1.0 - λ[λ_masks[1]]
            
            h_1 = 0.25 * (1.0 - 2.0 * λ_U) * np.sqrt(λ_U * (1.0 - λ_U)) / np.pi
            h_2 = 0.25 * (1.0 - 2.0 * λ_U) * ( (1.0 - 2.0 * λ_U) * np.acos(1.0 - 2.0 * λ_U) - 3.0 * np.sqrt(λ_U * (1.0 - λ_U))) / np.pi
            
            F = h_1 * (2.0 * ε*ε*np.log(ε)) + h_2 * (ε * ε)
            up_mask = np.zeros_like(masks[2], dtype = bool)
            up_mask[masks[2]] = λ_masks[1]                        # True only where masks[2] is True AND λ < 0.5
            E[up_mask] += 8.0 * F
            
        # Use exact formulae
        else:
            φ = np.asarray([φ_from_np(a_it, ε) for a_it in a[masks[2]]])
            E[masks[2]] = (np.pi - φ) * ε * ε * (np.log(ε*ε) - 1.0) / np.pi
            
            F = np.asarray([F_np(a_it,ε) for a_it in a[masks[2]]])
            E[masks[2]] += 8.0 * F
        
    # Scale E
    E *= 0.25 
    
    return E

__all__ = ["E_np"]