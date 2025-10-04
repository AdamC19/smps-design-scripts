
from enum import Enum
import numpy as np
import matplotlib.pyplot as plt
import control
from dataclasses import dataclass

import warnings
warnings.filterwarnings('ignore')

def pole(freq) -> control.TransferFunction:
    """Returns a transfer function for a pole at given angular freq (rad/s)"""
    return control.tf([1.0],[1.0/freq, 1])

def zero(freq) -> control.TransferFunction:
    """Returns a transfer function for a zero at given angular freq (rad/s)"""
    return control.tf([1.0/freq, 1.0], [1.0])

def rhp_zero(freq) -> control.TransferFunction:
    """Returns a transfer function for a Right-Half-Plane (inverted) zero at 
    given angular freq (rad/s)
    """
    return control.tf([-1.0/freq, 1.0], [1.0])


# def buck_output_filter(ind, cout, )
class SmpsType(Enum):
    BUCK_VOLTAGE_MODE = 'buck, voltage mode'
    BUCK_CURRENT_MODE = 'buck, current mode'
    BOOST_VOLTAGE_MODE = 'boost, voltage mode'
    BOOST_CURRENT_MODE = 'boost, current mode'

@dataclass
class SmpsInfo:
    conv_type: SmpsType
    fsw: float
    vin: float
    vout: float
    iout: float
    Lout: float
    Cout: float
    R_Cout_esr: float
    R_sns: float
    i_sns_gain: float
    vslope: float # slope value leading to 50% duty cycle
    eff: float


def buck_vmode(info: SmpsInfo) -> control.TransferFunction: # L_out, C_out, R_cout_esr, R_out: float = 10.0) -> control.TransferFunction:
    """Builds a transfer function for a voltage-mode buck."""
    w_zero = 1.0/(info.Cout * info.R_Cout_esr)
    w_0 = 1.0/(np.sqrt(info.Lout*info.Cout))
    R_out = info.vout / info.iout
    Q_0 = R_out/(np.sqrt(info.Lout/info.Cout))
    A_vc = info.vin/info.vslope
    return A_vc * control.tf([0.0, 1.0/w_zero, 1.0],[1.0/(w_0**2), 1.0/(Q_0*w_0), 1.0])

def buck_imode(info: SmpsInfo) -> control.TransferFunction: # L_out, C_out, R_cout_esr, R_out: float = 10.0, k_m: float = 4.0) -> control.TransferFunction:
    w_z = 1.0/(info.Cout * info.R_Cout_esr)
    R_out = info.vout / info.iout
    R_i = info.i_sns_gain * info.R_sns
    w_p = 1.0/(info.Cout * R_out)
    V_slope_opt = (info.vout * R_i)/(info.fsw * info.Lout)
    k_m = info.vin / V_slope_opt
    w_L = (k_m * R_i) / info.Lout

    A_vc = R_out/R_i

    xfer_zero = zero(w_z)
    xfer_pole_1 = pole(w_p)
    xfer_pole_2 = pole(w_L)

    return A_vc*xfer_zero*xfer_pole_1*xfer_pole_2


def boost_imode(info: SmpsInfo) -> control.TransferFunction:
    D = (info.vout - info.vin) / info.vout
    Dp = info.vin/info.vout
    R_out = info.vout / info.iout
    R_i = info.i_sns_gain * info.R_sns
    V_slope_opt = ((info.vout - info.vin) * R_i)/(info.fsw * info.Lout)
    k_m = info.vout / V_slope_opt

    A_vc = (R_out*Dp)/(2*R_i)
    w_p = 2/(info.Cout*R_out)
    w_L = (k_m * R_i)/info.Lout
    w_rhpz = (R_out * Dp**2)/info.Lout
    w_z = 1/(info.R_Cout_esr*info.Cout)

    rhpz = rhp_zero(w_rhpz)
    zero_1 = zero(w_z)
    pole_1 = pole(w_p)
    pole_2 = pole(w_L)

    return A_vc*rhpz*zero_1*pole_1*pole_2


def error_amp_type1(R_in, C_comp) -> control.TransferFunction:
    return control.tf([0.0, -1.0], [R_in*C_comp, 0.0])

def error_amp_type2(R_in, R_comp, C_comp, C_hf) -> control.TransferFunction:
    A_vm = R_comp/R_in # midband voltage gain
    w_zea = 1.0/(R_comp*C_comp)
    w_hf = 1.0/(R_comp*C_hf)

    return control.tf([1.0], [-1.0/(A_vm*w_zea), 0.0])*zero(w_zea)*pole(w_hf)

def error_amp_type3(R_in, R_comp, C_comp, C_hf, R_ff, C_ff) -> control.TransferFunction:
    A_vm = R_comp/R_in # midband voltage gain
    w_zea = 1.0/(R_comp*C_comp)
    w_fz = 1/(R_in*C_ff)
    w_fp = 1/(R_ff*C_ff)
    w_hf = 1.0/(R_comp*C_hf)

    return control.tf([0.0, -1.0*A_vm*w_zea], [1.0, 0.0]) * zero(w_zea) * zero(w_fz) * pole(w_fp) * pole(w_hf)


def calc_duty_cycle(info: SmpsInfo) -> float:
    """Calculate the steady state duty cycle for the given converter."""
    if info.conv_type == SmpsType.BUCK_CURRENT_MODE or info.conv_type == SmpsType.BUCK_VOLTAGE_MODE:
        return info.vout/(info.vin * info.eff)
    elif info.conv_type == SmpsType.BOOST_CURRENT_MODE or info.conv_type == SmpsType.BOOST_VOLTAGE_MODE:
        return 1.0 - (info.eff*info.vin)/info.vout
    

def calc_ind_ripple(info: SmpsInfo) -> float:
    """Calculate the peak-to-peak inductor ripple current in steady state for the given converter."""
    duty = calc_duty_cycle(info)
    Ton = duty/info.fsw
    if info.conv_type == SmpsType.BUCK_CURRENT_MODE or info.conv_type == SmpsType.BUCK_VOLTAGE_MODE:
        vind = info.vin - info.vout
    elif info.conv_type == SmpsType.BOOST_CURRENT_MODE or info.conv_type == SmpsType.BOOST_VOLTAGE_MODE:
        vind = info.vin
    
    return (Ton * vind) / (info.Lout)

def calc_ind_ipeak(info: SmpsInfo) -> float:
    """Calculate the peak (max) inductor ripple current in steady state for the given converter."""
    i_rip = calc_ind_ripple(info)

    if info.conv_type == SmpsType.BUCK_CURRENT_MODE or info.conv_type == SmpsType.BUCK_VOLTAGE_MODE:
        return info.iout + i_rip/2.0
    elif info.conv_type == SmpsType.BOOST_CURRENT_MODE or info.conv_type == SmpsType.BOOST_VOLTAGE_MODE:
        return info.iout/(1.0 - calc_duty_cycle(info)) + i_rip/2.0