######FLIGHT PATH ESTIMATOR######

#Units: kg, km, s, kN#

#ASSUMPTIONS#
#CONSTANT THRUST
#CONSTANT MASS FLOW RATE
#2D FLIGHT PATH
#INTERNATIONAL STANDARD ATMOSPHERIC MODEL

#Import Dependencies#
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import numpy as np
from scipy.integrate import solve_ivp

#CONSTANTS#

#Environmental Parameters#

g0 = 0.00981                #Initial gravitational acceleration (km/s^2)
rho0 = 1225000000           #Atmospheric density at sea level (kg/km^3)
p0 = 1                      #Static pressure at sea level (atm)
R_E = 6378                  #Radius of Earth in (km)

#MISSION PARAMETERS#
h_pitch = 0.01                              #Altitude to start gravity turn
gamma0 = np.radians(88)                     #Initial pitch angle at launch
v0 = x0 = h0 = vd0 = vg0 = 0                #Initial state-vector values
y0 = [v0, gamma0, x0, h0, vg0, vd0]         #Initial state-vector
t0 = 0                                      #Initial time
tf = 400                                    #Final time
t_span = (t0,tf)                            #Time span for RK45 function

#VEHICLE PARAMETERS#
stage_params = [

#dictionary of staging parameters:

    {"stage": 0, "m0":928.7,"mProp": 120,
    "radius": 0.00019, "t_burnout": 2.5, "thrustASL": 78.6, "thrustVac": 78.5},

    {"stage": 1, "m0":658,"mProp": 479.3,
    "radius": 0.00019, "t_burnout": 52, "thrustASL": 20.3, "thrustVac": 20.3},

    {"stage": 2, "m0":178.2,"mProp": 0.0,
    "radius": 0.00019, "t_burnout": 1000, "thrustASL": 0.0, "thrustVac": 0.0}
]

def get_grav(h):

    """
    Calculates gravitational acceleration based on altitude
    for use in calculating gravity losses

    :param h: altitude in km (float or int)
    :return: gravitational acceleration in km/s^2

    """
    if h < 0:
        raise ValueError("h must be non-negative")

    got_grav = g0 / ((1 + (h/R_E)) ** 2)
    return got_grav

def get_rho(h):

    """
    Calculates density based on ISA/NRLMSISE-00 atmospheric models using altitude(h) as an input
    for use in drag calculation

    :param h: altitude in km (int or float)
    :return: density in kg/km^3
    """

    if h < 0:
        raise ValueError("h must be non-negative")


    if h <= 100:
        got_rho = rho0 * np.exp((-h)/get_h_scale(h)) #exponential model
        return got_rho

    else:
        got_rho = get_rho(100) * np.exp(-(h-100)/60) #NRLMSISE-00 model
        return got_rho

def get_h_scale(h):

    """
    interpolates to find h_scale
    for use in density and static pressure calculation

    :param h: altitude in km (int or float)
    :return: h_scale
    """

    if h < 0:
        raise ValueError("h must be non-negative")

    values = [7.48, 6.33, 6.73, 7.48, 7.92, 5.86, 5.42, 6.5]
    altitudes = [0, 11, 20, 32, 47, 51, 71, 200]
    got_h_scale = interp1d(altitudes, values, kind='linear', fill_value='extrapolate')
    return got_h_scale(h)

def get_speed_of_sound(h):

    """
    uses interpolation to estimate speed of sound at a given altitude (h)
    for use in Mach number calculation

    :param h: altitude in km (int or float)
    :return: speed of sound in km/s
    """

    if h < 0:
        raise ValueError("h must be non-negative")

    values = [0.340, 0.320, 0.300, 0.395, 0.395, 0.305, 0.315, 0.320, 0.330, 0.335, 0.340, 0.340, 0.340]
    altitudes = [0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 200]
    got_speed_of_sound = interp1d(altitudes, values, kind='linear', fill_value='extrapolate')

    return got_speed_of_sound(h)

def get_stat_pressure(h):

    """
    gets static pressure using altitude (h) for use in determining current engine thrust
    :param h: altitude in km
    :return: static pressure (atm)
    """

    if h < 0:
        raise ValueError("h must be non-negative")

    if h <= 100:
        got_stat_pressure = p0 * np.exp((-h)/get_h_scale(h)) #exponential model
        return got_stat_pressure
    else:
        got_stat_pressure = get_stat_pressure(100) * np.exp(-(h-100)/60) #NRLMSISE-00 model
        return got_stat_pressure

def get_cd(mach):

    """
    uses bisect lists to find drag coefficient (c_d) for a given Mach number (M)
    for use in drag calculation

    :param mach: Mach number (int or float)
    :return: c_d (dimensionless)
    """

    mach_numbers = [0.00, 0.75, 1.0, 1.10, 1.25, 2.00, 2.25, 3.00, 4.00, 8.0]
    cd_values =    [0.55, 0.55, 0.7, 0.95, 0.95, 0.65, 0.50, 0.45, 0.41, 0.4]

    # Create interpolation function
    cd_interp = interp1d(mach_numbers, cd_values, kind='linear', fill_value='extrapolate')

    got_cd = cd_interp(abs(mach))

    return got_cd

def get_mach(h,v):

    """
    gets current Mach number based on velocity(v) and altitude(h)
    :param h: altitude (h) in km
    :param v: velocity (v) in km/s
    :return:
    """
    if h < 0:

        raise ValueError("h must be non-negative")

    got_mach = v / get_speed_of_sound(h)
    return got_mach

def get_stage(t):
    """
    gets current stage
    for use in stage dependent parameters (drag, thrust)
    :param t: time in seconds
    :return: current stage number
    """
    found_stage = False
    for p in stage_params:
        if not found_stage and t < p['t_burnout']:
            found_stage = True
            got_stage = p['stage']
            return got_stage


def get_mass(t, stage):
    """
    Defines mass at time t to account for propellant burn and staging
    for use in equations of motion

    Args:

    :param t: time in seconds (float or int)
    :return: mass in kg
    """

    if t < 0:
        raise ValueError("t must be non-negative")

    s = stage_params[stage]
    stage_time = t
    for i in range(stage - 1):
        stage_time = stage_time - stage_params[i]['t_burnout']

    got_mass = s['m0'] - (s['mProp'] / s['t_burnout']) * stage_time
    return got_mass


def get_drag(v, h, stage):
    """
    Defines Drag as a function of velocity (v), altitude (h), and drag coefficient (c_d)
    for use in equations of motion

    :param v: velocity (float or int)
    :param h: altitude (float or int)
    :param stage: current stage (float or int)
    :return: Drag force in kN
    """

    if h < 0:
        raise ValueError("h must be non-negative")

    m = get_mach(h, v)
    area = np.pi * stage_params[stage]['radius'] ** 2

    got_drag = 0.5 * get_cd(m) * get_rho(h) * v ** 2 * area
    return got_drag

def get_thrust(t, h, stage):

    """
    calculates thrust based on the thrust at sea level, thrust in a vacuum, and atmospheric pressure

    :param t: time in seconds
    :param h: altitude in km
    :param stage: current stage
    :return: thrust in kN
    """

    if t < 0:
        raise ValueError("t must be non-negative")
    if h < 0:
        raise ValueError("h must be non-negative")

    s = stage_params[stage]
    got_thrust = s['thrustASL'] + (s['thrustVac'] - s['thrustASL']) * (1-get_stat_pressure(h))
    return got_thrust

#####Simulation####

#Sim function to integrate state-vector (y) forward in time (t) using RK4(5) method#
#State vector: y = velocity (v), pitch angle (gamma), downstream displacement (x), altitude (h),
#drag losses (d_losses), and gravity losses (g_losses).

def rates(t, y):

    """
    Defines equations for state-vector rates in addition to drag and gravity losses.
    for use in integrator solver

    :param t: time in seconds (float or int)
    :param y: state-vector containing 6 elements - velocity, pitch angle, downstream displacement,
              altitude, drag loss, gravity loss (list, array, or tuple)
    :return: state-vector rates
    """

    if len(y) != 6:
        raise ValueError("state-vector y must contain 6 elements")

    if t < 0:
        raise ValueError("t must be non-negative")

    v,gamma,x, h, d_losses, g_losses = y

    stage_flag = get_stage(t)
    T = get_thrust(t, h, stage_flag)
    m = get_mass(t, stage_flag)
    drag = get_drag(v,h,stage_flag)
    g = get_grav(h)

    if h <= h_pitch:          #Pitching assumed to start at h_pitch. Until then, assume vertical flight

        dv_dt = T / m - drag / m - g  #Acceleration (from Newton's second law)
        dgamma_dt = 0                            #Rate of pitching due to gravity
        dx_dt = 0                                #Downstream velocity
        dh_dt = v                                #Vertical velocity component
        dg_losses_dt = g                      #Rate of losses due to gravity

    else:

        dv_dt = T / m - drag / m - g * np.sin(gamma)             #Acceleration (from Newton's second law)
        dgamma_dt = - (1 / v) * (g - (v ** 2 / (R_E + h))) * np.cos(gamma)   #Rate of pitching due to gravity
        dx_dt = R_E / (R_E + h) * v * np.cos(gamma)                             #Downstream velocity
        dh_dt = v * np.sin(gamma)                                               #Vertical velocity component
        dg_losses_dt = g * np.sin(gamma)                                     #Rate of losses due to gravity

    dd_losses_dt = drag / m              #Rate of drag losses

    return dv_dt, dgamma_dt, dx_dt, dh_dt, dd_losses_dt, dg_losses_dt

#####RK4(5) and solution######

sol = solve_ivp(rates, t_span, y0, method='RK45', rtol=1e-8, atol=1e-10, max_step=0.1)

v_final = sol.y[0, -1]                      #Final velocity
gamma_final = sol.y[1, -1]                  #Final pitch angle
x_final = sol.y[2, -1]                      #Final downstream displacement
h_final = sol.y[3, -1]                      #Final altitude
d_losses_final = sol.y[4, -1]               #Total drag losses
g_losses_final = sol.y[5, -1]               #Total gravity losses

time = sol.t                                #Time data
velocity = sol.y[0]                         #Velocity data
gamma = sol.y[1]                            #Gamma data
range_km = sol.y[2]                         #X_disp data
altitude = sol.y[3]                         #Altitude data
drag = sol.y[4]
grav = sol.y[5]

#Plotting Data

stage_data = [get_stage(t) for t in time]
thrust_plot_data = [get_thrust(t, h, stage) for t, h, stage in zip(time, altitude, stage_data)]
dynamic_pressure_data = [(0.5 * get_rho(h) * v ** 2) for h,v in zip(altitude, velocity)]
mass_data = [get_mass(t,stage) for t, stage in zip(time, stage_data)]
flight_path_angle_degs = np.rad2deg(gamma)

fig, results = plt.subplots(2, 2, figsize=(10, 8))

# Plots
results[0, 0].plot(time, altitude, label="Altitude")
results[0, 1].plot(time, range_km, label="Range")
results[1, 0].plot(time, velocity, label="Velocity")
results[1, 1].plot(time, flight_path_angle_degs, label="Flight Path Angle")

# Plot settings
titles = ["Altitude (km)", "Range (km)", "Velocity (km/s)", "Flight Path Angle (deg)"]
ylims = [[0,290], [0,80],[0,1.5], [-90,90]]
for ax, title, ylim in zip(results.flat, titles, ylims):
    ax.set_title(title)
    ax.set_ylim(ylim)
    ax.set_xlim([0, tf])

    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()

#PRINT FINAL VALUES#

print(f"Final velocity: {v_final:.2f} km/s")
print(f"Final altitude: {h_final:.2f} km")
print(f"Total gravity losses: {g_losses_final:.2f} km/s (Δv)")
print(f"Total drag losses: {d_losses_final:.2f} km/s (Δv)")
print(f"Final flight path angle: {np.degrees(gamma_final)} deg")
print(f"Max Altitude Reached: {max(altitude)} km")
print(f"Max Velocity Reached: {max(velocity)} km/s")












