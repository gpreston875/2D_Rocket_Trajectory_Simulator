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

g0 = 0.00981                #Initial gravitational acceleration km/s^2
rho0 = 1225000000           #Atmospheric density at sea level (kg/km^3)
p0 = 1                      #Static pressure at sea level (atm)
R_E = 6378                  #

#MISSION PARAMETERS#
h_pitch = 0.01                              #Altitude to start gravity turn
gamma0 = np.radians(90)                     #Initial pitch angle at launch
v0 = x0 = h0 = vd0 = vg0 = 0                #Initial state-vector values
y0 = [v0, gamma0, x0, h0, vg0, vd0]         #Initial state-vector
t0 = 0                                      #Initial time

#VEHICLE PARAMETERS#
stage_flag = 0
stage_params = [

#dictionary of staging parameters:

    {"stage": 0, "m0":623,"mProp": 66.2,
    "radius": 0.00015, "t_burnout": 0.53, "thrustASL": 222.4, "thrustVac": 244.4},

    {"stage": 1, "m0":342.8,"mProp": 191,
    "radius": 0.00015, "t_burnout": 47, "thrustASL": 6.7, "thrustVac": 7.7},

    {"stage": 2, "m0":55,"mProp": 0.0,
    "radius": 0.00015, "t_burnout": 1000, "thrustASL": 0.0, "thrustVac": 0.0}
]


#Stage - 1

area_s1 = np.pi * (0.00015 ** 2)     #Rocket tip area
m0_s1 = 623                                    #Initial mass
m_prop_s1 = 66.2                               #Propellant mass
mf_s1 = 623 - m_prop_s1                        #Final stage-1 mass (before decoupling)
ThrustASL_s1 = 222.4                           #Thrust at sea level (kN)
ThrustVac_s1 = 244.4                           #Thrust in a vacuum (kN)
t_bo_s1 = 0.63                                 #Burnout time stage1 (s)
m_dot_s1 = m_prop_s1/t_bo_s1                   #Propellant burn rate

#Stage - 2

area_s2 = np.pi * (0.00015 ** 2)               #Rocket tip area
m0_s2 = 623 - 214 - 66.2                       #Initial mass
m_prop_s2 = 191                                #Propellant mass
mf_s2 = 155                                    #Mass final
ThrustASL_s2 = 6.7                             #Thrust at sea level (kN)
ThrustVac_s2 = 7.7                             #Thrust in a vacuum (kN)
t_bo_s2 = 45                                   #Burnout time
m_dot_s2 = m_prop_s2/t_bo_s2

def get_stage(t):
    found_stage = False
    for p in stage_params:
        if not found_stage and t < p['t_burnout']:
            found_stage = True
            return p['stage']

def mass(t, num_stages):

    """
    Defines mass at time t to account for propellant burn and staging

    Args:

    :param t: time in seconds (float or int)
    :param num_stages: number of stages
    :return: mass at time t
    """

    if t < 0:
        raise ValueError("t must be non-negative")

    # TODO: Cascade get stages function throughout




    if num_stages == 1 or t <= t_bo_s1:
        return m0_s1 - m_dot_s1 * t
    elif num_stages == 1 and t > t_bo_s1:
        return m0_s1 - m_dot_s1 * t_bo_s1

    elif num_stages  == 2 and t <= t_bo_s2:
        return m0_s2 - m_dot_s2 * (t - t_bo_s1)
    elif num_stages == 2 and t > t_bo_s2:
        return m0_s2 - m_dot_s2 * t_bo_s2

def D(v, h, t, num_stages):

    """
    Defines Drag as a function of velocity (v), altitude (h), and drag coefficient (c_d)

    :param v: velocity (float or int)
    :param h: altitude (float or int)
    :param t: time in seconds
    :param num_stages: number of stages
    :return: Drag force in kN
    """


    if h < 0:
        raise ValueError("h must be non-negative")

    m = mach(h, v)

    if num_stages == 1 or t <= t_bo_s1:
        return 0.5 * c_d(m) * rho(h) * v ** 2 * area_s1
    elif num_stages > 1 and t <= tf:
        return 0.5 * c_d(m) * rho(h) * v ** 2 * area_s2


def g(h):

    """
    Calculates gravitational acceleration based on altitude

    :param h: altitude in km (float or int)
    :return: gravitational acceleration in km/s^2

    """


    if h < 0:
        raise ValueError("h must be non-negative")

    return g0 / ((1 + (h/R_E)) ** 2)

def rho(h):

    """
    Calculates density based on ISA/NRLMSISE-00 atmospheric models using altitude(h) as an input

    :param h: altitude in km (int or float)
    :return: density in kg/km^3
    """

    if h < 0:
        raise ValueError("h must be non-negative")


    if h <= 100:
        return rho0 * np.exp((-h)/get_h_scale(h)) #exponential model
    else:
        return rho(100) * np.exp(-(h-100)/60) #NRLMSISE-00 model

def get_h_scale(h):

    """
    interpolates to find h_scale for use in density calculation at a given altitude (h)

    :param h: altitude in km (int or float)
    :return: h_scale
    """

    if h < 0:
        raise ValueError("h must be non-negative")


    values = [7.48, 6.33, 6.73, 7.48, 7.92, 5.86, 5.42, 6.5]
    altitudes = [0, 11, 20, 32, 47, 51, 71, 200]
    h_s = interp1d(altitudes, values, kind='linear', fill_value='extrapolate')
    return h_s(h)

def get_speed_of_sound(h):

    """
    uses interpolation to estimate speed of sound at a given altitude (h)


    :param h: altitude in km (int or float)
    :return: speed of sound in km/s
    """

    if h < 0:
        raise ValueError("h must be non-negative")

    values = [0.340, 0.320, 0.300, 0.395, 0.395, 0.305, 0.315, 0.320, 0.330, 0.335, 0.340, 0.340, 0.340]
    altitudes = [0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 200]
    SoS = interp1d(altitudes, values, kind='linear', fill_value='extrapolate')

    return SoS(h)

def get_stat_pressure(h):



    if h < 0:
        raise ValueError("h must be non-negative")


    if h <= 100:
        return p0 * np.exp((-h)/get_h_scale(h)) #exponential model
    else:
        return get_stat_pressure(100) * np.exp(-(h-100)/60) #NRLMSISE-00 model

def mach(h,v):

    if h < 0:

        raise ValueError("h must be non-negative")

    return v / get_speed_of_sound(h)

def c_d(M):

    """
    uses bisect lists to find drag coefficient (c_d) for a given Mach number (M)

    :param M: Mach number (int or float)
    :return: c_d
    """



    mach_numbers = [0, 1, 1.25, 2.0, 3.5, 8]
    cd_values = [0.05, 0.010, 0.024, 0.034, 0.036, 0.036]

    # Create interpolation function
    cd_interp = interp1d(mach_numbers, cd_values, kind='linear', fill_value='extrapolate')


    cd = cd_interp(abs(M))

    return cd

def thrust(t, h, num_stages, TASL1, TVac1, TASL2, TVac2, tbo1, tbo2):
    if t < 0:
        raise ValueError("t must be non-negative")
    if h < 0:
        raise ValueError("h must be non-negative")
    if num_stages < 1:
        raise ValueError("number of stages must be greater than 0")

    if t < tbo1: #If propellant still available

        return TASL1 + (TVac1 - TASL1) * (1-get_stat_pressure(h))


    elif num_stages > 1 and t < tbo2:

        return TASL2 + (TVac2 - TASL2) * (1-get_stat_pressure(h))


    else:

        return 0


#####Simulation####

#Sim function to integrate state-vector (y) forward in time (t) using RK4(5) method#
#State vector: y = velocity (v), pitch angle (gamma), downstream displacement (x), altitude (h),
#drag losses (d_losses), and gravity losses (g_losses).

drag_test = []
def rates(t, y):

    """
    Defines equations for state-vector rates in addition to drag and gravity losses.


    :param t: ime in s (float or int)
    :param y: state-vector containing 6 elements - velocity, pitch angle, downstream displacement,
              altitude, drag loss, gravity loss (list, array, or tuple)
    :return: state-vector rates
    """



    if len(y) != 6:
        raise ValueError("state-vector y must contain 6 elements")

    if t < 0:
        raise ValueError("t must be non-negative")


    v,gamma,x, h, d_losses, g_losses = y

    T = thrust(t, h, stages, ThrustASL_s1, ThrustVac_s1, ThrustASL_s2, ThrustVac_s2, t_bo_s1, t_bo_s2)


    m = mass(t, stages)

    if h <= h_pitch:          #Pitching assumed to start at h_pitch. Until then, assume vertical flight

        dv_dt = T / m - D(v, h, t, stages) / m - g(h)  #Acceleration (from Newton's second law)
        dgamma_dt = 0                            #Rate of pitching due to gravity
        dx_dt = 0                                #Downstream velocity
        dh_dt = v                                #Vertical velocity component
        dg_losses_dt = g(h)                      #Rate of losses due to gravity

    else:

        dv_dt = T/m - D(v, h, t, stages) / m - g(h) * np.sin(gamma)             #Acceleration (from Newton's second law)
        dgamma_dt = - (1 / v) * (g(h) - (v ** 2 / (R_E + h))) * np.cos(gamma)   #Rate of pitching due to gravity
        dx_dt = R_E / (R_E + h) * v * np.cos(gamma)                             #Downstream velocity
        dh_dt = v * np.sin(gamma)                                               #Vertical velocity component
        dg_losses_dt = g(h) * np.sin(gamma)                                     #Rate of losses due to gravity


    dd_losses_dt = D(v,h, t, stages) / m              #Rate of drag losses
    drag_test.append(D(v, h, t, stages))
    return dv_dt, dgamma_dt, dx_dt, dh_dt, dd_losses_dt, dg_losses_dt,

if stages == 1:
    tf = t_bo_s1
else:
    tf = 200                                #Final time
t_span = (t0,tf)                            #Time span for RK45 function

#####RK4(5) and solution######

sol = solve_ivp(rates, t_span, y0, method='RK45', rtol=1e-8, atol=1e-10, max_step=0.1)

v_final = sol.y[0, -1]                      #Final velocity
gamma_final = sol.y[1, -1]                  #Final pitch angle
x_final = sol.y[2, -1]                      #Final downstream displacement
h_final = sol.y[3, -1]                      #Final altitude
g_losses_final = sol.y[4, -1]               #Total gravity losses
d_losses_final = sol.y[5, -1]               #Total drag losses

time = sol.t                                #Time data
velocity = sol.y[0]                         #Velocity data
gamma = sol.y[1]                            #Gamma data
x_disp = sol.y[2]                           #X_disp data
altitude = sol.y[3]                         #Altitude data
grav = sol.y[4]
drag = sol.y[5]

drag_output = [D(v, h, t, stages) for v, h, t in zip(velocity, altitude, time)]
rho_output = [rho(h) for h in altitude]
mach_output = [mach(h,v) for h, v in zip(altitude, velocity)]
mass_output = [mass(t,stages) for t in time]
thrust_output = [thrust(t, h, stages, ThrustASL_s1, ThrustVac_s1, ThrustASL_s2, ThrustVac_s2, t_bo_s1, t_bo_s2)
                 for t, h in zip(time, altitude)]
gravity_output = [g(h) for h in altitude]

fig, results = plt.subplots(2, 2, figsize=(10, 8))

# Plot each graph in its own subplot
results[0, 0].plot(time, thrust_output, label="Thrust")
results[0, 1].plot(time, drag_output, label="Drag")
results[1, 0].plot(time, mass_output, label="Mass")
results[1, 1].plot(time, gravity_output, label="Gravity")

# Add titles
titles = ["Thrust (kN)", "Drag (kN)", "Mass (kg)", "Gravity Acc (km/s^2)"]
ylims = [[0,9], [0,3],[0,400], [0,0.01]]
for ax, title, ylim in zip(results.flat, titles, ylims):
    ax.set_title(title)
    ax.set_ylim(ylim)
    ax.set_xlim([1.5, 200])

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












