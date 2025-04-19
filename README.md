Flight Path Estimator - 2D Rocket Trajectory Simulator

Physics based vertical launch vehicle trajectory simulator. Estimates flight path using RK45 integration of Newtonian equations of motion accounting for:

- variable drag
- thrust interpolation with altitude
- gravitational acceleration and density variability with altitude
- multi-staging.

Validated against Aerobee 150A flight data (Vought Astronautics, 1961) to within ~1.5% of altitude and ~5% of range (expected due to simplified 2D assumptions and no aerodynamic instability modeling). 

Outputs:

Final / max velocity
Final / max altitude 
Range
Flight path angle
Gravity losses
Drag losses

Dependencies:

scipy
numpy
matplotlib

Example Output:
![image](https://github.com/user-attachments/assets/40ab0f88-42c9-4cc0-8d26-fce29203ed00)
"
Final velocity: 1.77 km/s
Final altitude: 9.40 km
Final Downstream Range Reached: 97.00 km
Final flight path angle: -84.8 deg
Max Altitude Reached: 285.3 km
Max Velocity Reached: 2.15 km/s
Total gravity losses: 0.35 km/s (Δv)
Total drag losses: 0.90 km/s (Δv)
"

1 - Vought Astronautics. (1961). 'Performance Summary for the Aerobee 150A Sounding Rocket'. [Dallas, Texas]: Available at: https://www.rasaero.com/dloads/Aerobee%20150A%20-%20Vought%20Astronautics%20Report%20AST%20E1R-13319.pdf. (Accessed: 19/04/2025).
