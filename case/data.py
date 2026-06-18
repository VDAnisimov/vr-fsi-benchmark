import numpy as np
import math
from scipy.spatial.transform import Rotation as R
import re
import sys

# ============================================================
# 1. load params
# ============================================================
params_file = '0/params'
params = {}
with open(params_file, 'r') as f:
    content = f.read()
    content = re.sub(r'//.*', '', content)
    for key in ['Vel', 'gamma', 'beta', 'phi']:
        match = re.search(rf'{key}\s+([\d.-]+)', content)
        if match:
            params[key] = float(match.group(1))

required = ['gamma', 'beta', 'phi']
for r in required:
    if r not in params:
        sys.exit(f"Error: Parameter {r} not found in {params_file}")

gamma = params['gamma']
beta = params['beta']
phi = params['phi']

# ============================================================
# 2. load data
# ============================================================
data = np.loadtxt('./t_vs_cm')
data2 = np.loadtxt('./t_vs_orientation')
data_vel_lin = np.loadtxt('./t_vs_lv')
data_vel_ang = np.loadtxt('./t_vs_av')

time_data = data[1:-1, 0]
dy = data[1:-1, 2]
rot_raw = data2[1:-1, 1:10]
angles = np.zeros((len(rot_raw), 3))
for i in range(len(rot_raw)):
    r = R.from_matrix(rot_raw[i].reshape(3, 3))
    angles[i] = r.as_euler('zxy')
rot_signal = -angles[:, 0]

body_vel = data_vel_lin[:, 2]
body_ang_vel = data_vel_ang[:, 3]

min_len = min(len(time_data), len(body_vel), len(body_ang_vel))
time_data = time_data[:min_len]
dy = dy[:min_len]
rot_signal = rot_signal[:min_len]
body_vel = body_vel[:min_len]
body_ang_vel = body_ang_vel[:min_len]

# ============================================================
# 3. kappa and theta
# ============================================================
dt = time_data[1] - time_data[0]
N = int(4.0 / dt)
N = min(N, len(dy))

kappa = (np.max(dy[-N:]) - np.min(dy[-N:])) / 2.0
theta = (np.max(rot_signal[-N:]) - np.min(rot_signal[-N:])) / 2.0

# ============================================================
# 4. u_st
# ============================================================
coef_file = 'postProcessing/BCcontrol/0/coefficient.dat'
try:
    coef_data = np.loadtxt(coef_file)
    if coef_data.ndim == 1:
        u_raw = coef_data[2] if len(coef_data) > 2 else 0.0
    else:
        u_raw = coef_data[-1, 2]
    u_st = -u_raw * gamma / kappa if kappa != 0 else 0.0
except:
    u_st = 0.0

# ============================================================
# 5. phase
# ============================================================
T = 2.0 * np.pi * gamma
omega = 1.0 / gamma
IM_vel = -phi * omega * np.sin(omega * time_data)

def find_zeros(signal):
    zeros = []
    for i in range(len(signal)-1):
        if signal[i] < 0 and signal[i+1] > 0:
            zeros.append(i+1)
    return np.array(zeros)

last_period_start = time_data[-1] - T
last_idx = np.argmin(np.abs(time_data - last_period_start))
if last_idx < 0:
    last_idx = 0

time_last = time_data[last_idx:]
IM_last = IM_vel[last_idx:]
body_last = body_vel[last_idx:]
ang_last = body_ang_vel[last_idx:]

zeros_IM = find_zeros(IM_last)
zeros_body = find_zeros(body_last)
zeros_ang = find_zeros(ang_last)

phase_linear_deg = 0.0
phase_angular_deg = 0.0

if len(zeros_IM) > 0 and len(zeros_body) > 0:
    t_IM = time_last[zeros_IM[0]]
    t_body = time_last[zeros_body[0]]
    phase_linear_raw = 360 * (t_body - t_IM) / T
    phase_linear_deg = (phase_linear_raw + 180) % 360 - 180

if len(zeros_IM) > 0 and len(zeros_ang) > 0:
    t_IM = time_last[zeros_IM[0]]
    t_ang = time_last[zeros_ang[0]]
    phase_angular_raw = 360 * (t_ang - t_IM) / T
    phase_angular_deg = (phase_angular_raw + 180) % 360 - 180

phi_V_norm = phase_linear_deg / 180.0
phi_Omega_norm = phase_angular_deg / 180.0

# ============================================================
# 6. result
# ============================================================
print("\n" + "="*110)
print("SUMMARY TABLE")
print("="*110)
header = f"{'beta':>10} {'gamma':>10} {'phi':>10} {'u_st':>14} {'kappa':>14} {'theta':>14} {'phi_V/pi':>14} {'phi_Omega/pi':>14}"
print(header)
print("-"*110)
row = f"{beta:10.4f} {gamma:10.4f} {phi:10.4f} {u_st:14.4f} {kappa:14.4f} {theta:14.4f} {phi_V_norm:14.4f} {phi_Omega_norm:14.4f}"
print(row)
print("="*110)

with open('results_table.txt', 'w') as f:
    f.write("beta gamma phi u_st kappa theta phi_V/pi phi_Omega/pi\n")
    f.write(f"{beta:.4f} {gamma:.4f} {phi:.4f} {u_st:.4f} {kappa:.4f} {theta:.4f} {phi_V_norm:.4f} {phi_Omega_norm:.4f}\n")
