# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt, animation
from mpl_toolkits.mplot3d import Axes3D


AMB_TEMP = 298
INIT_COFFEE = 361
TIME_MAX = 3600
X_MAX = 100
Y_MAX = 100
Z_MAX = 100
#ITERATIONS = 30
ITERATIONS_X = 200
ITERATIONS_Y = 200
ITERATIONS_Z = 200
STEPS_T = 50000
X_STEP = X_MAX/ITERATIONS_X
Y_STEP = Y_MAX/ITERATIONS_Y
Z_STEP = Z_MAX/ITERATIONS_Z
TIME_STEP = TIME_MAX/STEPS_T
ALPHA = 0.143 #mm units for thermal diffusivity
R_X = ALPHA * TIME_STEP / (X_STEP)**2
R_Y = ALPHA * TIME_STEP / (Y_STEP)**2
R_Z = ALPHA * TIME_STEP / (Z_STEP)**2
BETA_1 = (40/0.6)*10**-3 #heat transfer coefficient divided thermal conductivity (mm units)
BETA_2 = BETA_1 * 0.94
BETA_3 = BETA_1 * 0.94
FRAMES = int(STEPS_T/30)
Z_INDEX = 3

data = np.genfromtxt("Data.csv", delimiter=",", skip_header=0, dtype=np.float16)

stability_condition = R_X
print("Stability condition (should be <= 0.5):", R_X, R_Y, R_Z)

def convective_boundary_cond(T_far, T_close, step, beta):

    return (T_far - 2*step*beta*(T_close - AMB_TEMP))


def heat_equation(T_mid, T_left_x, T_right_x, T_left_y, T_right_y, T_left_z,
                  T_right_z):
    
    val = (T_mid + R_X*(T_right_x - 2*T_mid + T_left_x) +
            R_Y*(T_right_y - 2*T_mid + T_left_y) + 
            R_Z*(T_right_z - 2*T_mid + T_left_z))
    
    
    
    return val


def main():
    #x = np.linspace(0, X_MAX, ITERATIONS + 1)
    #y = np.linspace(0, Y_MAX, ITERATIONS + 1)
    #z = np.linspace(0, Z_MAX, ITERATIONS + 1)
    #time = np.linspace(0, TIME_MAX, STEPS_T + 1)
    stencil_current = np.full((ITERATIONS_X + 1, ITERATIONS_Y + 1, ITERATIONS_Z + 1), np.nan)
    stencil_new = np.full((ITERATIONS_X + 1, ITERATIONS_Y + 1, ITERATIONS_Z + 1), np.nan)
    temp_av = np.zeros(STEPS_T )
    #sum_ = 0
    #array_size_full = (ITERATIONS_X + 1)*(ITERATIONS_Y + 1)*(ITERATIONS_Z + 1)
    #array_size_part = (array_size_full - (ITERATIONS_X + 1)*4
                       #- (ITERATIONS_Z + 1 -2)*4 - (ITERATIONS_Y + 1 - 2)*4 - (ITERATIONS_X + 1)*(ITERATIONS_Y + 1)*2 - (ITERATIONS_X+1)*(ITERATIONS_Z+1)*2 - (ITERATIONS_Y+1)*(ITERATIONS_Z+1)*2)
    #print(array_size_part)
    #print(len(stencil[0, 0, :, 0]))

    #X, Y, Z = np.meshgrid(x, y, z, indexing="xy")
    stencil_current[:, :, :] = INIT_COFFEE

    
    
    stencil_current[0, 1:-1, 1:-1] = convective_boundary_cond(stencil_current[2, 1:-1, 1:-1] , stencil_current[1, 1:-1, 1:-1], X_STEP, BETA_2)
    stencil_current[-1, 1:-1, 1:-1] = convective_boundary_cond(stencil_current[-3, 1:-1, 1:-1] , stencil_current[-2, 1:-1, 1:-1], X_STEP, BETA_2)
    stencil_current[1:-1, 0, 1:-1] = convective_boundary_cond(stencil_current[1:-1, 2, 1:-1], stencil_current[1:-1, 1, 1:-1], Y_STEP, BETA_2)
    stencil_current[1:-1, -1, 1:-1] = convective_boundary_cond(stencil_current[1:-1, -3, 1:-1], stencil_current[1:-1, -2, 1:-1], Y_STEP, BETA_2)
    stencil_current[1:-1, 1:-1, 0] = convective_boundary_cond(stencil_current[1:-1, 1:-1, 2], stencil_current[1:-1, 1:-1, 1], Z_STEP, BETA_3)
    stencil_current[1:-1, 1:-1, -1] = convective_boundary_cond(stencil_current[1:-1, 1:-1, -3], stencil_current[1:-1, 1:-1, -2], Z_STEP, BETA_1)
        

    
    temp_av[0] = np.nanmean(stencil_current[1:-1, 1:-1, 1:-1])

    
    for i in range(1, STEPS_T):
        


        stencil_new[1:-1, 1:-1, 1:-1] = heat_equation(stencil_current[1:-1, 1:-1, 1:-1],
                                                     stencil_current[:-2, 1:-1, 1:-1],
                                                     stencil_current[2:, 1:-1, 1:-1],
                                                     stencil_current[1:-1, :-2, 1:-1],
                                                     stencil_current[1:-1, 2:, 1:-1],
                                                     stencil_current[1:-1, 1:-1, :-2],
                                                     stencil_current[1:-1, 1:-1, 2:])
        
        stencil_new[0, 1:-1, 1:-1] = convective_boundary_cond(stencil_new[2, 1:-1, 1:-1] , stencil_new[1, 1:-1, 1:-1], X_STEP, BETA_2)
        stencil_new[-1, 1:-1, 1:-1] = convective_boundary_cond(stencil_new[-3, 1:-1, 1:-1] , stencil_new[-2, 1:-1, 1:-1], X_STEP, BETA_2)
        stencil_new[1:-1, 0, 1:-1] = convective_boundary_cond(stencil_new[1:-1, 2, 1:-1], stencil_new[1:-1, 1, 1:-1], Y_STEP, BETA_2)
        stencil_new[1:-1, -1, 1:-1] = convective_boundary_cond(stencil_new[1:-1, -3, 1:-1], stencil_new[1:-1, -2, 1:-1], Y_STEP, BETA_2)
        stencil_new[1:-1, 1:-1, 0] = convective_boundary_cond(stencil_new[1:-1, 1:-1, 2], stencil_new[1:-1, 1:-1, 1], Z_STEP, BETA_3)
        stencil_new[1:-1, 1:-1, -1] = convective_boundary_cond(stencil_new[1:-1, 1:-1, -3], stencil_new[1:-1, 1:-1, -2], Z_STEP, BETA_1)

        temp_av[i] = np.nanmean(stencil_new[1:-1, 1:-1, 1:-1])
        stencil_current = stencil_new
        stencil_new = np.full((ITERATIONS_X + 1, ITERATIONS_Y + 1, ITERATIONS_Z + 1), np.nan)

    plt.plot(np.arange(STEPS_T)*TIME_STEP/60, temp_av, label="heat equation")
    plt.scatter(data[:,0], data[:,1], label="data", marker=".")
    plt.xlabel("time")
    plt.ylabel("temperature")
    plt.legend()

    plt.show()
    
    return temp_av

x = main()