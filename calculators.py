import numpy as np
from numpy import linalg as la
from scipy.optimize import curve_fit
from scipy.spatial.distance import euclidean

A_m=10e10
h_bar=1.0545718e-34
eV=1.6e-19
me=9.10938356e-31
'''
unit conversion relationship:
E(eV)*1.6e-19=(1.0545e-34)^2 * (10^10)^2 * k2=^2 / 2m
where the k is in unit of 1/A
mass in unit of kg
'''


def fitting_functions(data,a,b,c,d,e,f,g,h,i,k):
    x, y, z = data
    return a * x ** 2 + b * y ** 2 + c * z ** 2 + d * x * y + e * x * z + f * y * z + g * x + h * y + i * z + k
    
def effectiveMassCalcu_CF(e_grid,kx_grid,ky_grid,kz_grid):
    # energy in unit of eV
    # kz in unit of 1/A
    e_fit=e_grid.flatten()
    x_fit=kx_grid.flatten()
    y_fit=ky_grid.flatten()
    z_fit=kz_grid.flatten()
    data_fit=np.vstack((x_fit,y_fit,z_fit))
    result,varance=curve_fit(fitting_functions,data_fit,e_fit,)
    second_order=np.array([[result[0],result[3],result[4]],
                           [result[3],result[1],result[5]],
                           [result[4],result[5],result[2]]])
    e_mass_tensor=1/ (second_order*(2*eV)/((h_bar*A_m)**2))
    return e_mass_tensor/me

def effectiveMassCalcu_FD(e_grid,kx_grid,ky_grid,kz_grid):
    second_order=np.zeros((3,3))
    #find the finite difference
    k_m=np.array([kx_grid[1,1,1],ky_grid[1,1,1],kz_grid[1,1,1]])
    k_x=np.array([kx_grid[2,1,1],ky_grid[2,1,1],kz_grid[2,1,1]])
    k_y=np.array([kx_grid[1,2,1],ky_grid[1,2,1],kz_grid[1,2,1]])
    k_z=np.array([kx_grid[1,1,2],ky_grid[1,1,2],kz_grid[1,1,2]])
    d=[euclidean(k_m,k_x),euclidean(k_m,k_y),euclidean(k_m,k_z)]

    for i in range(3):
        for j in range(3):
            tmp=1.0/(d[i]*d[j])
            if i==j:
                o = [1, 1, 1]
                p1 = [1, 1, 1]
                p2 = [1, 1, 1]
                p1[i]+=1
                p2[i]+=-1
                deff=tmp*(e_grid[p1[0],p1[1],p1[2]]+e_grid[p2[0],p2[1],p2[2]]-2*e_grid[o[0],o[1],o[2]])
            else:
                p1 = [1, 1, 1]
                p2 = [1, 1, 1]
                p3 = [1, 1, 1]
                p4 = [1, 1, 1]
                p1[i]+=1
                p1[j]+=1
                p2[i] += 1
                p2[j] += -1
                p3[i] += -1
                p3[j] += 1
                p4[i] += -1
                p4[j] += -1
                deff=tmp*(e_grid[p1[0]][p1[1]][p1[2]]-e_grid[p2[0]][p2[1]][p2[2]]-e_grid[p3[0]][p3[1]][p3[2]]+e_grid[p4[0]][p4[1]][p4[2]])
            second_order[i][j]=deff
    e_mass_tensor=1/ (second_order*(2*eV)/((h_bar*A_m)**2))
    return e_mass_tensor/me