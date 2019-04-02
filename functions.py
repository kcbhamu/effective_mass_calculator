import numpy as np
from scipy.signal import argrelextrema

def locate_extreme(energy_grid,kx_grid,ky_grid,kz_grid,option,energy_window):
    #expand the original grid to each direction by 1 using the periodic boundary condition
    e=1 #how much to expand on each direction
    mp=energy_grid.shape
    new_mp=[mp[0]+2*e,mp[1]+2*e,mp[2]+2*e]
    expand_energy=np.zeros((new_mp[0],new_mp[1],new_mp[2]))
    expand_kx=np.zeros((new_mp[0],new_mp[1],new_mp[2]))
    expand_ky=np.zeros((new_mp[0],new_mp[1],new_mp[2]))
    expand_kz=np.zeros((new_mp[0],new_mp[1],new_mp[2]))

    for i in range(new_mp[0]):
        for j in range(new_mp[1]):
            for k in range(new_mp[2]):
                #map back to the old position
                if i<e:
                    old_x=i+mp[0]-e
                elif i>=mp[0]+e:
                    old_x=i-mp[0]-e
                else:
                    old_x=i-e

                if j<e:
                    old_y=j+mp[1]-e
                elif j>=mp[1]+e:
                    old_y=j-mp[1]-e
                else:
                    old_y=j-e

                if k<e:
                    old_z=k+mp[2]-e
                elif k>=mp[2]+e:
                    old_z=k-mp[2]-e
                else:
                    old_z=k-e

                expand_energy[i][j][k]=energy_grid[old_x][old_y][old_z]
                expand_kx[i][j][k]=kx_grid[old_x][old_y][old_z]
                expand_ky[i][j][k]=ky_grid[old_x][old_y][old_z]
                expand_kz[i][j][k]=kz_grid[old_x][old_y][old_z]

    #find the extreme
    if option=='VBM':
        is_vbm=True
    elif option=='CBM':
        is_vbm=False
    else:
        print("wrong options")
        return
    
    if is_vbm:
        maximum=energy_grid.max()
    else:
        maximum=energy_grid.min()

    matrix_extreme=[]
    position_extreme=[]

    positions=[]
    for axis_num in range(3):
        if is_vbm:
            ax1,ax2,ax3=argrelextrema(expand_energy,np.greater,axis=axis_num)
        else:
            ax1,ax2,ax3=argrelextrema(expand_energy,np.less,axis=axis_num)
        
        if axis_num==0:
            #append every position found
            for i in range(len(ax1)):
                positions.append([ax1[i],ax2[i],ax3[i]])
        else:
            #select the points that meet the criteria still
            tmp_position=[]
            for i in range(len(ax1)):
                if [ax1[i],ax2[i],ax3[i]] in positions:
                    tmp_position.append([ax1[i],ax2[i],ax3[i]])
            positions=tmp_position

    for pos in positions:
        xx=pos[0]
        yy=pos[1]
        zz=pos[2]
        if xx in range(e,mp[0]+e) and yy in range(e,mp[1]+e) and zz in range(e,mp[2]+e):
            extreme=expand_energy[xx-e:xx+1+e,yy-e:yy+1+e,zz-e:zz+1+e]
            extreme_x=expand_kx[xx-e:xx+1+e,yy-e:yy+1+e,zz-e:zz+1+e]
            extreme_y=expand_ky[xx-e:xx+1+e,yy-e:yy+1+e,zz-e:zz+1+e]
            extreme_z=expand_kz[xx-e:xx+1+e,yy-e:yy+1+e,zz-e:zz+1+e]
            if is_vbm:
                if extreme[1,1,1]==extreme.max() and extreme[1,1,1]>maximum-energy_window:
                    matrix_extreme.append(extreme)
                    position_extreme.append([extreme_x,extreme_y,extreme_z])
            else:
                if extreme[1,1,1]==extreme.min() and extreme[1,1,1]<maximum-energy_window:
                    matrix_extreme.append(extreme)
                    position_extreme.append([extreme_x,extreme_y,extreme_z])
    return matrix_extreme,position_extreme
    #the return is matrix_extreme, which is the matrix 3*3*3
    #position_extreme is a list of (array_x,array_y,array_z)

