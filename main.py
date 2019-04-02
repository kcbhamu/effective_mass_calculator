from input_converter.qe_converter import read_out_put
from functions import locate_extreme
from calculators import effectiveMassCalcu_CF,effectiveMassCalcu_FD
from numpy.linalg import eig

def main():
    mode="VBM"
    
    e_grid,kx_grid,ky_grid,kz_grid=read_out_put("./testdata/data-file-schema.xml",mode)

    local_extremes,extreme_positions=locate_extreme(e_grid,kx_grid,ky_grid,kz_grid,mode,energy_window=0.2)

    for i in range(len(local_extremes)):
        pos=extreme_positions[i]
        print("local extreme "+str(i)+':')
        print("positions: ")
        print((pos[0][1,1,1],pos[1][1,1,1],pos[2][1,1,1]))
        print("effective mass:")
        result=effectiveMassCalcu_FD(local_extremes[i],pos[0],pos[1],pos[2])
        value,vector=eig(result)
        print(value)

main()