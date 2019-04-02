"""
this file converts the Quantum espresso output and sort the data out into the input for the Effective mass calculator.
In particule, the file it reads is the XML file after the dos calculation. 

The process does:
read the result of dos calculation performed on a relatively dense k grid.
for each k point in reciprocal space, it finds the highest energy below the Ef or the lowest energy above Ef
(eg: in the case of 2D k space, the result would be the surface of band maximum)
"""

import xml.etree.ElementTree as ET
import numpy as np


Ry = 13.6057
Bohr = 0.52917721

def read_out_put(xml_file_name,option):
    #band structure node is the xml node that is found by root.find('output').find('bandstructure')
    #option is either VBM or CBM
    #the return is a list of k in unit (2Pi/Bohr) and energy in eV
    if option=='VBM':
        is_next_one=False
    elif option=='CBM':
        is_next_one=True
    else:
        print("wrong options")

    output=ET.parse(xml_file_name)

    bandstructure=output.find('output').find('band_structure')

    e_fermi=float(bandstructure.find('fermi_energy').text) * 2 * Ry
    #fermi energy in eV

    alat=float(output.find('output').find('atomic_structure').get('alat')) * Bohr
    twopi_over_alat=np.pi*2/alat
    #alat is the lattice parameter, which is coverted to a.
    #two_pi_over_alat is the 2pi/alat.  the unit for recipicroal lattice

    reciprocal_vector=output.find('output').find('basis_set').find('reciprocal_lattice')
    tmp=reciprocal_vector.find('b1').text.split()
    b1=np.array([float(tmp[0])*twopi_over_alat,float(tmp[1])*twopi_over_alat,float(tmp[2])*twopi_over_alat])
    tmp=reciprocal_vector.find('b2').text.split()
    b2=np.array([float(tmp[0])*twopi_over_alat,float(tmp[1])*twopi_over_alat,float(tmp[2])*twopi_over_alat])
    tmp=reciprocal_vector.find('b3').text.split()
    b3=np.array([float(tmp[0])*twopi_over_alat,float(tmp[1])*twopi_over_alat,float(tmp[2])*twopi_over_alat])
    
    print("check reciprocal lattice vector(in unit 1/A):")
    print("b1= "+str(round(b1[0],8))+', '+str(round(b1[1],8))+', '+str(round(b1[2],8)))
    print("b2= "+str(round(b2[0],8))+', '+str(round(b2[1],8))+', '+str(round(b2[2],8)))
    print("b3= "+str(round(b3[0],8))+', '+str(round(b3[1],8))+', '+str(round(b3[2],8)))
    print('')

    tmp=bandstructure.find('starting_k_points').find('monkhorst_pack')
    mp=(int(tmp.get('nk1')),int(tmp.get('nk2')),int(tmp.get('nk3')))
    print("monkhorst_pack grid used is:"+' '+str(mp[0])+' '+str(mp[1])+' '+str(mp[2])+'\n')

    #the properties read by the above part: 
    #   e_fermi (eV)
    #   alat    (A)
    #   twopi_over_alat (1/A)
    #   mp (three integral)
    #   b1,b2,b3 (1/A)
    
    position_kx=np.zeros(mp)
    position_ky=np.zeros(mp)
    position_kz=np.zeros(mp)
    k_energy=np.zeros(mp)

    i=0
    j=0
    k=0

    reduced_fermi=e_fermi/(2*Ry)

    for kpoints in bandstructure.iter('ks_energies'):
        kxyz=kpoints.find('k_point').text.split()
        kx=float(kxyz[0])*twopi_over_alat
        ky=float(kxyz[1])*twopi_over_alat
        kz=float(kxyz[2])*twopi_over_alat
        #kx,ky,kz is in unit of 1/A
        energys=kpoints.find('eigenvalues').text.split()
        i_edge=0
        while i_edge<=len(energys):
            if float(energys[i_edge])>reduced_fermi:
                if is_next_one:
                    e_band=float(energys[i])*2*Ry
                else:
                    e_band=float(energys[i-1])*2*Ry
                break
            else:
                i_edge+=1
        #E_band is in unit of eV
        position_kx[i,j,k]=kx
        position_ky[i,j,k]=ky
        position_kz[i,j,k]=kz
        k_energy[i,j,k]=e_band

        if k==mp[2]-1:
            k=0
            if j==mp[1]-1:
                j=0
                i+=1
            else:
                j+=1
        else:
            k+=1
    # the E_grid and position_grid is created
    
    
    if is_next_one:
        print("the conduction band minimum at: " + str(k_energy.max()))
    else:
        print("the valence band maximum at: " + str(k_energy.max()))
    
    return k_energy,position_kx,position_ky,position_kz
    


