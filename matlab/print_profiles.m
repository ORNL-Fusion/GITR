dlmwrite('Bfield_x_out.txt',Bfield3D.x,'delimiter','\t','precision',3)
dlmwrite('Bfield_y_out.txt',Bfield3D.y,'delimiter','\t','precision',3)
dlmwrite('Bfield_z_out.txt',Bfield3D.z,'delimiter','\t','precision',3)
dlmwrite('Bfield_mag_out.txt',Bfield3D.mag,'delimiter','\t','precision',3)

dlmwrite('Efield_x_out.txt',Efield3D.x,'delimiter','\t','precision',3)
dlmwrite('Efield_y_out.txt',Efield3D.y,'delimiter','\t','precision',3)
dlmwrite('Efield_z_out.txt',Efield3D.z,'delimiter','\t','precision',3)


dlmwrite('flowVelocity_x_out.txt',flowVelocity_ms.x,'delimiter','\t','precision',3)
dlmwrite('flowVelocity_y_out.txt',flowVelocity_ms.y,'delimiter','\t','precision',3)
dlmwrite('flowVelocity_z_out.txt',flowVelocity_ms.z,'delimiter','\t','precision',3)

dlmwrite('perDiffCoeff_out.txt',perDiffusionCoeff,'delimiter','\t','precision',3)
dlmwrite('density_m3_out.txt',density_m3,'delimiter','\t','precision',3)
dlmwrite('temp_eV_out.txt',temp_eV,'delimiter','\t','precision',3)



%