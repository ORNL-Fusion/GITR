import gitr
import shutil
import solps
import hpic
import gitrParticleSource

gitr.make_gitr_geometry_from_solps(gitr_geometry_filename='gitr_geometry.cfg', \
                                    solps_mesh_extra='mesh.extra.iter', \
                                    solps_geom = 'b2fgmtry')

shutil.copyfile('gitr_geometry.cfg', '../input/gitr_geometry.cfg')

solps.process_solps_output_for_gitr(dakota_filename = 'dakota', \
                                   nR = 500, nZ = 1000, plot_variables=1, \
                                   b2fstate_filename = 'b2fstate')

shutil.copyfile('profiles.nc', '../input/profiles.nc')

solps.readEquilibrium(filename='Baseline2008-li0.70.x4.equ')
shutil.copyfile('bField.nc', '../input/bField.nc')


solps.make_solps_targ_file(gitr_geom_filename='gitr_geometry.cfg', \
     solps_geom = 'b2fgmtry', \
     right_target_filename= 'rightTargOutput')

solps.make_solps_targ_file_txt(solps_geom='b2fgmtry',b_field_file = 'Baseline2008-li0.70.x4.equ')

hpic.plot_hpic_iead(solps_path='solpsTarg.txt', \
        HpicDataFolder = ['hPIC_IEAD_solps_conditions/hPIC_IEAD_DATA', \
        'hPIC_IEAD_He2_final/hPIC_IEAD_He2'],\
        solpsIndex = [3,4])
hpic.computeSputtYld(plots=False,hpic_zipfile='hpic_ieads',solps_inds = [3,4],nEdistPts=40)

gitrParticleSource.particleSource2d(nParticles = int(1e3),spylFile = 'Yields.txt', coordsFile='right_target_coordinates.txt',edist_file = 'Edists.txt',adist_file = 'Adists.txt')
shutil.copyfile('particleSource.nc', '../input/particleSource.nc')
