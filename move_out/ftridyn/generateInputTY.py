#!/usr/bin/python
import numpy as np
import subprocess as subp
import getopt, sys, glob, shutil, os
import FGen2
#sys.path.append("/home/kyle/TRIM_CPMI/trunk/utils/FTRIM/")

# The input File for Tridyn has the following structure:

#NAMExxxx.IN
#n
# IFOUT
# NH   IDOUT   IQOUT   NCP   IDREL   IQ0  IRC0   IRAND  JCP1  JCP2  JFRP  JNRM
# FLC    INEL    IWC    IDIFF	CSBE	ANGMN	ANGMX
# TT      TTDYN      NQX  DSF        IQXN  IQXX    IMCP
# ZZ(1)     M(1)     BE(1)     ED(1)    EF(1)     QU(1)     DNS0(1) CK(1)
# E0(1)  ALPHA0(1)  QUBEAM(1)  QUMAX(1)  SBV(1,1) ... (1,NCP)
# ZZ(2)     M(2)     BE(2)     ED(2)    EF(2)     QU(2)     DNS0(2) CK(2)
# E0(2)  ALPHA0(2)  QUBEAM(2)  QUMAX(2)  SBV(2,1) ... (2,NCP)
# ZZ(3)...
# ZZ(NCP)...

# Example: Helium implanting on Berillium
#HeBexxxx.IN
#n
#      5000
#    100000  20000  20000   2   0   0  -1   127484897  0  1  1  1
# 0.100E+01    1    1    0	0	0	90
# 0.100E+05 0.100E+05  250  100.        0  250    0
# 2.00     4.00     0.00     20.0    0.100     0.00     0.189E-01 1.0
#  30.  0.  1.  0.  0.  0.
# 4.00     9.01     0.00     20.0    0.100     1.00     0.120     1.0
#  0.  0.  0.  1.  3.38000011  3.38000011

# The input File for FTridyn (Fractal Tridyn) will have the following structure

#NAMExxxx.IN
# IFOUT
# NH   IDOUT   IQOUT   NCP   IDREL   IQ0  IRC0   IRAND  JCP1  JCP2  JFRP  JNRM
# FLC    INEL    IWC    IDIFF	CSBE	ANGMN	ANGMX
# TT      TTDYN      NQX  DSF        IQXN  IQXX    IMCP
# surfaceFile1.surf surfaceFile2.surf
# sbeadj spread
# ZZ(1)     M(1)     BE(1)     ED(1)    EF(1)     QU(1)     DNS0(1) CK(1)
# E0(1)  ALPHA0(1)  QUBEAM(1)  QUMAX(1)  SBV(1,1) ... (1,NCP)
# ZZ(2)     M(2)     BE(2)     ED(2)    EF(2)     QU(2)     DNS0(2) CK(2)
# E0(2)  ALPHA0(2)  QUBEAM(2)  QUMAX(2)  SBV(2,1) ... (2,NCP)
# ZZ(3)...
# ZZ(NCP)...

# define energy, angle, dimenstion arrays...


def main(R,s,E,a,d):

  IS_ROUGH = True
#  IS_ROUGH = False
  is_rough_inline = not IS_ROUGH
  ANALYZE_SURF = False
#  ANALYZE_SURF = False
  analyze_surf_inline = not ANALYZE_SURF

  #try:
  #  options, arguments = getopt.getopt( sys.argv[1:], 'hs:R:E:a:d:')
  #  print(options)
  #  print(arguments)
  #except getopt.GetoptError:
  #  print 'Usage:   generateInput.py [-h] -R isRoughBooleanInteger -s analyzeSurf'
  #  sys.exit(2)
  #for opt, arg in options:  #options is a list of tuples (keywords with -/--)
  #  if opt == '-h':
  #    print 'Usage:   generateInput.py [-h] -R isRoughBooleanInteger -s analyzeSurf'
  #  if opt == '-R':
  #    is_rough_inline = int(arg)
  is_rough_inline = R
  #  if opt == '-s':
  #    analyze_surf_inline = int(arg)
  analyze_surf_inline = s
  #  if opt == '-E':
  #    print(arg)
  #    E_input = float(arg)
  E_input = E
  #  if opt == '-a':
  #    print(arg)
  #    angle_input = float(arg)
  angle_input = a
  #  if opt == '-d':
  #    print(arg)
  #    dimension_input = float(arg)
  dimension_input = d

  if(IS_ROUGH!=is_rough_inline):
    print 'Command line choice and scipt choice for surface roughness are not the same!'
    print 'Aborting...'
    sys.exit(2)

  if(ANALYZE_SURF!=analyze_surf_inline):
    print 'Command line choice and scipt choice for analyzing surface are not the same!'
    print 'Aborting...'
    sys.exit(2)

  (iter_char,iter_var) = pick_iteration(dimension_input)

  useDim = True
  spacing = 2.0
  print 'iter_var size = ',np.size(iter_var,0)
  for ii in range(0,np.size(iter_var,0)):
    NAME  = 'He_W';           # Four-chars tag name
  #  SURF1 = '200.surf1';      # file name for surface 1
  #  SURF2 = '200.surf2';      # file name for surface 2
    if(not IS_ROUGH):         # after IQOUT projectile histories each.
      Y_N = 'n'               # Dynamic screen display of profiles
                              #   only needed for original tridyn (not rough)
    NH    = 1E4;              # Number of projectile histories
    IFOUT = NH/20;            # Projectile history performance information on screen
                              # after IFOUT histories each
    IDOUT = NH;             # IDOUT>0  : integral data output after each IDOUT'th history
                              # IDOUT<=0 : output suppressed
    IQOUT = NH/5;             # Additional profile output after each IQOUT'th history
    NCP   = 2;                # Number of target components (including projectile)
    IDREL = 1;                # IDREL>0 : Suppression of dynamic relaxation
                              # IDREL<0 : Suppression of dynamic relaxation and cascades
    IQ0   = 0;                # IQ0<>0  : Initial composition variable according to layer input
                              # IQ0=0   : Initial composition homogeneous
    IRC0  = -1;               # IRC0<0  : Subthreshold recoil atoms free
                              # IRC0>=0 : Subthreshold recoil atoms bound
    IRAND = 12855897;        # IRAND Initial Random Number
    JSP1  = 0;                # JSP1,JSP2 Suppression of recoils of type JSP1,...,JSP2
    JSP2  = 1;                # JSP1=0  All recoils traced
    JFRP  = 1;                # Frenkel pair generation for components JFRP,...,NCP
                              # (default: 1) (for print output (NAMEOUT.DAT) only)
    JNRM  = 1;                # JNRM : Normalization of print output (NAMEOUT.DAT) to partial
                              #        density of components JNRM,...,NCP (default: 1)


    FLC   = 1.0e-16;           #   FLC       implanted fluence [1e16/cm2)] for complete run
    INEL  = 1;                #   INEL =1   inelastic interaction nonlocal
                              #         2   local
                              #         3   equipart.
    IWC   = 3;                #   IWC       max. order of weak projectile-target collisions
    IDIFF = 1;                #   IDIFF     treatment of excess atoms (see QUMAX below)
                              #        =1   reemission
                              #        =0   "diffusion" to nearest unsaturated depth interval,
                              #	   	        reemission if whole bulk is saturated
    CSBE  = 0;                #   CSBE      surface binding energy model for compounds
                              #        =0   continuous variation (matrix model)
                              #        =x   discrete variation for a compound CPT2-CPT3(x)
                              #             in substoichiometry
    ANGMN = 0;                #   ANGMN,    list mode output of sputtered energy/angle distributions for
    ANGMX = 90;               #   ANGMX     polar emission angles between ANGMN and ANGMAX (degree)

    if(IS_ROUGH):
      spread = 1E6            #   spread    if > nset(=xsrf.size), then incident particles spread across entire target surface
      sbeadj = 0.0
      track_list = ''
      track_all    = 0        # track / create dump files for all particles?
      track_list = track_list + track_all*' ALL'
      track_target = 1        # track / create dump files for target particles?
      track_list = track_list + track_all*' TARGET'
      track_pska   = 0        # track / create dump files for pska particles?
      track_list = track_list + track_all*' PSKA'
      track_pka    = 0        # track / create dump files for pka particles?
      track_list = track_list + track_all*' PKA'
      track_ska    = 0
      track_list = track_list + track_all*' SKA'
      track_inc    = 0
      track_list = track_list + track_all*' INC'
      track_6bit   = track_inc*(2**5)+track_pka*(2**4)+track_ska*(2**3)+track_pska*(2**2)+track_target*(2**1)+track_all*(2**0)

      if(ANALYZE_SURF and (track_6bit<1) ):
        track_6bit = 2**0
        track_list = 'ALL'

    #end if(is_rough)
    TT    = 200;              #   TT        target thickness [A]
    TTDYN = 10.0*TT;           #   TTDYN     depth range [A] for dynamic relaxation
    NQX   = 500;              #   NQX       nr. of depth intervals within TTDYN
    DSF   = 100.0;            #   DSF       averaging depth [A] for surface composition
    IQXN  = 0;                #   IQXN,IQXX limiting depth intervals for profile output
    IQXX  = 250;              #             (IQXN = 0: all intervals)
    IMCP  = 0;                #   IMCP      component for which moments shall be calculated
                              #             (= 0: no moment calculation) (not used for present version)

    # Species
    ZZ    = [ 74.00, 74.00];#,1.00];    #   ZZ        Atomic number
    M     = [ 183.84, 183.84 ];#i4.002603,2.014];    #   M         Atomic mass
    BE    = [ 0.00, 0.00];#,0.00];    #   BE        Bulk binding energy [eV]
    ED    = [ 20.0, 20.0];#,20.0];    #   ED        Displacement energy [eV]
    EF    = [ 0.10, 0.10];#,0.10];    #   EF        Cutoff energy [eV]
    QU    = [ 0.00, 1.00];#,0.40];    #   QU        Initial atomic fraction [QU(1) ignored, obtained from normaliz.]
    DNS0  = [ 0.06306, 0.06306];#0.01878,0.0427];  #   DNS0      Atomic density of pure component [10^24 atoms/cm3]
    CK    = [ 1.00, 1.00];#,1.00];    #   CK        Electronic stopping correction factor

    E0     = [ E_input, 0.00];   #   E0        Incident energy [eV]; If <0 than en. is tabulated from input file
    ALPHA0 = [ angle_input, 0.00];#,0.00];   #   ALPHA0    Angle of incidence [deg] with respect to surface normal
    QUBEAM = [ 1.00, 0.00];#,0.00];   #   QUBEAM    Fraction of component in incident beam [QUBEAM(1) ignored]
    numBeam = QUBEAM.index(0.0);
    QUMAX  = [ 0.00, 1.00];#,1.00];   #   QUMAX     Max allowed bulk atomic fraction
    SBV    = [[8.79, 8.79], #[0.00, 0.00] for He
              [8.79, 8.79]];   #   SBV(I,J)     Surface binding energy vector (NCP components per species)
              #SBV is a list of lists, therefore individual entries are accessed with SBV['list/row i']['entry/col j']
    if(NCP!=np.size(ZZ)):
      print 'NCP != size(ZZ)'

    if(IS_ROUGH):
      surfaceType = 'c'
      dim = 2.0 - 1.0 #   DIM       Fractal Dimension. 2.0 is flat. Max 3.0

      beta = beta = 60.0*np.pi/180.
      Range = 2*np.array([TT])#np.array([2.1*TT])

    # switch for iterating variable

    if   iter_char=='NH':
      NH = iter_var[ii]
      IFOUT = NH/20
      IDOUT = NH/5
      IQOUT = NH/5
      pass
    elif iter_char=='IFOUT':
      IFOUT = iter_var[ii]
      pass
    elif iter_char=='IDOUT':
      IDOUT = iter_var[ii]
      pass
    elif iter_char=='IQOUT':
      IQOUT = iter_var[ii]
      pass
    elif iter_char=='IDREL':
      IDREL = iter_var[ii]
      pass
    elif iter_char=='IQ0':
      IQ0 = iter_var[ii]
      pass
    elif iter_char=='IRC0':
      IRC0 = iter_var[ii]
      pass
    elif iter_char=='IRAND':
      IRAND = iter_var[ii]
      pass
    elif iter_char=='FLC':
      FLC = iter_var[ii]
      pass
    elif iter_char=='INEL':
      INEL = iter_var[ii]
      pass
    elif iter_char=='IWC':
      IWC = iter_var[ii]
      pass
    elif iter_char=='IDIFF':
      IDIFF = iter_var[ii]
      pass
    elif iter_char=='TT':
      TT = iter_var[ii]
      TTDYN = 1.2*TT
      pass
    elif iter_char=='TTDYN':
      TTDYN = iter_var[ii]
      pass
    elif iter_char=='TT_TTDYN':
      TT    = TT
      TTDYN = TTDYN
    elif iter_char=='NQX':
      NQX = iter_var[ii]
      pass
    elif iter_char=='DSF':
      DSF = iter_var[ii]
      pass
    elif iter_char=='ZZ':
      ZZ = iter_var[ii]
      pass
    elif iter_char=='M':
      M = iter_var[ii]
      pass
    elif iter_char=='BE':
      BE = iter_var[ii]
      pass
    elif iter_char=='ED':
      ED = iter_var[ii]
      pass
    elif iter_char=='EF':
      EF = iter_var[ii]
      pass
    elif iter_char=='QU':
      QU = iter_var[ii]
      pass
    elif iter_char=='DNS0':
      DNS0 = iter_var[ii]
      pass
    elif iter_char=='CK':
      CK = iter_var[ii]
      pass
    elif iter_char=='E0':
      E0 = iter_var[ii]
      pass
    elif iter_char=='ALPHA0':
      ALPHA0 = iter_var[ii]
      pass
    elif iter_char=='QUBEAM':
      QUBEAM = iter_var[ii]
      numBeam = QUBEAM.index(0.0);
      pass
    elif iter_char=='QUMAX':
      QUMAX = iter_var[ii]
      pass
    elif iter_char=='SBV':
      SBV = iter_var[ii]
      pass
    elif iter_char=='dim':
      dim = iter_var[ii]
      pass
    elif iter_char=='beta':
      beta = iter_var[ii]
      pass
    elif iter_char=='surfaceType':
      surfaceType = iter_var[ii]
      pass
    elif iter_char=='sbeadj':
      sbeadj = iter_var[ii]
      pass
    elif iter_char=='spread':
      spread = iter_var[ii]
      pass
    elif iter_char=='JNRM':
      JNRM = iter_var[ii]
      pass
    elif iter_char=='surfaceType':
      surfaceType = iter_var[ii]
      pass
    #end switch

    if(IS_ROUGH):
      betaDeg = beta*180./np.pi
      if(useDim):
        beta = FGen2.FAngle(dim,surfaceType)
        beta = beta[0]
      else:
        dim = FGen2.FDimd(betaDeg,surfaceType)
#        dim = FGen2.FDim(beta,surfaceType)
      surfFileName,srfseg = write_surf_files(dim,beta,Range,spacing,surfaceType)
      if(IQ0==0):
        dns = 0.0
        quSum = 0.0
        for i in range(1, NCP):
          quSum = quSum + QU[i]
          dns = dns + QU[i]/DNS0[i]
        dns = 1./(dns + (1.-quSum)/DNS0[0])
      fscale = dns**(-1./3.)/srfseg
      #print 'orig TTDYN = ',TTDYN
      #print 'fscale     = ',fscale
      #TTDYN = TTDYN*fscale
      #print 'curr TTDYN = ',TTDYN
      #TT    = TT*fscale

    print 'E0 = ',E0
    # write input file to disk
    FileName = NAME + '%04i' % (ii+1) + '.IN' #sprintf('%03.0f',ii) '.IN'];
    fid = open(FileName, 'w')
    fid.write(FileName + '\n')
    if(not IS_ROUGH):
      fid.write(Y_N + '\n') #pplot no longer needed in FTridyn
    fid.write('%i\n' % IFOUT)
    fid.write('%i %i %i %i %i %i %i %i %i %i %i %i\n' % \
               (NH, IDOUT, IQOUT, NCP, IDREL, IQ0, IRC0, IRAND, JSP1, JSP2, JFRP, JNRM))
    fid.write('%e %i %i %i %i %i %i\n' % \
               (FLC, INEL, IWC, IDIFF, CSBE, ANGMN, ANGMX))
    fid.write('%e %e %i %e %i %i %i\n' % \
               (TT, TTDYN, NQX, DSF, IQXN, IQXX, IMCP))
    if(IS_ROUGH):
      fid.write('%9s ' % (surfFileName))#,SURF2,RMS))
      #fid.write('%10.4f %6i\n' % (sbeadj,spread))
      #fid.write('%3i %3i %3i\n' % (track_6bit,track_sp_3bit,OUT_ROUTINE))
      fid.write('1 1 1\n')
    for jj in range(0, NCP):
      fid.write('%4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e \n' % \
               (ZZ[jj], M[jj], BE[jj], ED[jj], EF[jj], QU[jj], DNS0[jj], CK[jj] ))
      fid.write('%4.3e %4.3e %4.3e %4.3e ' % \
               (E0[jj], ALPHA0[jj], QUBEAM[jj], QUMAX[jj]))
      for kk in range(0, NCP):
        fid.write('%4.3e ' % (SBV[jj][kk]))
      fid.write('\n')
    fid.close()
#    fid.write('%4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e \n' % \
#               (ZZ[0], M[0], BE[0], ED[0], EF[0], QU[0], DNS0[0], CK[0] ))
#    fid.write('%4.3e %4.3e %4.3e %4.3e %4.3e %4.3e \n' % \
#               (E0[0], ALPHA0[0], QUBEAM[0], QUMAX[0], SBV[0][0], SBV[0][1] ))#, SBV[0][2], SBV[0][3] ))
#    fid.write('%4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e %4.3e \n' % \
#               (ZZ[1], M[1], BE[1], ED[1], EF[1], QU[1], DNS0[1], CK[1] ))
#    fid.write('%4.3e %4.3e %4.3e %4.3e %4.3e %4.3e \n' % \
#               (E0[1], ALPHA0[1], QUBEAM[1], QUMAX[1], SBV[1][0], SBV[1][1] ))#, SBV[1][2], SBV[1][3] ))
#    fid.close()

    for ie in range(0,np.size(E0)):
      E = E0[ie]
      if E < 0:
        EnDist = abs( np.array( [E/6., 2*E/6., 3*E/6., 4*E/6., 5*E/6., E]     ) )
        FnDist =      np.array( [0.0 , 0.375 , 0.50  , 0.0   , 0.0   , 0.125] )
        EnDist = EnDist.reshape((EnDist.size,1))
        FnDist = FnDist.reshape((EnDist.size,1))
        EdisFileName = NAME + '%04i' % (ii+1) + '.ED%1.0i' % (ie+1)
  #      fid = open(EdistFileName, 'w')
        np.savetxt(EdisFileName, np.concatenate( (EnDist,FnDist), axis=1), fmt='%4.3f,%4.3f')
  #end for loop

  if(IS_ROUGH and ANALYZE_SURF):
    use_numBins_or_sizes = 2    # 1 to use the number of bins
                                # 2 to use the sizes of bins (default)
    if( (use_numBins_or_sizes<0) or (use_numBins_or_sizes>2) ):
      print 'PICK ONE CORRECT BIN METHOD (PITY THE FOOL!)'
      use_numBins_or_sizes = 2
    if(use_numBins_or_sizes==1):
      bin_array = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
    elif(use_numBins_or_sizes==2):
      bin_array = [.25,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,\
                          3.,3.25,3.5,3.75,4.,4.25,4.5,4.75,5.,6.,7.,8.,9.,10.0]
    bin_array = np.array( bin_array )

    FileName = NAME + '.binsurfIN'
    fid = open(FileName, 'w')
    fid.write(FileName+'\n')
    fid.write('%2i\n' % track_6bit)
    fid.write('%2i %2i\n' % (use_numBins_or_sizes,len(bin_array)) )
    for i in range( len(bin_array) ):
      fid.write('%f ' % bin_array[i] )
    #enddo

    track_list = track_list.split()
    track_list.sort()
    np.savetxt('dim.out',[],header=' '.join(track_list),comments='')

  #endif IS_ROUGH

  zz = np.array([ZZ]) #make into shape (1,n) row vector
  np.savetxt('SputYld.out',zz,fmt='%i') #save atomic numbers so that:
                                        #   -file is not empty
                                        #   -a psuedo-way to identify columns
  np.savetxt('Reflect.out',zz[0,0:numBeam],fmt='%i')
  np.savetxt('ReflectedEnergies.out',zz[0,0:numBeam],fmt='%i')

  np.savetxt('Transmitted.out',zz[0,0:numBeam],fmt='%i')
  np.savetxt('TransmittedEnergies.out',zz[0,0:numBeam],fmt='%i')

  if iter_char=='surfaceType':
    np.savetxt('Iterating_'+iter_char+'.out',iter_var[:,0], fmt='%c')
  else:
    np.savetxt('Iterating_'+iter_char+'.out',iter_var[:] )
#------------------------------ end main -------------------------------
#-----------------------------------------------------------------------



#-----------------------------------------------------------------------
#                          pick_iteration()
#-------------------- PICK VARIABLE TO ITERATE OVER --------------------
def pick_iteration(dimension_input):

  # NAME  = 'D_Be';           # Four-chars tag name
  #  Y_N   = 'n';              # Dynamic screen display of profiles
  #                            # after IQOUT projectile histories each.
  # NH    = 1E5;              # Number of projectile histories
  # iter_char = 'NH'           # iterating this may need to change those values dependent on NH
  # iter_list = 10**(np.arange(2,7,1).reshape(-1, 1))

  # IFOUT = NH/20;            # Projectile history performance information on screen                             # after IFOUT histories each
  # iter_char= 'IFOUT'
  # iter_list = np.array([NH/40, NH/20, NH/10, NH/5]).reshape(-1,1)
  # IDOUT = NH/5;             # IDOUT>0  : integral data output after each IDOUT'th history
                            # IDOUT<=0 : output suppressed
  #iter_char = 'IDOUT'
  #iter_list = np.arange(0, NH+1, NH/10).reshape(-1,1)

  # IQOUT = NH/5;             # Additional profile output after each IQOUT'th history
  # iter_char = 'IQOUT'
  # iter_list = np.arange(NH/25,NH+1, NH/25).reshape(-1,1)

  # NCP   = 2;                # Number of target components (including projectile)
  #Iterating this would only work if we increased the number of types of particles in the system
  # IDREL = 0;                # IDREL>0 : Suppression of dynamic relaxation
                            # IDREL<0 : Suppression of dynamic relaxation and cascades
  # iter_char = 'IDREL'
  # iter_list = np.array([-1, 0, 1]).reshape(-1, 1)

  # IQ0   = 0;                # IQ0=-1  : Initial composition variable according to layer input
                            # IQ0=0   : Initial composition homogeneous
#  iter_char = 'IQ0'
#  iter_list = np.array([0, -1]).reshape(-1, 1)
  # for IQ0 = -1, another file is needed (.LAY)

  # IRC0  = -1;               # IRC0<0  : Subthreshold recoil atoms free
                            # IRC0>=0 : Subthreshold recoil atoms bound
  # iter_char = 'IRC0'
  # iter_list = np.array([-1, 1]).reshape(-1,1)

  # IRAND = 127484897;        # IRAND Initial Random Number

  # JSP1  = 0;                # JSP1,JSP2 Suppression of recoils of type JSP1,...,JSP2
  # JSP2  = 1;                # JSP1=0  All recoils traced

  # iter_char = 'JSP1'
  # iter_list = np.arange(0,JSP2+1,1)
  # iter_char = 'JSP2'
  # iter_list = np.arange(1,5,1).reshape(-1,1)

  # JFRP  = 1;                # Frenkel pair generation for components JFRP,...,NCP
                            # (default: 1) (for print output (NAMEOUT.DAT) only)
  # iter_char = 'JFRP'
  # iter_list = np.arange(1, 6, 1).reshape(-1,1)
  # JNRM  = 1;                # JNRM : Normalization of print output (NAMEOUT.DAT) to partial
                            #        density of components JNRM,...,NCP (default: 1)
  # iter_char = 'JNRM'
  # iter_list = np.array([1, 2, 3, 4, 5]).reshape(-1, 1)
  # FLC   = 1.0e-9;           #   FLC       implanted fluence [1e16/cm2)] for complete run
  # INEL  = 1;                #   INEL =1   inelastic interaction nonlocal
  # iter_char = 'INEL'                          #         2   local
  # iter_list = np.array([1, 2, 3]).reshape(-1, 1)                         #         3   equipart.
  # IWC   = 3;                #   IWC       max. order of weak projectile-target collisions
  # iter_char = 'IWC'
  # iter_list = np.array([1, 2, 3]).reshape(-1,1)
  # IDIFF = 0;                #   IDIFF     treatment of excess atoms (see QUMAX below)
                            #        =1   reemission
                            #        =0   "diffusion" to nearest unsaturated depth interval,
                            #	   	        reemission if whole bulk is saturated
  # CSBE  = 0;                #   CSBE      surface binding energy model for compounds
                            #        =0   continuous variation (matrix model)
                            #        =x   discrete variation for a compound CPT2-CPT3(x)
                            #             in substoichiometry
  # ANGMN = 0;                #   ANGMN,    list mode output of sputtered energy/angle distributions for
  # ANGMX = 90;               #   ANGMX     polar emission angles between ANGMN and ANGMAX (degree)

  # spread = 2*xsrf.size      #   spread    if > nset(=xsrf.size), then incident particles spread across entire target surface
  # sbeadj = 0.0

  # TT    = Range[0];         #   TT        target thickness [A]
  # iter_char = 'TT'
  # iter_list = 10**np.arange(2, 6, 1).reshape(-1,1)

  # TTDYN = Range[0];         #   TTDYN     depth range [A] for dynamic relaxation
  # TT = 1E3
  # iter_char = 'TTDYN'
  # iter_list = TT*np.arange(0.2, 2.001, 0.2).reshape(-1,1)

  # NQX   = 500;              #   NQX       nr. of depth intervals within TTDYN
  # iter_char = 'NQX'
  # iter_list = np.arange(100, 501, 100).reshape(-1, 1)

  # DSF   = 100.0;            #   DSF       averaging depth [A] for surface composition
  # TTDYN = 1.2E3
  # iter_char = 'DSF'
  # iter_list = TTDYN*np.arange(0,1,0.1).reshape(-1,1)

  # IQXN  = 0;                #   IQXN,IQXX limiting depth intervals for profile output
  # iter_char = 'IQXN'
  # iter_list = np.arange(0, 250, 50).reshape(-1,1)

  # IQXX  = 250;              #             (IQXN = 0: all intervals)
  # iter_char = 'IQXX'
  # iter_list = np.arange(0, 500, 50).reshape(-1,1)

  # IMCP  = 0;                #   IMCP      component for which moments shall be calculated
                            #             (= 0: no moment calculation) (not used for present version)

  # Species

  # ZZ    = [ 1.00, 4.00];#,1.00];    #   ZZ        Atomic number
  # M     = [ 2.014, 9.012];#,2.014];    #   M         Atomic mass
  # BE    = [ 0.00, 0.00];#,0.00];    #   BE        Bulk binding energy [eV]
  # ED    = [ 20.0, 20.0];#,20.0];    #   ED        Displacement energy [eV]
  # EF    = [ 0.10, 0.10];#,0.10];    #   EF        Cutoff energy [eV]
  # QU    = [ 0.00, 1.00];#,0.40];    #   QU        Initial atomic fraction [QU(1) ignored, obtained from normaliz.]
  # DNS0  = [ 0.0427, 0.1235];#,0.0427];  #   DNS0      Atomic density of pure component [10^24 atoms/cm3]
  # CK    = [ 1.00, 1.00];#,1.00];    #   CK        Electronic stopping correction factor

#  E0     = [ Energy[ii], 0.00];   #   E0        Incident energy [eV]; If <0 than en. is tabulated from input file
#  iter_char = 'E0'
#  a = np.arange(100,101,5)
#  a = -1*a                    #for Energy distribution
#  a = a.reshape((a.size,1))
#  b = np.zeros((a.size,1))
#  iter_list = np.hstack((a,b))

  # ALPHA0 = [ 0.00, 0.00];#,0.00];   #   ALPHA0    Angle of incidence [deg] with respect to surface normal
  #  iter_char = 'ALPHA0'
  #  a = np.arange(0.0,86.0,5.)
  #  a = a.reshape((a.size,1))
  #  b = np.zeros((a.size,1))
  #  iter_list = np.hstack((a,b)) #angle (from 0.0 to 85.0 incrementing by 5 degrees)

  # QUBEAM = [ 1.00, 0.00];#,0.00];   #   QUBEAM    Fraction of component in incident beam [QUBEAM(1) ignored]
  # numBeam = QUBEAM.index(0.0);
  # QUMAX  = [ 0.40, 1.00];#,1.00];   #   QUMAX     Max allowed bulk atomic fraction
  # SBV    = [[2.00, 2.00],
  #           [3.38, 3.38]];   #   SBV(I,J)     Surface binding energy vector (NCP components per species)
  iter_char = 'dim'
  iter_list = np.arange(1.0,1.51,0.1).reshape(-1,1) #fractal dimension
  iter_list = np.array([dimension_input])
#  surfaceType = 'A' OR 'B' OR 'C'      # surfaceType rough surface character that defines the type of fractal on the sufrace
#  iter_char = 'surfaceType'
#  iter_list = np.array(['A','B','C']).reshape(-1,1)

  return iter_char,iter_list
#------------------------- end pick_iteration --------------------------



#-----------------------------------------------------------------------
#------------------------- WRITE SURFACE FILES -------------------------
def write_surf_files(dim,beta,Range,spacing,surfaceType):
  (xgen,ygen) = FGen2.FGen2(Range,beta,Range/2.,surfaceType)
  (xsrf,ysrf) = FGen2.FGen2(Range,beta,spacing,surfaceType)
  numGens = np.log(xsrf.size-1)/np.log(xgen.size-1)
  srfseg = np.sqrt( (xsrf[1] - xsrf[0])**2. + (ysrf[1] - ysrf[0])**2. )
  surfFileName = 'surface.surf'#'1p' + '%03i' % (np.round(1E3*(dim-1))) +'.surf'
  fid = open(surfFileName, 'w')
  fid.write(' %8i%8i%8i\n'% (xgen.size,numGens,xsrf.size) )
  fid.write(' %6.4f\n' % (dim+1) )
  for i in range(0,xgen.size):
    fid.write(' %10.4f%10.4f\n' % (xgen[i],ygen[i]) )
  for i in range(0,xsrf.size):
    fid.write(' %10.4f%10.4f\n' % (xsrf[i],ysrf[i]) )
  fid.close()

  SURF1 = '%03i' % (np.round(1E3*(dim-1))) + '.surf1'
  SURF2 = '%03i' % (np.round(1E3*(dim-1))) + '.surf2'

  #shutil.copy(surfFileName, SURF1)
  #shutil.copy(surfFileName, SURF2)

  return surfFileName, srfseg

#------------------------- end write_surf_files ------------------------
#-----------------------------------------------------------------------



#if script is called, use main ... something like that
if __name__ == "__main__":
  main()
