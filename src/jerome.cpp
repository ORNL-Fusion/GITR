/*
 * jerome.cpp
 *
 *  Created on: Dec 15, 2021
 *      Author: guterl
 */


 /* Start BLOCK: write_connection_length(...) - easy - netcdf  */
void write_lcs(int nTracers, int nR_Lc, int nY_Lc, int nZ_Lc, sim::Array<gitr_precision> &gridZLc, sim::Array<gitr_precision> &gridYLc, sim::Array<gitr_precision> &gridZLc)
{
    NcFile ncFileLC("LcS.nc", NcFile::replace);
    vector<NcDim> dims_lc;
    NcDim nc_nTracers = ncFileLC.addDim("nTracers", nTracers);
    NcDim nc_nRLc = ncFileLC.addDim("nR", nR_Lc);
    dims_lc.push_back(nc_nRLc);
// THIS WILL BREAK IF nY_LC is not declared
#if USE3DTETGEOM
    NcDim nc_nYLc = ncFileLC.addDim("nY", nY_Lc);
    dims_lc.push_back(nc_nYLc);
#endif

    NcDim nc_nZLc = ncFileLC.addDim("nZ", nZ_Lc);
    dims_lc.push_back(nc_nZLc);

    NcVar nc_Lc = ncFileLC.addVar("Lc", netcdf_precision, dims_lc);
    NcVar nc_s = ncFileLC.addVar("s", netcdf_precision, dims_lc);
    NcVar nc_ftx = ncFileLC.addVar("fx", netcdf_precision, dims_lc);
    NcVar nc_fty = ncFileLC.addVar("fy", netcdf_precision, dims_lc);
    NcVar nc_ftz = ncFileLC.addVar("fz", netcdf_precision, dims_lc);
    NcVar nc_btx = ncFileLC.addVar("bx", netcdf_precision, dims_lc);
    NcVar nc_bty = ncFileLC.addVar("by", netcdf_precision, dims_lc);
    NcVar nc_btz = ncFileLC.addVar("bz", netcdf_precision, dims_lc);
    NcVar nc_nI = ncFileLC.addVar("noIntersection", netcdf_precision, dims_lc);
    NcVar nc_gridRLc = ncFileLC.addVar("gridR", netcdf_precision, nc_nRLc);
#if USE3DTETGEOM
    NcVar nc_gridYLc = ncFileLC.addVar("gridY", netcdf_precision, nc_nYLc);
#endif
    NcVar nc_gridZLc = ncFileLC.addVar("gridZ", netcdf_precision, nc_nZLc);
   //FIXME - commented these because of disrupted workflow compile errors
    //nc_Lc.putVar(&Lc[0]);
    //nc_s.putVar(&s[0]);
    //nc_ftx.putVar(&forwardTracerX[0]);
    //nc_fty.putVar(&forwardTracerY[0]);
    //nc_ftz.putVar(&forwardTracerZ[0]);
    //nc_btx.putVar(&backwardTracerX[0]);
    //nc_bty.putVar(&backwardTracerY[0]);
    //nc_btz.putVar(&backwardTracerZ[0]);
    //nc_nI.putVar(&noIntersectionNodes[0]);
    //nc_gridRLc.putVar(&gridRLc[0]);
#if USE3DTETGEOM
    nc_gridYLc.putVar(&gridYLc[0]);
#endif
    nc_gridZLc.putVar(&gridZLc[0]);
    ncFileLC.close();
}
    /* END BLOCK: write_connection_length(...) - easy - netcdf  */


/* START BLOCK: write_flow(...) - easy - netcdf - return struct with xp and vp values -  */
 void write_flow(int n_flowV, int nR_flowV, int nY_flowV, int nZ_flowV, sim::Array<gitr_precision> &flowVr, sim::Array<gitr_precision> &flowVt, sim::Array<gitr_precision> &flowVz)
 { NcFile ncFileFlow("flowV.nc", NcFile::replace);
 NcDim nFlowV = ncFileFlow.addDim("n_flowV", n_flowV);
 NcDim nc_nRflow = ncFileFlow.addDim("nR", nR_flowV);
 NcDim nc_nYflow = ncFileFlow.addDim("nY", nY_flowV);
 NcDim nc_nZflow = ncFileFlow.addDim("nZ", nZ_flowV);
 vector<NcDim> dimsFlowV;
 dimsFlowV.push_back(nc_nZflow);
 dimsFlowV.push_back(nc_nYflow);
 dimsFlowV.push_back(nc_nRflow);
 NcVar nc_flowVr = ncFileFlow.addVar("flowVr", netcdf_precision, dimsFlowV);
 NcVar nc_flowVt = ncFileFlow.addVar("flowVt", netcdf_precision, dimsFlowV);
 NcVar nc_flowVz = ncFileFlow.addVar("flowVz", netcdf_precision, dimsFlowV);
 nc_flowVr.putVar(&flowVr[0]);
 nc_flowVt.putVar(&flowVt[0]);
 nc_flowVz.putVar(&flowVz[0]);
 ncFileFlow.close();
 /* remove that dumping into matlab format. useless*/
/*  std::string outnameFlowVr = "flowVr.m";
 std::string outnameFlowVz = "flowVz.m";
 std::string outnameFlowVt = "flowVt.m";
#if LC_INTERP == 3
 OUTPUT3d(profiles_folder, outnameFlowVr, nR_flowV, nY_flowV, nZ_flowV,
          &flowVr.front());
 OUTPUT3d(profiles_folder, outnameFlowVz, nR_flowV, nY_flowV, nZ_flowV,
          &flowVz.front());
 OUTPUT3d(profiles_folder, outnameFlowVt, nR_flowV, nY_flowV, nZ_flowV,
          &flowVt.front());
#else
 OUTPUT2d(profiles_folder, outnameFlowVr, nR_flowV, nZ_flowV, &flowVr.front());
 OUTPUT2d(profiles_folder, outnameFlowVz, nR_flowV, nZ_flowV, &flowVz.front());
 OUTPUT2d(profiles_folder, outnameFlowVt, nR_flowV, nZ_flowV, &flowVt.front());
#endif
*/
}
/* END BLOCK: write_flow(...) - easy - netcdf - return struct with xp and vp values -  */


 /* START BLOCK: read_particle(...) - medium - netcdf - return struct with xp and vp values -  */
 void read_particles(libconfig::Config cfg, vector<gitr_precision> xpfile, vector<gitr_precision> ypfile, vector<gitr_precision> zpfile, vector<gitr_precision> vxpfile, vector<gitr_precision>vypfile,vector<gitr_precision>vzpfile)
     {getVariable(cfg, "particleSource.ncFileString", ncParticleSourceFile);
     std::cout << "About to try to open NcFile ncp0 " << std::endl;
     // Return this in event of a problem.
     static const int NC_ERR = 2;
     try {
 	    netCDF::NcFile ncp0("input/" + ncParticleSourceFile, netCDF::NcFile::read);
     } catch (netCDF::exceptions::NcException &e) {
       e.what();
       cout << "FAILURE*************************************" << endl;
       return NC_ERR;
     }
     std::cout << "finished NcFile ncp0 starting ncp" << std::endl;
     netCDF::NcFile ncp("input/" + ncParticleSourceFile, netCDF::NcFile::read);
     std::cout << "getting dim nP" << std::endl;
     netCDF::NcDim ps_nP(ncp.getDim("nP"));

     nPfile = ps_nP.getSize();
     xpfile.resize(nPfile);
     ypfile.resize(nPfile);
     zpfile.resize(nPfile);
     vxpfile.resize(nPfile);
     vypfile.resize(nPfile);
     vzpfile.resize(nPfile);
     // std::cout << "nPfile "<< nPfile << std::endl;
     netCDF::NcVar ncp_x(ncp.getVar("x"));
     netCDF::NcVar ncp_y(ncp.getVar("y"));
     netCDF::NcVar ncp_z(ncp.getVar("z"));
     netCDF::NcVar ncp_vx(ncp.getVar("vx"));
     netCDF::NcVar ncp_vy(ncp.getVar("vy"));
     netCDF::NcVar ncp_vz(ncp.getVar("vz"));
     std::cout << "got through NcVar " << std::endl;
     ncp_x.getVar(&xpfile[0]);
     ncp_y.getVar(&ypfile[0]);
     ncp_z.getVar(&zpfile[0]);
     ncp_vx.getVar(&vxpfile[0]);
     ncp_vy.getVar(&vypfile[0]);
     ncp_vz.getVar(&vzpfile[0]);
     std::cout << "defined file vectors " << std::endl;
     ncp.close();
     std::cout << "closed ncp " << std::endl;
     }
     /* END BLOCK: read_particle(...) - medium - netcdf - return struct with xp and vp values -  */


 /* START BLOCK: write_particle(...) - easy - netcdf  */

 void write_particles(int nP, int pNP, sim::Array<gitr_precision> &pSurfNormX, sim::Array<gitr_precision> &pSurfNormY, sim::Array<gitr_precision> &pSurfNormZ, sim::Array<gitr_precision> &pvx, sim::Array<gitr_precision> &pvy, sim::Array<gitr_precision> &pvz, sim::Array<gitr_precision> &px, sim::Array<gitr_precision> &py, sim::Array<gitr_precision> &pz)
 {
 	std::cout << "writing particles out file" << std::endl;
     netCDF::NcFile ncFile_particles("output/particleSource.nc", netCDF::NcFile::replace);
     netCDF::NcDim pNP = ncFile_particles.addDim("nP", nP);
     netCDF::NcVar p_surfNormx = ncFile_particles.addVar("surfNormX", netcdf_precision, pNP);
     netCDF::NcVar p_surfNormy = ncFile_particles.addVar("surfNormY", netcdf_precision, pNP);
     netCDF::NcVar p_surfNormz = ncFile_particles.addVar("surfNormZ", netcdf_precision, pNP);
     netCDF::NcVar p_vx = ncFile_particles.addVar("vx", netcdf_precision, pNP);
     netCDF::NcVar p_vy = ncFile_particles.addVar("vy", netcdf_precision, pNP);
     netCDF::NcVar p_vz = ncFile_particles.addVar("vz", netcdf_precision, pNP);
     netCDF::NcVar p_x = ncFile_particles.addVar("x", netcdf_precision, pNP);
     netCDF::NcVar p_y = ncFile_particles.addVar("y", netcdf_precision, pNP);
     netCDF::NcVar p_z = ncFile_particles.addVar("z", netcdf_precision, pNP);
     p_surfNormx.putVar(&pSurfNormX[0]);
     p_surfNormy.putVar(&pSurfNormY[0]);
     p_surfNormz.putVar(&pSurfNormZ[0]);
     p_vx.putVar(&pvx[0]);
     p_vy.putVar(&pvy[0]);
     p_vz.putVar(&pvz[0]);
     p_x.putVar(&px[0]);
     p_y.putVar(&py[0]);
     p_z.putVar(&pz[0]);
     ncFile_particles.close();
     std::cout << "finished writing particles out file" << std::endl;
 }
     /* END BLOCK: write_particle(...) - easy - netcdf  */


 /* BLOCK: write_force(...) - easy - netcdf  */
 void write_force(int nR_force,int nZ_force, sim::Array<gitr_precision> &forceR, sim::Array<gitr_precision> &forceZ,sim::Array<gitr_precision>  &tIon,sim::Array<gitr_precision>  &tRecomb,sim::Array<gitr_precision>  &dvEr,sim::Array<gitr_precision>  &dvEz,sim::Array<gitr_precision> &dvEt,sim::Array<gitr_precision>  &dvBr,sim::Array<gitr_precision> &dvBz,sim::Array<gitr_precision> &dvBt,sim::Array<gitr_precision> &dvCollr,sim::Array<gitr_precision> &dvCollz,sim::Array<gitr_precision> &dvCollt,sim::Array<gitr_precision> &dvITGr,sim::Array<gitr_precision> &dvITGz,sim::Array<gitr_precision> &dvITGt,sim::Array<gitr_precision> &dvEr,sim::Array<gitr_precision> &dvEz,sim::Array<gitr_precision> &dvEt)
 { std::cout << " about to write ncFile_forces " << std::endl;


    netCDF::NcFile ncFile_force("output/forces.nc", netCDF::NcFile::replace);
    netCDF::NcDim nc_nRf = ncFile_force.addDim("nR", nR_force);
    netCDF::NcDim nc_nZf = ncFile_force.addDim("nZ", nZ_force);
    vector<netCDF::NcDim> forceDims;
    forceDims.push_back(nc_nZf);
    forceDims.push_back(nc_nRf);
    netCDF::NcVar forceRf = ncFile_force.addVar("r", netcdf_precision, nc_nRf);
    netCDF::NcVar forceZf = ncFile_force.addVar("z", netcdf_precision, nc_nZf);
    netCDF::NcVar nction = ncFile_force.addVar("tIon", netcdf_precision, forceDims);
    netCDF::NcVar nctrec = ncFile_force.addVar("tRec", netcdf_precision, forceDims);
    netCDF::NcVar dvErf = ncFile_force.addVar("dvEr", netcdf_precision, forceDims);
    netCDF::NcVar dvEzf = ncFile_force.addVar("dvEz", netcdf_precision, forceDims);
    netCDF::NcVar dvEtf = ncFile_force.addVar("dvEt", netcdf_precision, forceDims);
    netCDF::NcVar dvBrf = ncFile_force.addVar("dvBr", netcdf_precision, forceDims);
    netCDF::NcVar dvBzf = ncFile_force.addVar("dvBz", netcdf_precision, forceDims);
    netCDF::NcVar dvBtf = ncFile_force.addVar("dvBt", netcdf_precision, forceDims);
    netCDF::NcVar dvCollrf = ncFile_force.addVar("dvCollr", netcdf_precision, forceDims);
    netCDF::NcVar dvCollzf = ncFile_force.addVar("dvCollz", netcdf_precision, forceDims);
    netCDF::NcVar dvColltf = ncFile_force.addVar("dvCollt", netcdf_precision, forceDims);
    netCDF::NcVar dvITGrf = ncFile_force.addVar("dvITGr", netcdf_precision, forceDims);
    netCDF::NcVar dvITGzf = ncFile_force.addVar("dvITGz", netcdf_precision, forceDims);
    netCDF::NcVar dvITGtf = ncFile_force.addVar("dvITGt", netcdf_precision, forceDims);
    netCDF::NcVar dvETGrf = ncFile_force.addVar("dvETGr", netcdf_precision, forceDims);
    netCDF::NcVar dvETGzf = ncFile_force.addVar("dvETGz", netcdf_precision, forceDims);
    netCDF::NcVar dvETGtf = ncFile_force.addVar("dvETGt", netcdf_precision, forceDims);
    forceRf.putVar(&forceR[0]);
    forceZf.putVar(&forceZ[0]);
    nction.putVar(&tIon[0]);
    nctrec.putVar(&tRecomb[0]);
    dvErf.putVar(&dvEr[0]);
    dvEzf.putVar(&dvEz[0]);
    dvEtf.putVar(&dvEt[0]);
    dvBrf.putVar(&dvBr[0]);
    dvBzf.putVar(&dvBz[0]);
    dvBtf.putVar(&dvBt[0]);
    dvCollrf.putVar(&dvCollr[0]);
    dvCollzf.putVar(&dvCollz[0]);
    dvColltf.putVar(&dvCollt[0]);
    dvITGrf.putVar(&dvITGr[0]);
    dvITGzf.putVar(&dvITGz[0]);
    dvITGtf.putVar(&dvITGt[0]);
    dvETGrf.putVar(&dvETGr[0]);
    dvETGzf.putVar(&dvETGz[0]);
    dvETGtf.putVar(&dvETGt[0]);
    ncFile_force.close();
 }
    /* END BLOCK: write_force(...) */

 /* START BLOCK: write_position() - easy - netcdf  */
 void write_position(int nP, sim::Array<gitr_precision> &xGather,sim::Array<gitr_precision> &yGather,sim::Array<gitr_precision> &zGather,sim::Array<gitr_precision> &vxGather,sim::Array<gitr_precision> &vyGather,sim::Array<gitr_precision> &vzGather,sim::Array<gitr_precision> &hitWallGather,sim::Array<gitr_precision> &surfaceHitGather,sim::Array<gitr_precision> &weightGather,sim::Array<gitr_precision> &chargeGather,sim::Array<gitr_precision> &hasLeakedGather,Particles &particlearray)
     /* should add a subsampling option. Dumping and plotting a billion of particles is kind of expansive... Parameters should be a list of indexes of particles that need be dumped */
     /* write_particle_position(int nP, )*/
     {
	 netCDF::NcFile ncFile0("output/positions.nc", netCDF::NcFile::replace);
     netCDF::NcDim nc_nP0 = ncFile0.addDim("nP", nP);
     vector<netCDF::NcDim> dims0;
     dims0.push_back(nc_nP0);

     netCDF::NcVar nc_x0 = ncFile0.addVar("x", netcdf_precision, dims0);
     netCDF::NcVar nc_y0 = ncFile0.addVar("y", netcdf_precision, dims0);
     netCDF::NcVar nc_z0 = ncFile0.addVar("z", netcdf_precision, dims0);
     netCDF::NcVar nc_vx0 = ncFile0.addVar("vx", netcdf_precision, dims0);
     netCDF::NcVar nc_vy0 = ncFile0.addVar("vy", netcdf_precision, dims0);
     netCDF::NcVar nc_vz0 = ncFile0.addVar("vz", netcdf_precision, dims0);
     netCDF::NcVar nc_trans0 = ncFile0.addVar("transitTime", netcdf_precision, dims0);
     netCDF::NcVar nc_impact0 = ncFile0.addVar("hitWall", netcdf_precision, dims0);
     netCDF::NcVar nc_surfHit0 = ncFile0.addVar("surfaceHit", netCDF::ncInt, dims0);
     netCDF::NcVar nc_weight0 = ncFile0.addVar("weight", netcdf_precision, dims0);
     netCDF::NcVar nc_charge0 = ncFile0.addVar("charge", netcdf_precision, dims0);
     netCDF::NcVar nc_leak0 = ncFile0.addVar("hasLeaked", netCDF::ncInt, dims0);
     netCDF::NcVar nc_dist0 = ncFile0.addVar("distTraveled", netcdf_precision, dims0);
     netCDF::NcVar nc_time0 = ncFile0.addVar("time", netcdf_precision, dims0);
     netCDF::NcVar nc_dt0 = ncFile0.addVar("dt", netcdf_precision, dims0);
 #if USE_MPI > 0
     // could use a structure to collect data! simplicity would suggest to have one instance of structure (non-allocated) per mpi process then allocate each on each process to limit memory issue.
      * then use an operator overloading to broadcast and sum all collector array together and dump them.
      * could also iterate over the attributes with boost instead of writing tons of lines to dump the data!
      * see https://stackoverflow.com/questions/17660095/iterating-over-a-struct-in-c
      */
     nc_x0.putVar(&xGather[0]);
     nc_y0.putVar(&yGather[0]);
     nc_z0.putVar(&zGather[0]);
     nc_vx0.putVar(&vxGather[0]);
     nc_vy0.putVar(&vyGather[0]);
     nc_vz0.putVar(&vzGather[0]);
     nc_trans0.putVar(&particleArray->transitTime[0]);
     nc_impact0.putVar(&hitWallGather[0]);
     nc_surfHit0.putVar(&surfaceHitGather[0]);
     nc_weight0.putVar(&weightGather[0]);
     nc_charge0.putVar(&chargeGather[0]);
     nc_leak0.putVar(&hasLeakedGather[0]);
 #else
   std::cout << "not using mpi output" << std::endl;
   nc_x0.putVar(&particleArray->xprevious[0]);
   nc_y0.putVar(&particleArray->yprevious[0]);
   nc_z0.putVar(&particleArray->zprevious[0]);
   nc_vx0.putVar(&particleArray->vx[0]);
   nc_vy0.putVar(&particleArray->vy[0]);
   nc_vz0.putVar(&particleArray->vz[0]);
   nc_trans0.putVar(&particleArray->transitTime[0]);
   nc_impact0.putVar(&particleArray->hitWall[0]);
   nc_surfHit0.putVar(&particleArray->surfaceHit[0]);
   nc_weight0.putVar(&particleArray->weight[0]);
   nc_charge0.putVar(&particleArray->charge[0]);
   nc_leak0.putVar(&particleArray->hasLeaked[0]);
   nc_dist0.putVar(&particleArray->distTraveled[0]);
   nc_time0.putVar(&particleArray->time[0]);
   nc_dt0.putVar(&particleArray->dt[0]);
 #endif
     ncFile0.close();
     }
     /* END BLOCK: write_position(...) - easy - netcdf  */
 /* START BLOCK: write_surface_data(...) - easy - netcdf  */
 void write_surface_data(int nSurfaces, int nEdist, int nAdist, int dimsSurf,sim::Array<gitr_precision>  &grossDeposition,sim::Array<gitr_precision>  &surfaceNumbers,sim::Array<gitr_precision>  &grossErosion,sim::Array<gitr_precision>  &aveSputtYld,sim::Array<gitr_precision>  &sputtYldCount,sim::Array<gitr_precision>  &sumParticlesStrike,sim::Array<gitr_precision>  &sumWeightStrike,sim::Array<gitr_precision>  &energyDistribution,sim::Array<gitr_precision>  &reflDistribution,sim::Array<gitr_precision>  &sputtDistribution)
 {
     netCDF::NcFile ncFile1("output/surface.nc", netCDF::NcFile::replace);
     netCDF::NcDim nc_nLines = ncFile1.addDim("nSurfaces", nSurfaces);
     vector<netCDF::NcDim> dims1;
     dims1.push_back(nc_nLines);

     vector<netCDF::NcDim> dimsSurfE;
     dimsSurfE.push_back(nc_nLines);
     netCDF::NcDim nc_nEnergies = ncFile1.addDim("nEnergies", nEdist);
     netCDF::NcDim nc_nAngles = ncFile1.addDim("nAngles", nAdist);
     dimsSurfE.push_back(nc_nAngles);
     dimsSurfE.push_back(nc_nEnergies);
     netCDF::NcVar nc_grossDep = ncFile1.addVar("grossDeposition", netcdf_precision, nc_nLines);
     netCDF::NcVar nc_grossEro = ncFile1.addVar("grossErosion", netcdf_precision, nc_nLines);
     netCDF::NcVar nc_aveSpyl = ncFile1.addVar("aveSpyl", netcdf_precision, nc_nLines);
     netCDF::NcVar nc_spylCounts = ncFile1.addVar("spylCounts", netCDF::ncInt, nc_nLines);
     netCDF::NcVar nc_surfNum = ncFile1.addVar("surfaceNumber", netCDF::ncInt, nc_nLines);
     netCDF::NcVar nc_sumParticlesStrike =
         ncFile1.addVar("sumParticlesStrike", netCDF::ncInt, nc_nLines);
     netCDF::NcVar nc_sumWeightStrike =
         ncFile1.addVar("sumWeightStrike", netcdf_precision, nc_nLines);
     nc_grossDep.putVar(&grossDeposition[0]);
     nc_surfNum.putVar(&surfaceNumbers[0]);
     nc_grossEro.putVar(&grossErosion[0]);
     nc_aveSpyl.putVar(&aveSputtYld[0]);
     nc_spylCounts.putVar(&sputtYldCount[0]);
     nc_sumParticlesStrike.putVar(&sumParticlesStrike[0]);
     nc_sumWeightStrike.putVar(&sumWeightStrike[0]);
     // NcVar nc_surfImpacts = ncFile1.addVar("impacts",netcdf_precision,dims1);
     // NcVar nc_surfRedeposit = ncFile1.addVar("redeposit",netcdf_precision,dims1);
     // NcVar nc_surfStartingParticles =
     // ncFile1.addVar("startingParticles",netcdf_precision,dims1); NcVar nc_surfZ =
     // ncFile1.addVar("Z",netcdf_precision,dims1);
     netCDF::NcVar nc_surfEDist = ncFile1.addVar("surfEDist", netcdf_precision, dimsSurfE);
     netCDF::NcVar nc_surfReflDist = ncFile1.addVar("surfReflDist", netcdf_precision, dimsSurfE);
     netCDF::NcVar nc_surfSputtDist =
         ncFile1.addVar("surfSputtDist", netcdf_precision, dimsSurfE);
     // nc_surfImpacts.putVar(impacts);
     //#if USE3DTETGEOM > 0
     // nc_surfRedeposit.putVar(redeposit);
     //#endif
     // nc_surfStartingParticles.putVar(startingParticles);
     // nc_surfZ.putVar(surfZ);
     nc_surfEDist.putVar(&energyDistribution[0]);
     nc_surfReflDist.putVar(&reflDistribution[0]);
     nc_surfSputtDist.putVar(&sputtDistribution[0]);
     // NcVar nc_surfEDistGrid = ncFile1.addVar("gridE",ncDouble,nc_nEnergies);
     // nc_surfEDistGrid.putVar(&surfaces->gridE[0]);
     // NcVar nc_surfADistGrid = ncFile1.addVar("gridA",ncDouble,nc_nAngles);
     // nc_surfADistGrid.putVar(&surfaces->gridA[0]);
     ncFile1.close();
     }
     /* END BLOCK: write_surface_data(...) - easy - netcdf  */


 /* START BLOCK: write_history_data(...) - easy - netcdf  */
    // Write netCDF output for histories
 void write_history(int nHistoriesPerParticle, int nP,sim::Array<gitr_precision>  &positionHistoryXgather,sim::Array<gitr_precision>  &positionHistoryYgather,sim::Array<gitr_precision>  &positionHistoryZgather,sim::Array<gitr_precision>  &velocityHistorygather,sim::Array<gitr_precision>  &velocityHistoryXgather,sim::Array<gitr_precision>  &velocityHistoryYgather,sim::Array<gitr_precision>  &velocityHistoryZgather,sim::Array<gitr_precision> &chargeHistory)
    {
	 netCDF::NcFile ncFile_hist("output/history.nc", netCDF::NcFile::replace);
    netCDF::NcDim nc_nT = ncFile_hist.addDim("nT", nHistoriesPerParticle);
    netCDF::NcDim nc_nP = ncFile_hist.addDim("nP", nP);
    vector<netCDF::NcDim> dims_hist;
    dims_hist.push_back(nc_nP);
    dims_hist.push_back(nc_nT);
    // NcDim nc_nPnT = ncFile_hist.addDim("nPnT",nP*nT/subSampleFac);
    // dims_hist.push_back(nc_nPnT);
    netCDF::NcVar nc_x = ncFile_hist.addVar("x", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_y = ncFile_hist.addVar("y", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_z = ncFile_hist.addVar("z", netCDF::ncDouble, dims_hist);

    netCDF::NcVar nc_v = ncFile_hist.addVar("v", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_vx = ncFile_hist.addVar("vx", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_vy = ncFile_hist.addVar("vy", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_vz = ncFile_hist.addVar("vz", netCDF::ncDouble, dims_hist);

    netCDF::NcVar nc_charge = ncFile_hist.addVar("charge", netCDF::ncDouble, dims_hist);
    netCDF::NcVar nc_weight = ncFile_hist.addVar("weight", netCDF::ncDouble, dims_hist);
#if USE_MPI > 0
    // if(world_rank ==0)
    //{
    // for(int i=0;i<401;i++)
    //{
    //  std::cout << "Rank " << world_rank << "z " << positionHistoryZgather[i]
    //  << std::endl;
    //}
    //}
    nc_x.putVar(&positionHistoryXgather[0]);
    nc_y.putVar(&positionHistoryYgather[0]);
    nc_z.putVar(&positionHistoryZgather[0]);

    nc_v.putVar(&velocityHistorygather[0]);
    nc_vx.putVar(&velocityHistoryXgather[0]);
    nc_vy.putVar(&velocityHistoryYgather[0]);
    nc_vz.putVar(&velocityHistoryZgather[0]);

    nc_charge.putVar(&chargeHistoryGather[0]);
    nc_weight.putVar(&weightHistoryGather[0]);
#else
    nc_x.putVar(&positionHistoryX[0]);
    nc_y.putVar(&positionHistoryY[0]);
    nc_z.putVar(&positionHistoryZ[0]);

    nc_vx.putVar(&velocityHistoryX[0]);
    nc_vy.putVar(&velocityHistoryY[0]);
    nc_vz.putVar(&velocityHistoryZ[0]);

    nc_charge.putVar(&chargeHistory[0]);
#endif
    ncFile_hist.close();
    }
    /* END BLOCK: write_history_data(...) - easy - netcdf  */
