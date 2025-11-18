#include "mcmodel.h"

void MeshQuantities::GetSurfScRate()
{
  double *roughnessRate, *surfPhononRate;

  if (par_type == PELEC) {
    c_e_roughnessRate->ExtractView(&roughnessRate);
    c_e_surfPhononRate->ExtractView(&surfPhononRate);

    SurfRoughnessScRate = roughnessRate[C_LINDEX_GHOST_ONE(icell,jcell,kcell)]; 
    
    SurfPhononScRate = surfPhononRate[C_LINDEX_GHOST_ONE(icell, jcell, kcell)]; 

    SurfScGamma = SurfRoughnessScRate + SurfPhononScRate;
  } else {
    c_h_roughnessRate->ExtractView(&roughnessRate);
    c_h_surfPhononRate->ExtractView(&surfPhononRate);

    SurfRoughnessScRate = roughnessRate[C_LINDEX_GHOST_ONE(icell,jcell,kcell)]; 
    
    SurfPhononScRate = surfPhononRate[C_LINDEX_GHOST_ONE(icell, jcell, kcell)]; 

    SurfScGamma = SurfRoughnessScRate + SurfPhononScRate;
  }
  
   if(Flag_calSurfscatt==1 && mpi_rank==0)
  {
     ofstream ofile;
	 string filename;
	 filename="./data/test_SurfSc.txt";
	 ofile.open(filename.c_str());
	 ofile<<"SurfRoughnessScRate="<<SurfRoughnessScRate*pow(field0,2)<<endl; 
	 ofile<<"SurfPhononScRate="<<SurfPhononScRate*pow(field0,2)<<endl; 
	 ofile<<"SurfScGamma="<<SurfScGamma*pow(field0,2)<<endl; 
	 ofile<<"SurfSc_delta"<<SurfSc_delta*spr0<<endl;
	 ofile<<"SurfSc_ail"<<SurfSc_ail*spr0<<endl;
	 ofile<<"icell="<<icell<<endl;
	 ofile<<"jcell="<<jcell<<endl;
	 ofile<<"kcell="<<kcell<<endl;
	 ofile.close();
	 exit(0);
  }
  
}


void MeshQuantities::ParticleSurfaceScatter()
{
	if(Random() * SurfScGamma < SurfRoughnessScRate){
	        idir = SurfaceDir[(*c_nearestSurf)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)]];
		Diffuse();
	}
      	else
      	{
   		if(par_type==PELEC)
              ElectronSurfacePhononScatter();
      		else if(par_type==PHOLE)
              HoleSurfacePhononScatter();
      	}
}


void MeshQuantities::HoleSurfacePhononScatter() 
{
    //Purpose:   performs hole Surface scattering//S/fSc1ugamma_wvce\R
	//kHolePhononScatter()\NgRvce\$Re 
	//V:Nhole surface phonon scattering rateNesQ (WN[eN[VQ/f8^pe
    
    //local variables
    int iscat,ibold;
    double drand,dscat;
    double sca[MNScaHole*NBH];

    //only holes are allowed in this routine
    if(par_type!=PHOLE)
    {

        cout<<"error HoleSurfacePhononScatter: Wrong particle type";exit(0);
    }
    //save intial band
    ibold=iband;

    //scattering statistics

    Flag_SelfScatter =false;

        //calculate scattering-rate for all prozesses at energy energy
    band.CALSCATTH(energy,sca,iband);

        //select the scattering prozess with AcRejection
    iscat=-1;
    dscat=0.0;
    while(drand>dscat)
      {
        
        iscat=iscat+1;
        if(iscat>=(band.scprh*band.nband[PHOLE]))
          {
            
            Flag_SelfScatter =true;
            break;
            
          }
        dscat=dscat+sca[iscat];
      }

        //iscat is number of chosen scatt.pro.
    if(!Flag_SelfScatter)
      if(((iscat%band.scprh)+1)==band.scprh)
        HSCATII(iscat);
      else
        HPSCAT(iscat);
}


void MeshQuantities::ElectronSurfacePhononScatter()
{
    //Purpose:   performs electron surface phonon scattering
	//kElectronPhononScatter()\NgRvce\$Re 
	//V:Nelectron surface phonon scattering rateNesQ (WN[eN[VQ/f8^pe
    //local variables
    int iscat,ibold;
    double drand,dscat;
    double sca[MNScaEle*NBE];

    //only electrons are allowed in this routine
    if(par_type!=PELEC)
    {
	cout<<"error ElectronSurfacePhononScatter: Wrong particle type";
	exit(0);
    }
    //save inital band
    ibold=iband;

    //scattering statistics
        
    //first check for fictious scattering
    
    Flag_SelfScatter=false;
    
    //calculate scattering-rate for all prozesses at energy energy
    band.CALSCATTE(energy,sca,iband);

        //select the scattering prozess with AcRejection
    iscat=-1;
    dscat=0.0;
    
    while(drand>dscat)
      {
        iscat=iscat+1;
        if(iscat>=(band.scpre*band.nband[PELEC]))
          {
            
            Flag_SelfScatter=true;
            break;
            
          }
        dscat=dscat+sca[iscat];
        
      }
        if(!Flag_SelfScatter)
          {
            //iscat is number of chosen scatt.pro.
            
            if(((iscat%band.scpre)+1)==band.scpre)
            {
              ESCATII(iscat);
            }
            //Phonon
            else
            {
                if(band.jacophfl)
                {
                    EPSCAT(iscat);
                }
                else
                {
                    if(band.fiscphfl)
                    {
                        EFPSCAT(iscat);
                    }
                    else
                    {

                        cout<<"error Warning: ESCTR: Specify phonon type";exit(0);
                    }
                }
            }
        }
    
}

