
#include "mcmodel.h"

void MeshQuantities::GetSurfRoughnessPhononScRate()
{
   double Wcl,Wq,Pf,Wclq;
   double effmass;
   int i,j,iptype,k;
   int * EeffDirection;
   double *xfield, *yfield, *zfield;
   double *roughnessRate, *surfPhononRate;
   double Eeff;
   double *par_charge;

   c_EeffDirection->ExtractView(&EeffDirection);

   c_field_x->ExtractView(&xfield);
   c_field_y->ExtractView(&yfield);
   c_field_z->ExtractView(&zfield);
   c_par_charge->ExtractView(&par_charge);
    
   for(iptype=0;iptype< 2 ;iptype++){
        
     effmass=MASSeff[iptype];

     if (iptype == PELEC) {
       c_e_roughnessRate->ExtractView(&roughnessRate);
       c_e_surfPhononRate->ExtractView(&surfPhononRate);
     } else {
       c_h_roughnessRate->ExtractView(&roughnessRate);
       c_h_surfPhononRate->ExtractView(&surfPhononRate);
     }

     for (i = c_ibegin; i <= c_iend; i ++)
       for (k = c_kbegin; k <= c_kend; k ++)
	 for (j = c_jbegin_ghost; j <= c_jend_ghost; j ++){

	   Rho = par_charge[C_LINDEX_GHOST_ONE(i, j, k)];

	   switch(EeffDirection[C_LINDEX_GHOST_ONE(i,j,k)])
	      {	
	      case 0:Eeff = fabs(xfield[C_LINDEX_GHOST_ONE(i,j,k)]);break;
	      case 1:Eeff = fabs(yfield[C_LINDEX_GHOST_ONE(i,j,k)]);break;
	      case 2:Eeff = fabs(zfield[C_LINDEX_GHOST_ONE(i,j,k)]);break;
	      }
                
                  if(Eeff<1e-20) Eeff=1e-20;
							//  surface roughness scattering rate
                    if(!Flag_SurfaceRoughnessScatter)
                       roughnessRate[C_LINDEX_GHOST_ONE(i,j,k)] = 0;
                    else
                      roughnessRate[C_LINDEX_GHOST_ONE(i,j,k)] = 2 * PI * effmass *
                        pow(SurfSc_delta* SurfSc_ail * Eeff * eps[SILICON]/(eps[SILICON]+eps[OXIDE]),2);

                    //	surface phonon scattering rate
                    if(!Flag_SurfacePhononScatter)	
		                  surfPhononRate[C_LINDEX_GHOST_ONE(i,j,k)] = 0;
                    else
                      {
                        Wcl = pow(Tn, SurfSc_Nbmod) * BOLTZ*T0 / (Eeff*field0);
                        Wq = pow(12.0*effmass*em0*(Eeff*field0)/ec0/hq0/hq0, -1.0/3.0);
                        Pf = SurfSc_Pft*Tn + SurfSc_Pfn * pow(Rho * conc0/(1e17*1e6), SurfSc_Npf) / Tn;
                        Wclq=(SurfSc_gama*Wcl+SurfSc_theta*Wq)/Pf;
                        surfPhononRate[C_LINDEX_GHOST_ONE(i,j,k)] = effmass*em0*BOLTZ*T0*SurfSc_XIph*SurfSc_XIph
                          /(pow(hq0,3)*sirho*pow(9050.02,2)*Wclq)/scrt0
                        +ec0/(effmass*em0*SurfSc_Kpha*pow(Tn,2.24))/scrt0;
                      }
      }
  }
    
}
