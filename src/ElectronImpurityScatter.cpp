
#include "mcmodel.h"

/*find the maximum value of scattering rate in the tetrahedron
 * */
double MeshQuantities::GetImpScGamma()
{
	//Gamma
    //local variables
    double eemax,gammax,rv;
    
    /* the scattering rate f(e) as a function of energy has a maximum point
     * located at (md * eebeta / (4 * mt) , which can be obtained by solving f'(e) = 0 
     * analyticaly. -- Wei Zhang
     * */
      if(band.bhmrtfl)
    {
      /*1. do not miss the 4 * PI
       *2. the Rho should be density of carrier, not charge.
       */
        eemax=frickel* 4 * PI * Rho/(8.0*band.meld*eps[SILICON])*0.56875;
    }
    else
    {
      /*do not miss the 4 * PI 
       * notice we use mt here*/
        eemax=frickel*4 * PI * Rho/(8.0*band.melt*eps[SILICON]);
    }
    /* f(e) decreases when e > eemax, so maximum point is at
     * band.eek[band.tet[0][itet]]] */

    if(eemax<band.eek[band.tet[0][itet]])
    {
       gammax=GetImpScRate(band.eek[band.tet[0][itet]]);
    }
    else
    {

    /* f(e) decreases when e > eemax, so maximum point is at
     * band.eek[band.tet[3][itet]]] */

        if(eemax>band.eek[band.tet[3][itet]])
        {
          gammax=GetImpScRate(band.eek[band.tet[3][itet]]);
         }
        else
        {
          gammax=GetImpScRate(eemax);
         }
    }
    if(band.bhmrtfl)
    {
        gammax=1.01*gammax;
     }

    rv=Max(gammax,1.0/scrt0);

    return rv;
}


double MeshQuantities::GetImpScRate(double eel)
{
    double eebeta,rvscrt,zzz,Rsh,f2,f3;

    /* the scattering event will not happen when:
     * 1. the cell has a low dopping density, or 
     * 2. carrier' energy is higher than 0.12 eV.
     *
     * so I guess that impurity scattering is not important
     * in channel.  ----Wei Zhang
     * */

    if((eel*eV0<0.120)&&(DA * conc0>=1e22))
    {
      /*notice that there should have a coefficient: 4 * PI.
       * (if frickel equals 1), by Wei Zhang
       * */
        eebeta=frickel *  4 * PI * Rho / (2.0*band.meld * eps[SILICON]);
	
        if(band.bhmrtfl)
        {
            eel=Max(eel,(1e-5));
            rvscrt=DA/(16.0*PI*pow(eps[SILICON],2.0)
                *sqrt(2.0*band.meld*pow(eel,3.0)))
                * (log((eebeta+4.0*eel)/eebeta)
                - 4.0*eel/(eebeta+4.0*eel));
        }
        else 
        {
          if(Flag_SurfaceScatter &&  Flag_SurfaceImpurityScatter && ((*c_InSurfRegion)[C_LINDEX_GHOST_ONE(icell, jcell,kcell)])) 
	        {
		  /*
              int CellNearestSurface_tmp = (*c_CellNearestSurface)[C_LINDEX_GHOST_ONE(icell, jcell,kcell)];
              
              switch(SurfaceType[CellNearestSurface_tmp])
			{	
            case 0:zzz=fabs(x-SurfacePosition[CellNearestSurface_tmp]) * spr0 * 1e9;break;
            case 1:zzz=fabs(y-SurfacePosition[CellNearestSurface_tmp]) * spr0 * 1e9;break;
			}
              Rsh=SurfSc_Rshmin+SurfSc_Rshmax/(1+exp((1e18*1e6-DA*conc0)/(SurfSc_GAMMAn*1e18*1e6)));
              f2=-Rsh*pow(zzz,6.0)/((pow(zzz,4.0)+pow(1.0,4.0))*pow(1.0,2.0));
              f3=1-exp(f2);
              if(fabs(f3)<1e-20)f3=1e-20;
              rvscrt=SurfSc_Fsimp*DA*band.meld*sqrt(2.0*band.meld*eel)/(PI*pow(frickel*Rho/f3,2.0)
                                                                        *(1.0+4.0*eel*band.melt/(eebeta*band.meld)));
									*/
            }	
          else
	    /*I have no idea what the if sentence does, but for the else part
	     * refer to (1) on page 682 , equation (3.100), which describe the  
	     * scattering rate of ionized impurity as a function of energy.
	     * the Brooks and Herring(BH) model is used and deduced by
	     * parabolic band. ---- Wei Zhang
	     * */
            rvscrt=DA*band.meld*sqrt(2.0*band.meld*eel)/(PI*pow(frickel*Rho,2.0)
	                	*(1.0+4.0*eel*band.melt/(eebeta*band.meld)));

        }
        
    }
    else
      rvscrt=0.0;

    return rvscrt;
}

void MeshQuantities::ElectronImpurityScatter()
{

    //Purpose:   performs electron impurity scattering

    //local variables
    int itetold,isymold;
    double gamimp;

    //only electrons are allowed in this routine
    if(par_type != PELEC)
    {

    }
    //scattering statistics
    //check maximum scattering rate
    
    gamimp=GetImpScRate(energy);

    if(gamimp>(ImpScGamma + 1e-12))
    {
      
    	cout << "----------------------------------------------------------"<<endl;
        cout << "error EBHSCTR: Variable Gamma wrong"<<endl;
        cout << "error EBHSCTR: Variable Gamma wrong";
        //exit(0);
    }
    //first check for fictious scattering
    if(gamimp<Random()*ImpScGamma)
    {
        Flag_SelfScatter=true;
    }
    else
    {
        //Chose final state
        if(band.bhmrtfl)
        {
            isymold=isym;
            itetold=itet;
            TALIND(1);
            get_state();
            if(Flag_SelfScatter)
            {
                isym=isymold;
                itet=itetold;
            }
            else
            {
                //GETGVE(band);
            }
        }
        else
        {
            EleImpScUpdateK();
        }
    }

}


void MeshQuantities::EleImpScUpdateK()
{
    //Purpose:   performs electron impurity scattering
    //ElectronImpurityScatter
    //local variables
    double xkl,ykl,zkl,ca,cb;
    double xkll,ykll,zkll;
    double xklll,yklll,zklll;
    double cosphi,sinphi,costheta,sintheta;
    double cospr,sinpr,costr,sintr;
    double eebeta,pr,alfa,betaq;
    double dn,dnmax;

    //Transform kcell-vector into irreducible wedge
    in_wedge(xkl,ykl,zkl,kx,ky,kz);

    //calculate kcell-vector relative to band minimum
    
    xkl=xkl-0.85*a0pi;

    //apply BH-transformation
    xkl=xkl*sqrt(band.meld/band.mell);
    ykl=ykl*sqrt(band.meld/band.melt);
    zkl=zkl*sqrt(band.meld/band.melt);

    //screening energy
    
    /* the betaq miss the coefficient 4 * PI, I added it.
     * ----Wei Zhang */
    betaq=frickel*4 * PI * Rho/eps[SILICON];
    eebeta=betaq/(2.0*band.meld);

    /* the following 3 lines coincide with the cos expression 
     * of scattering angle theta , given in the review article
     * (1) on page 683 , equation (3.104).  ----Wei Zhang
     */
    //cal. angle between old and new kcell-vector
    alfa=eebeta*band.meld/(2.0*energy*band.melt);
    pr=Random()/(1.0+0.5*alfa);
    costr=1.0-alfa*pr/(1.0-pr);

    sintr=sqrt(1.0-costr*costr);

    //cal. angle of the new kcell-vector in the plane perpendicular to
    //the old kcell-vector
    
    /*"the azimuthal angle of k' around k is chosen at random" 
     *   ----cited from (1) on page 683  */
    pr=2.0*PI*Random();
    cospr=cos(pr);
    sinpr=sin(pr);

    //cal. new kcell-vector in a coordinate frame, where the old kcell-vector points
    //into the positve z-direction
    /*correct, calculate the new vector in a coordinate frame whose z-axis is
     * parallel to old kcell-vector ----Wei Zhang */ 
    ca=sqrt(xkl*xkl+ykl*ykl);
    cb=sqrt(ca*ca+zkl*zkl);
    xkll=cb*sintr*cospr;
    ykll=cb*sintr*sinpr;
    zkll=cb*costr;

    //transform new kcell-vector into the BH-transformed coordinate frame
    cosphi=xkl/ca;
    sinphi=ykl/ca;

    costheta=zkl/cb;
    sintheta=ca/cb;

    xklll=cosphi*costheta*xkll-sinphi*ykll+cosphi*sintheta*zkll;
    yklll=sinphi*costheta*xkll+cosphi*ykll+sinphi*sintheta*zkll;
    zklll=-sintheta*xkll+costheta*zkll;


    dn=pow(betaq+band.mell/band.meld*(xkl-xklll)*(xkl-xklll)+band.melt/band.meld*
        ((ykl-yklll)*(ykl-yklll)+(zkl-zklll)*(zkl-zklll)),-2);

    dnmax=pow(betaq+band.melt/band.meld* (xkl-xklll)*(xkl-xklll)
        +band.melt/band.meld*((ykl-yklll)*(ykl-yklll)+(zkl-zklll)*(zkl-zklll)),-2);

    if(dn>dnmax)
    {
        cout<<"error EleImpScUpdateK: AcRe - Error";exit(0);
    }
    if(dn>Random()*dnmax)
    {
        //inverse BH-transformation
        xkl=xklll*sqrt(band.mell/band.meld);
        ykl=yklll*sqrt(band.melt/band.meld);
        zkl=zklll*sqrt(band.melt/band.meld);

        //move origin of the coordinate frame back into the Gamma point
       
        xkl=xkl+0.85*a0pi;

        //transform vector back into the BZ
        out_wedge(kx,ky,kz,xkl,ykl,zkl);

        //find the new symmetry operation and tetrahedron index
        band.CONFBZ(kx,ky,kz);
        
        Getisym();

        //cal. new energy (since the assumed elliptical parabolic band strucutre
        //only approximates the real band strucutre)
        //GETGVE(band);
        GetEnergy();

        //real scattering event
        Flag_SelfScatter=false;

        //statistic of scattering-process
    }
    else
    {
        //self scattering event
        Flag_SelfScatter=true;
    }
}
