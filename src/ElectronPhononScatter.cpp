
#include "mcmodel.h"


void MeshQuantities::ElectronPhononScatter()
{
    //Purpose:   performs electron scattering

    //local variables
    int iscat,ibold;
    double drand,dscat;
    double sca[MNScaEle*NBE],scasum;

    //only electrons are allowed in this routine
    if(par_type!=PELEC)
    {

    }
    //save inital band
    ibold=iband;

    //scattering statistics
        
    //first check for fictious scattering
    scasum=band.CALSCATTSUM(energy,iband);

    if(scasum>band.gamtet[itet])
    {
        drand=Random()*scasum;
    }
    else
    {
        drand=Random()*band.gamtet[itet];
    }

    if(scasum<drand)
    {
        Flag_SelfScatter=true;
    }
    else
    {
        Flag_SelfScatter=false;

        //calculate scattering-rate for all processes at energy energy
        band.CALSCATTE(energy,sca,iband);

        //select the scattering process with AcRejection
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
}


void MeshQuantities::EPSCAT(int iscat)
{
    //Purpose: Calculation of the electron state just after scattering
    //Parameter :  iscat = number of actual scattering-prozess

    //local variables
    int ibold,itetold,isymold;
    double eeold,xkold,ykold,zkold;

    //save particle state
    eeold=energy;
    xkold=kx;
    ykold=ky;
    zkold=kz;
    ibold=iband;
    itetold=itet;
    isymold=isym;

    //calculate final band and type of scattering
    iband=iscat/band.scpre+band.bandof[par_type];

    int ntemp=(iscat%band.scpre);

    if(!Flag_QuantumCorrection)
		switch(ntemp)
		{
									 //elastic intravalley
			case 0:energy=energy;TALIND(1);break;
									 //transversal acoustic g-Phonon
			case 1:energy=energy+band.temptag;TALIND(2);break;
									 //transversal acoustic g-Phonon
			case 2:energy=energy-band.temptag;TALIND(2);break;
									 //longitudinal acoustic g-Phonon
			case 3:energy=energy+band.templag;TALIND(2);break;
									 //longitudinal acoustic g-Phonon
			case 4:energy=energy-band.templag;TALIND(2);break;
									 //longitudinal optical g-Phonon
			case 5:energy=energy+band.templog;TALIND(2);break;
									 //longitudinal optical g-Phonon
			case 6:energy=energy-band.templog;TALIND(2);break;
									 //transversal acoustic f-Phonon
			case 7:energy=energy+band.temptaf;TALIND(3);break;
									 //transversal acoustic f-Phonon
			case 8:energy=energy-band.temptaf;TALIND(3);break;
									 //longitudinal acoustic f-Phonon
			case 9:energy=energy+band.templaf;TALIND(3);break;
									 //longitudinal acoustic f-Phonon
			case 10:energy=energy-band.templaf;TALIND(3);break;
									 //transversal optical f-Phonon
			case 11:energy=energy+band.temptof;TALIND(3);break;
									 //transversal optical f-Phonon
			case 12:energy=energy-band.temptof;TALIND(3);break;
			case 13:exit(0);break;   //Impact Ionization
			default:break;
		}
	else
		switch(ntemp){
			  case 0:energy=energy;TALIND(1);break;//elastic intravalley
			  case 1:energy=energy+band.temptag;TALIND(2);break;//transversal acoustic g-Phonon
			  case 2:energy=energy-band.temptag;TALIND(2);break;//transversal acoustic g-Phonon
			  case 3:energy=energy+band.templag;TALIND(2);break;//longitudinal acoustic g-Phonon
			  case 4:energy=energy-band.templag;TALIND(2);break;//longitudinal acoustic g-Phonon
			  case 5:energy=FEKLOAB(energy);TALIND(2);break;//longitudinal optical g-Phonon
			  case 6:energy=FEKLOEM(energy);TALIND(2);break;//longitudinal optical g-Phonon
			  case 7:energy=energy+band.temptaf;TALIND(3);break;//transversal acoustic f-Phonon
			  case 8:energy=energy-band.temptaf;TALIND(3);break;//transversal acoustic f-Phonon
			  case 9:energy=FEKLAAB(energy);TALIND(3);break;//longitudinal acoustic f-Phonon
			  case 10:energy=FEKLAEM(energy);TALIND(3);break;//longitudinal acoustic f-Phonon
			  case 11:energy=FEKTOAB(energy);TALIND(3);break;//transversal optical f-Phonon
			  case 12:energy=FEKTOEM(energy);TALIND(3);break;//transversal optical f-Phonon
			  case 13:exit(0);break;//Impact Ionization
			  default:break;
		}
	

    //calculate the final state for energy
    //energy and bandindex iband on an aeuquienergy
    //surface

    get_state();
    
    if(Flag_SelfScatter)
    {
        energy=eeold;
        kx=xkold;
        ky=ykold;
        kz=zkold;
        iband=ibold;
        itet=itetold;
        isym=isymold;
    }

   if (flag_heat) 
   {
     heat_to_point((energy - eeold) * fabs(charge), p_electron_heat);
//     heat_to_point(Rho, p_e_heat_weight);
   }
}


void MeshQuantities::EFPSCAT(int iscat)
{

    //Purpose: Calculation of the electron state just after scattering
    //	  Fischetti's phonon system
    //Parameter :  iscat = number of actual scattering-prozess

    //local variables
    bool aefl;
    int ibold,itetold,isymold;
    double eeold,xkold,ykold,zkold;
    double qq,qact,eep,xkctl,ykctl,zkctl;
    double koef,koefmax,aovtet;
    double xqct,yqct,zqct,eepmax;

    //save particle state
    eeold=energy;
    xkold=kx;
    ykold=ky;
    zkold=kz;
    ibold=iband;
    itetold=itet;
    isymold=isym;

    Flag_SelfScatter=false;

    iband=iscat/band.scpre+1+band.bandof[par_type];

    int cntemp=(iscat%band.scpre);
    switch(cntemp)
    {
        case 0:energy=energy+band.efopee;
        isym=(int(48*Random())%48);
        get_state();
        if(!Flag_SelfScatter)
        {
            xqct=kx-xkold;
            yqct=ky-ykold;
            zqct=kz-zkold;
            //confine the phonon vector to the first BZ
            band.CONFBZ(xqct,yqct,zqct);
            qq=sqrt(xqct*xqct+yqct*yqct+zqct*zqct);
            if(Random()>EOLINT(qq))
            {
                Flag_SelfScatter=true;
            }
        }
        break;
        case 1:energy=energy-band.efopee;
                                 //+1;mini=0
        isym=(int(48*Random())%48);
        get_state();
        if(!Flag_SelfScatter)
        {
            xqct=kx-xkold;
            yqct=ky-ykold;
            zqct=kz-zkold;
            //confine the phonon vector to the first BZ
            band.CONFBZ(xqct,yqct,zqct);
            qq=sqrt(xqct*xqct+yqct*yqct+zqct*zqct);
            if(Random()>EOLINT(qq))
            {
                Flag_SelfScatter=true;
            }
        }
        break;
        case 2:case 3:
        case 4:case 5:
            switch(cntemp)
            {
                case 2:FPFSAB();
                //select maximum phonon energy
                eepmax=band.efapeet;
                //select absorption
                aefl=true;
                break;
                case 3:FPFSEM();
                eepmax=band.efapeet;
                aefl=false;
                break;
                case 4:FPFSAB();
                eepmax=band.efapeel;
                aefl=true;
                break;
                case 5:FPFSEM();
                eepmax=band.efapeel;
                aefl=false;
                break;
                default:break;
            }
            if(!Flag_SelfScatter)
            {
                //select wedge
                                 //+1;minim=0;
                isym=(int(48*Random())%48);
                //transform center of chosen tetrahedron into wedge
                out_wedge(xkctl,ykctl,zkctl,band.xkct[itet],band.ykct[itet],band.zkct[itet]);
                //calculata momentum transfer (approximation)
                xqct=xkctl-xkold;
                yqct=ykctl-ykold;
                zqct=zkctl-zkold;
                //confine the phonon vector to the first BZ
                band.CONFBZ(xqct,yqct,zqct);
                //calculate phonon wave vector
                qact=sqrt(xqct*xqct+yqct*yqct+zqct*zqct)*sia0;
                //calculate phonon energy
                if(qact>TWOPI)
                {
                    eep=eepmax;
                }
                else
                {
                    if(qact<0.1)
                    {
                        eep=eepmax*0.176776695*qact;
                    }
                    else
                    {
                        eep=eepmax*sqrt(1.0-cos(0.25*qact));
                    }
                }
                //calculate final particle energy (absorption/emisson)
                if(aefl)
                {
                    energy=energy+eep;
                }
                else
                {
                    energy=energy-eep;
                }
                //tetrahedra which do not contain the energy are rejected
                if((band.eek[band.tet[0][itet]]<energy)||(band.eek[band.tet[3][itet]]>energy))
                {
                    //calculate the DOS within the selected tetrahedron (area
                    //over velocity) for the final energy
                    aovtet=band.FASTSURF(energy,itet)/band.vgt[itet];
                    //check the upper bound of the DOS
                    if(aovtet>band.maxaovtet[itet])
                    {
                    }
                    //selfscattering due to the upper bound of DOS
                    if(aovtet>Random()*band.maxaovtet[itet])
                    {
                        //calculate the coupling factor and the upper bound
                        if(aefl)
                        {
                            koefmax=35.0/(eepmax*eepmax);
                            koef=qact*qact/eep/(exp(eep)-1.0);
                        }
                        else
                        {
                            koefmax=6.0*PI*PI/eepmax*(1.0/(exp(eepmax)-1.0)+1.0);
                            koef=qact*qact/eep*(1.0/(exp(eep)-1.0)+1.0);
                        }
                        //check the upper bound
                        if(koef>koefmax)
                        {

                        }
                        //selfscattering to account for the upper bound of coupling factor
                        if(Random()*koefmax>koef)
                        {
                            Flag_SelfScatter = true;
                        }
                        else
                        {
                            //calculate final kcell-vector in the tetrahedron
                            getK();
                            //selscattering to account for the overlap integral
                            if(Random()>EOLINT(qact/sia0))
                            {
                                Flag_SelfScatter = true;
                            }
                        }
                    }
                    else
                    {
                        Flag_SelfScatter = true;
                    }
                }
                else
                {
                    Flag_SelfScatter = true;
                }
            }
            break;
        case 6:
            cout<<"error EFPSCAT: No II in EFPSCAT";exit(0);
            break;
        default:break;
    }
    //If selscattering restore inital particle state

    if(Flag_SelfScatter)
    {
        energy=eeold;
        kx=xkold;
        ky=ykold;
        kz=zkold;
        iband=ibold;
        itet=itetold;
        isym=isymold;
    }
    //GETGVE(band);

}


void MeshQuantities::ESCATII(int iscat)
{

    //Purpose: Calculation of the electron state just after scattering
    //=======  for II

    //local variables
    int ibold,itetold,isymold;
    double eeold,xkold,ykold,zkold;

    //save particle state
    eeold=energy;
    xkold=kx;
    ykold=ky;
    zkold=kz;
    ibold=iband;
    itetold=itet;
    isymold=isym;

    iband=iscat / band.scpre + band.bandof[par_type];

    //Impact Ionization
    //calculate the final state for energy
    //energy and bandindex iband on an aeuquienergy
    //surface
    energy = (energy-sieg) / 3.0;
    isym = (int(48*Random())%48);
    get_state();

    if(Flag_SelfScatter)
    {
        energy = eeold;
        kx = xkold;
        ky = ykold;
        kz = zkold;
        iband = ibold;
        itet = itetold;
        isym = isymold;
    }
    else
    {
        if(band.seciifl)
        {//band.seciifl 
            //generate secondary electron
          Particle tmp_par;

          list<Particle>::iterator
            tmp_iter = current_par_list->insert(current_par_list->end(), tmp_par);
          
          OutPar(tmp_iter);

          tmp_iter->par_id = -2;
          
          gen_par ++;
          
          //generate secondary hole
          par_type = PHOLE;
          
          iband=band.bandof[par_type] - 1;
  
          Flag_SelfScatter=true;
  
          while(Flag_SelfScatter)
            {
              iband=iband+1;
              if((iband-band.bandof[par_type])>=band.nband[par_type])
                {
                  
                }
              isym=(int(48*Random())%48);
              get_state();
            }
          
          charge = - charge;

          tmp_iter = current_par_list->insert(current_par_list->end(), tmp_par);
          
          OutPar(tmp_iter);

          tmp_iter->par_id = -3;
          gen_par ++;
          //get state of primary electron
          InPar(par_iter);
          //primary electron has negative kcell-vector of secondary electron
          kx=-kx;
          ky=-ky;
          kz=-kz;
          isym=band.indmat[-band.matsym[3][isym]+1][-band.matsym[4][isym]+1][-band.matsym[5][isym]+1]
            [band.matsym[0][isym]-1][band.matsym[1][isym]-1][band.matsym[2][isym]-1];
        }
    }
}

 void MeshQuantities::FPFSAB()
 {

    //Purpose: - calculate the final state for energy and iband on
    //    an aequienergy surface for Fischetti
    //    acoustic phonon scattering (absorbtion)

    //local variables
    int itab,nted,count2;
    double maxdostet;

    //energy index of tetraheder list
    itab=int(Max(0.0,energy)*band.dlistfp);
    if(itab>=MWLE)
    {
        itab=MWLEFP-1;           //maxim=MWLEFP-1
    }

    //get maximum dos of tetraheda
    //for the given energy
    maxdostet=band.maxaovfpab[iband][itab];

    //Number of tetraheder in list
    nted=band.ntlistfpab[itab][iband];
    count2=0;
    if(nted==0)
    {
        //no final state possible for this band and energy

        Flag_SelfScatter=true;

    }
    else
    {
        //chose random tetraheder from list
        itet=band.tlistfpab[int(nted*Random())+band.ptlistfpab[itab][iband]];
        while(((Random()*maxdostet)>band.maxaovtet[itet])&&(count2<10000))
        {
            itet=band.tlistfpab[int(nted*Random())+band.ptlistfpab[itab][iband]];

            count2=count2+1;

            //acception/rejection due to dos in tetrahedron

        }
        //check band index
        if(iband!=band.ibt[itet])
        {

            cout<<"error FPFSAB: iband <> band.ibt";exit(0);
        }
        //no final state found
        if(count2>=10000)
        {

            cout<<"error FPFSAB: count2 > 10000 iband = ";exit(0);
        }
    }

}


void MeshQuantities::FPFSEM()
{

    //Purpose: - calculate the final state for energy and iband on
    //    an aequienergy surface for Fischetti
    //    acoustic phonon scattering (emisson)

    //local variables
    int itab,nted,count2;
    double maxdostet;

    //energy index of tetraheder list
    itab=int(Max(0.0,energy)*band.dlistfp);
    if(itab>=MWLE)
    {
        itab=MWLEFP-1;           //maxim=MWLEFP-1
    }
    //get maximum dos of tetraheda
    //for the given energy
    maxdostet=band.maxaovfpem[iband][itab];

    //Number of tetraheder in list
    nted=band.ntlistfpem[itab][iband];
    count2=0;
    if(nted==0)
    {
        //no final state possible for this band and energy

        Flag_SelfScatter=true;

    }
    else
    {
        //chose random tetraheder from list
        itet= band.tlistfpem[int(nted*Random())+band.ptlistfpem[itab][iband]];
        while(((Random()*maxdostet)>band.maxaovtet[itet])&&(count2<10000))
        {

            itet= band.tlistfpem[int(nted*Random())+band.ptlistfpem[itab][iband]];
            count2=count2+1;
        }
        //acception/rejection due to dos in tetrahedron

        //check band index
        if(iband!=band.ibt[itet])
        {

            cout<<"error FPFSEM: iband <> band.ibt";exit(0);
        }
        //no final state found
        if(count2>=10000)
        {

            cout<<"error FPFSEM: count2 > 10000 iband = ";exit(0);
        }
    }
}


inline double MeshQuantities::EOLINT(double qq)
{
    //Purpose: Calculation of the electron overlap integral

    double qa;
    qa=qq*sia0*0.3907963210;
    if(qa>0.01)
        return pow((3.0*(sin(qa)-qa*cos(qa))/pow(qa,3.0)),2.0);
    return 1.0;
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double MeshQuantities::FEKLOEM(double eet)
{

//	local variables
	int i,j,k,it,ei;
	double ie,edeta,ran,ksum,efinal;
  {
	ie=eet/band.templog;
	if(ie<1.0e0)
	{
		efinal=1.0e-10;
		return efinal;
	}
	if(ie<3.05e0)
	{
		it=int((ie-1.0e0)*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=1e0;
	}
	else if(ie<5.625e0)
	{
		it=int((ie-3.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=2e0;
	}
	else if(ie<13.25e0)
	{
		it=int((ie-5.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=30+it;
		edeta=4e0;
	}
	else if(ie<28.5e0)
	{
		it=int(ie-13e0+0.5e0);
		if(it<1)	it=1;
		ei=45+it;
		edeta=8e0;
	}
	else
	{
		it=int((ie-28.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=60;
		ei=60+it;
		edeta=12e0;
	}
	ran=Random();
	i=0;
	ksum=0.0e0;
	while(ksum<ran)
	{
		i++;
		if(i==201)
		{
			ksum=1.1e0;
			i=1;
		}
		ksum+=kloem[i-1][ei-1];
	}
	j=i%2;
	k=int(i/2);
	if(j==1)	efinal=eet-band.templog*(1.0e0-0.01*edeta*k);
	else	efinal=eet-band.templog*(1.0e0+0.01*edeta*k);
  }
  
	return efinal;
}

double MeshQuantities::FEKLOAB(double eet)
{
//	local variables
	int i,j,k,it,ei;
	double ie,edeta,ran,ksum,efinal;
	it=0;
  
  {
	ie=eet/band.templog;
	if(ie<1.05e0)
	{
		it=int(ie * 10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=1e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=2e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=4e0;
	}
	else if(ie<26.5e0)
	{
		it=int((ie-11e0)+0.5e0);
		if(it<1)	it=1;
		ei=35+it;
		edeta=8e0;
	}
	else
	{
		it=int((ie-26.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=60;
		ei=50+it;
		edeta=12e0;
	}
	ran=Random();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=kloab[i-1][ei-1];
		if(i==201)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=i%2;
	k=int(i/2);
	if(j==1)	efinal=eet+band.templog*(1.0e0+0.01*edeta*k);
	else	efinal=eet+band.templog*(1.0e0-0.01*edeta*k);
  }
    
	return efinal;
}

double MeshQuantities::FEKLAAB(double eet)
{
//	local variables
	int i,j,k,it,ei;
	double ie,edeta,ran,ksum,efinal;
	it=0;
  
  {
	ie=eet/band.templaf;
	if(ie<1.05e0)
	{
		it=int(ie*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-11.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=115;
		ei=35+it;
		edeta=1.0e0;
	}
	ran=Random();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=klaab[i-1][ei-1];
		if(i==51)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=i%2;
	k=int(i/2);
	if(j==1)	efinal=eet+band.templaf*(1.0e0+0.04*edeta*k);
	else	efinal=eet+band.templaf*(1.0e0-0.04*edeta*k);
	if(efinal<=0e0)	efinal=(1.0e-5)*band.templaf;
  }
    
	return efinal;
}

double MeshQuantities::FEKLAEM(double eet)
{
//	local variables
	int i,j,k,it,ei;
	double ie,edeta,ran,ksum,efinal;
	it=0;
  
  {
	ie=eet/band.templaf;
	if(ie<1e0)
	{
		efinal=-1e-10;
		return efinal;
	}
	if(ie<3.05e0)
	{
		it=int((ie-1)*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<5.625e0)
	{
		it=int((ie-3.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else if(ie<13.25e0)
	{
		it=int((ie-5.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=30+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-13.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>115)	it=115;
		ei=45+it;
		edeta=1.0e0;
	}
	ran=Random();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=klaem[i-1][ei-1];
		if(i==51)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=i%2;
	k=int(i/2);
	if(j==1)	efinal=eet-band.templaf*(1.0e0 - 0.04*edeta*k);
	else	efinal=eet-band.templaf*(1.0e0 + 0.04*edeta*k);
	if(efinal<=0e0)	efinal=(1.0e-5)*band.templaf;
  }
  
  
	return efinal;
}

double MeshQuantities::FEKTOAB(double eet)
{
//	local variables
	int i,j,k,it,ei;
	double ie,edeta,ran,ksum,efinal;
	it=0;

  {
	ie=eet/band.temptof;
	if(ie<1.05e0)
	{
		it=int(ie*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<3.625e0)
	{
		it=int((ie-1.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=10+it;
		edeta=0.25e0;
	}
	else if(ie<11.25e0)
	{
		it=int((ie-3.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-11.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>60)	it=115;
		ei=35+it;
		edeta=1.0e0;
	}
	ran=Random();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=ktoab[i-1][ei-1];
		if(i==51)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=i%2;
	k=int(i/2);
	if(j==1)	efinal=eet+band.temptof*(1.0e0+0.04*edeta*k);
	else	efinal=eet+band.temptof*(1.0e0-0.04*edeta*k);
  }
  
  
	return efinal;
}

double MeshQuantities::FEKTOEM(double eet)
{
//	local variables
	int i,j,k,it,ei;
	double ie,edeta,ran,ksum,efinal;

  {
	ie=eet/band.temptof;
	if(ie<1.0e0)
	{
		efinal=1.0e-10;
		return efinal;
	}
	else if(ie<3.05e0)
	{
		it=int((ie-1.0e0)*10e0+0.5e0);
		if(it<1)	it=1;
		ei=it;
		edeta=0.25e0;
	}
	else if(ie<5.625e0)
	{
		it=int((ie-3.0e0)*4e0+0.5e0);
		if(it<1)	it=1;
		ei=20+it;
		edeta=0.5e0;
	}
	else if(ie<13.25e0)
	{
		it=int((ie-5.5e0)*2e0+0.5e0);
		if(it<1)	it=1;
		ei=30+it;
		edeta=0.5e0;
	}
	else
	{
		it=int((ie-13.0e0)+0.5e0);
		if(it<1)	it=1;
		if(it>115)	it=115;
		ei=45+it;
		edeta=1.0e0;
	}
	ran=Random();
	i=0;
	ksum=0e0;
	while(ksum<ran)
	{
		i++;
		ksum+=ktoem[i-1][ei-1];
		if(i==51)
		{
			ksum=1.1e0;
			i=1;
		}
	}
	j=i%2;
	k=int(i/2);
	if(j==1)	efinal=eet-band.temptof*(1.0e0-0.04*edeta*k);
	else	efinal=eet-band.temptof*(1.0e0+0.04*edeta*k);
  }
    
	return efinal;
}
