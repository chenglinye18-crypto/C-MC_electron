
#include "mcmodel.h"

void MeshQuantities::HolePhononScatter()
{
    //Purpose:   performs hole scattering

    //local variables
    int iscat, ibold;
    double drand, dscat;
    double sca[MNScaHole*NBH], scasum;

    //only holes are allowed in this routine
    if(par_type != PHOLE){
        cout<<"error HSCTR: Wrong particle type";
        exit(0);
    }
    //save intial band
    ibold = iband;

    //scattering statistics
    //first check for fictious scattering
    scasum = band.CALSCATTSUM(energy, iband);
    if(scasum > band.gamtet[itet]){
        drand = Random() * scasum;
    }
    else{
        drand = Random() * band.gamtet[itet];
    }
    if(scasum < drand){
        Flag_SelfScatter = true;
    }
    else{
        Flag_SelfScatter = false;

        //calculate scattering-rate for all prozesses at energy energy
        band.CALSCATTH(energy,sca,iband);
        //select the scattering prozess with AcRejection
        iscat = -1;
        dscat = 0.0;
        while(drand > dscat){
            iscat = iscat + 1;
            if(iscat >= (band.scprh*band.nband[PHOLE])){
                Flag_SelfScatter =true;
                break;

            }
            dscat = dscat + sca[iscat];
        }

        //iscat is number of chosen scatt.pro.
        if(!Flag_SelfScatter){
            if(((iscat % band.scprh) + 1) == band.scprh){
                HSCATII(iscat);
            }
            else{
                HPSCAT(iscat);
            }
        }
            
    }

}


void MeshQuantities::HPSCAT(int iscat)
{

    //Purpose: Calculation of the hole state just after scattering
    //Parameter :  iscat = number of actual scattering-prozess

    //local variables
    int ibold,itetold,isymold;
    double eeold,xkold,ykold,zkold;

    //save inital particle state
    eeold=energy;
    xkold=kx;
    ykold=ky;
    zkold=kz;
    ibold=iband;
    itetold=itet;
    isymold=isym;

    //calculate final band and type of scattering
    iband=iscat / band.scprh + band.bandof[par_type];

    int ntemp=(iscat % band.scprh);
    switch(ntemp)
    {
                                 //elastic intravalley
        case 0:energy=energy;break;
                                 //absorption
        case 1:energy=energy+band.temphop;break;
                                 //emission
        case 2:energy=energy-band.temphop;break;
        default:break;           //Impact Ionization

    }

    //calculate the final state for energy
    //energy and bandindex iband on an aeuquienergy
    //surface

    isym=(int(48*Random())%48);
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
      heat_to_point((energy - eeold) * fabs(charge) , p_hole_heat);
      //heat_to_point(Rho , p_h_heat_weight);
    }
    //GETGVE(band);

}


void MeshQuantities::HSCATII(int iscat)
{
    //Purpose: Calculation of the hole state just after scattering
    // for Impact Ionization scatter
    //Parameter :  iscat = number of actual scattering-prozess

    //local variables
    int ibold,itetold,isymold;
    double eeold,xkold,ykold,zkold;

    //save inital particle state
    eeold=energy;
    xkold=kx;
    ykold=ky;
    zkold=kz;
    ibold=iband;
    itetold=itet;
    isymold=isym;

    //calculate final band and type of scattering
    iband=iscat/band.scprh+band.bandof[par_type];

    //Impact Ionization
    //calculate the final state for energy
    //energy and bandindex iband on an aeuquienergy surface
    energy=(energy-sieg)/3.0;

    isym=(int(48*Random())%48);
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
    else if(band.seciifl)
    {
      Particle tmp_par;
      
      list<Particle>::iterator
        tmp_iter = current_par_list->insert(current_par_list->end(), tmp_par);
        
        OutPar(tmp_iter);
        tmp_iter->par_id = -4;
        gen_par ++;
        //generate secondary electron
        par_type=PELEC;
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
        tmp_iter->par_id = -4;
        gen_par ++;
        //get state of primary hole

        InPar(par_iter);
        //primary hole has negative kcell-vector of secondary hole
        kx=-kx;
        ky=-ky;
        kz=-kz;
        isym=band.indmat[-band.matsym[3][isym]+1][-band.matsym[4][isym]+1][-band.matsym[5][isym]+1]
            [band.matsym[0][isym]-1][band.matsym[1][isym]-1][band.matsym[2][isym]-1];

    }
    //GETGVE(band);

}
