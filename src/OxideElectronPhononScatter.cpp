
#include "MCModel.h"

void MCModel::OxideElectronPhononScatter()
{
    //Purpose:   performs oxide electron scattering

    //local variables
    int iscat,ibold;
    double drand,dscat;
    double sca[MNScaOxEle*NBOE],scasum;

    //only oxide electrons are allowed in this routine

    if(mc_ParType!=POXEL)
    {
		cout<<"error OESCTR: Wrong particle type";exit(0);
    }

    //scattering statistics

    mc_sttt.ntotp[mc_iband]++;

    mc_sttt.dtotp[mc_iband]+=fabs(mc_Charge)/mc_band.gamtet[mc_itet];
    //save band index
    ibold=mc_iband;

    //first check for fictious scattering
    scasum=mc_band.CALSCATTSUM(mc_Energy,mc_iband);

    if(scasum>mc_band.gamtet[mc_itet])
    {

        drand = Random()*scasum;
    }
    else
    {
        drand = Random()*mc_band.gamtet[mc_itet];
    }

    if(scasum<drand)
    {
        mc_Flag_SelfScatter=true;
    }
    else
    {
        mc_Flag_SelfScatter=false;

        //calculate scattering-rate for all prozesses at energy mc_Energy
        mc_band.CALSCATTOE(mc_Energy,sca,mc_iband);

        //select the scattering prozess with AcRejection
        iscat=-1;
        dscat=0.0;

        while(drand>dscat)
        {

            iscat++;
            if(iscat>=mc_band.scproe*mc_band.nband[POXEL])
            {

                mc_Flag_SelfScatter=true;
                break;
            }
            else
                dscat=dscat+sca[iscat];

        }
        //iscat is number of chosen scatt.pro.
        if(!mc_Flag_SelfScatter)OEPSCAT(iscat);
    }

    //statistic of scattering-process
    if(mc_Flag_SelfScatter)
    {
        mc_sttt.nslfp[ibold]++;
    }
    else
    {
        mc_sttt.nreap[ibold]++;
        mc_sttt.nsctypoe[iscat][ibold-mc_band.bandof[mc_ParType]]++;
    }
}


void MCModel::OEPSCAT(int iscat)
{

    //Calculation of the oxide electron state just after scattering
    //Parameter :  iscat = number of actual scattering-prozess

    //local variables

    int ibold, itetold, mc_isymold;

    double eeold, xkold, ykold, zkold;
    //save old particle state
    eeold=mc_Energy;
    xkold=mc_Kx;
    ykold=mc_Ky;
    zkold=mc_Kz;
    ibold=mc_iband;
    itetold=mc_itet;
    mc_isymold=mc_isym;

    //calculate new band and type of scattering
    mc_iband=iscat/mc_band.scproe+ mc_band.bandof[mc_ParType];

    //elastic intravalley

    int ntemp=(iscat%mc_band.scproe);
    switch(ntemp)
    {
                                 //elastic intravalley
        case 0:mc_Energy=mc_Energy;break;
                                 //absorption
        case 1:mc_Energy=mc_Energy+mc_band.tempoeop;break;
                                 //emission
        case 2:mc_Energy=mc_Energy-mc_band.tempoeop;break;
        default:break;
    }

    //calculate the final state for energy
    //mc_Energy and bandindex mc_iband on an aeuquienergy
    //surface

    mc_isym=(int(48*Random())%48);
    //CALL GetState
    GetState();
    if(mc_Flag_SelfScatter)
    {
        mc_Energy=eeold;
        mc_Kx=xkold;
        mc_Ky=ykold;
        mc_Kz=zkold;
        mc_iband=ibold;
        mc_itet=itetold;
        mc_isym=mc_isymold;
    }
    //GETGVE(band);

}
