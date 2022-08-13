
#include "mcmodel.h"

double MeshQuantities::TetTime()
{
    //purpose:calculate time when the particle reaches the boarder of tetrahedron
    //
    double  xfl,yfl,zfl,xkl,ykl,zkl,tfp1,tfp2,tfp3,tfp4,nf1,nf2,nf3,nf4,tettime;

    //transform kcell-vector into irreduzible wedge
    in_wedge(xkl,ykl,zkl,kx,ky,kz);

    //transform electric field into irreduzible wedge
    in_wedge(xfl,yfl,zfl,Ex,Ey,Ez);//E
    //datantlin[a][b][c]ctetbtet(,,)
    //b=0~3   a=0/1/2<=>x/y/x a=3ctetbk
    nf1=band.datantlin[0][0][itet]*xfl+band.datantlin[1][0][itet] * yfl
        +band.datantlin[2][0][itet]*zfl;
    nf2=band.datantlin[0][1][itet]*xfl+band.datantlin[1][1][itet] * yfl
        +band.datantlin[2][1][itet]*zfl;
    nf3=band.datantlin[0][2][itet]*xfl+band.datantlin[1][2][itet] * yfl
        +band.datantlin[2][2][itet]*zfl;
    nf4=band.datantlin[0][3][itet]*xfl+band.datantlin[1][3][itet] * yfl
        +band.datantlin[2][3][itet]*zfl;

	//tettime<0 kitet() 
	//kitetEnf>0 tettime>0
	//,,
	//tettime<0tettime0
    if(nf1<0)
	{	
        tfp1=(band.datantlin[3][0][itet]-(band.datantlin[0][0][itet]*xkl
            +band.datantlin[1][0][itet]*ykl+band.datantlin[2][0][itet]*zkl))/nf1;
		if(tfp1<0){
		  cout << "tfp1 = " << tfp1 << ' ' << tfp1 * nf1 << ' ' << nf1 << endl;
		  tfp1=0;
		}
	}
    else
        tfp1=1e99;

    if(nf2<0)
	{
        tfp2=(band.datantlin[3][1][itet]-(band.datantlin[0][1][itet]*xkl
            +band.datantlin[1][1][itet]*ykl+band.datantlin[2][1][itet]*zkl))/nf2;
		if(tfp2<0){
		  cout << "tfp2 = " << tfp2 << ' ' << tfp2 * nf2 << ' ' << nf2 << endl;
		  tfp2=0;
		}
	}
    else
        tfp2=1e99;

    if(nf3<0)
	{
        tfp3=(band.datantlin[3][2][itet]-(band.datantlin[0][2][itet]*xkl
            +band.datantlin[1][2][itet]*ykl+band.datantlin[2][2][itet]*zkl))/nf3;
		if(tfp3<0){
		  cout << "tfp3 = " << tfp3 << ' ' << tfp3 * nf3 << ' ' << nf3 << endl;
		  tfp3=0;
		}
	}
    else
        tfp3=1e99;

    if(nf4<0)
	{
        tfp4=(band.datantlin[3][3][itet]-(band.datantlin[0][3][itet]*xkl
			+band.datantlin[1][3][itet]*ykl+band.datantlin[2][3][itet]*zkl))/nf4;
		if(tfp4<0){
		  cout << "tfp4 = " << tfp4 << ' ' << tfp4 * nf4 << ' ' << nf4 << endl;
		  tfp4=0;
		}
	}
    else
        tfp4=1e99;

    tettime=tfp1;
    itetdir=0;

    if(tfp2<tettime)
    {
        itetdir=1;
        tettime=tfp2;
    }

    if(tfp3<tettime)
    {
        itetdir=2;
        tettime=tfp3;
    }

    if(tfp4<tettime)
    {
        itetdir=3;
        tettime=tfp4;
    }
    if (tettime < MY_ZERO) {
      cout << "tettime = " << tettime << endl;
      band.output_tet(itet);
    }

    return tettime;
}


void MeshQuantities::HitTet()
{
    //Perform change of tetrahedron

    //tunnel in oxide

    //only change of tetrahedron

    if(band.NeibrTet[itetdir][itet]>=0)
        itet=band.NeibrTet[itetdir][itet];
    else if((band.NeibrTet[itetdir][itet]==-5)||(band.NeibrTet[itetdir][itet]==-6))
    {
        //34 
        //change of tetrahedron, wedge and BZ
        int ii = -band.NeibrTet[itetdir][itet]-2;
        kx = kx + band.xbz[ii][isym];//[5] !!!
        ky = ky + band.ybz[ii][isym];//ybzkisymiiykcell
        kz = kz + band.zbz[ii][isym];
                                 //newisymkisymiikcellisym  Band::SETNEIBSYM
        isym=band.newisym[ii][isym];
        GetisymBZ();              //kitet
        if (Flag_Catch) 
          return ;

        //       ..linear interpolation destroys symmetry and energy must be
        //         calculated anew
        GetEnergy();             // 
    }
    else
        //change of tetrahedron and wedge
        //125  wedge
        isym=band.newisym[-band.NeibrTet[itetdir][itet]-2][isym];
		//sym itet
}
