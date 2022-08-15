
#include "mcmodel.h"

void MeshQuantities::get_state()
{
  bool smallfl;
  int  itab,nted,count2;
  double dostet;
  double maxdostet;

  /*all bands have energies larger than 0eV
    energies smaller than 0eV are sometimes chosen due
    to discretization errors in the scattering rates*/
  if(energy <= 0.0){
    Flag_SelfScatter = true;
  }
  else
    {
      /*get maximum dos of tetraheda for the given energy
        energy index of tetraheder list*/
      maxdostet = band.CALDOSTETMAX(energy, iband);
      
      /*fine list below 120meV*/
      itab = (int)(energy * band.dlists);
      if(itab < MWLES)
        smallfl = true;
      else
        {
          /*otherwise coarse list*/
          itab = (int)(energy * band.dlist);
          if(itab >= MWLE)
            itab = MWLE-1;
          smallfl = false;
        }
      /*number of tetraheder in list*/
      if(smallfl)
        nted = band.ntlists[itab][iband];
      else
        nted = band.ntlist[itab][iband];
      
      count2=0;

      if(nted == 0)
        Flag_SelfScatter = true;
      else
        {
          bool nfindfl=true;
          while(nfindfl)
            {/*choose a tetraheder randomly */
              if(smallfl)
                itet=band.tlists[(int)(nted * Random()) + band.ptlists[itab][iband]];
              else
                itet=band.tlist[(int)(nted * Random()) + band.ptlist[itab][iband]];

              count2++;
              /*tetrahedra which do not contain the energy are rejected)
               for more detail see PhysicsDoc.pdf page26 */
              if((count2 < 10000)&&!((band.eek[band.tet[0][itet]] < energy)&&(band.eek[band.tet[3][itet]] > energy)))continue;
              /*dos of this tetrahedra */
              dostet=band.FASTSURF(energy, itet) / band.vgt[itet];
              /*acception/rejection due to dos in tetrahedron*/
              if(((Random() * maxdostet) > dostet) && (count2 < 10000))continue;
              nfindfl=false;
            }
          /*no final state found*/
          if(count2 >= 10000)
            Flag_SelfScatter = true;
          else{
              /*at this point a tetrahedron has been selected
                now we are looking for a point on the surface within the tetrahedron*/
            Flag_SelfScatter = false;
            getK();
          }
        }
    }
}

void MeshQuantities::getK() {
  double rv1,rv2,rfac,xkl,ykl,zkl;
  double ax,ay,az,bx,by,bz,cx,cy,cz;
  double dx,dy,dz,a1,a2,a3,a4,ar;

  /*get pair of random numbers, which is homogeneously distributed in
    a triangle */
  rv1=1.0-sqrt(Random());
  rv2=Random()*(1.0-rv1);

  if((energy < band.eek[band.tet[0][itet]]) || (energy > band.eek[band.tet[3][itet]]))
	{
      cout<<"make up new energy in GetK\n";
      energy = band.eek[band.tet[0][itet]] + (band.eek[band.tet[3][itet]]-band.eek[band.tet[0][itet]])*Random();
	}
  /*Case 1: equi-energy surface is a  triangular intersection with energy lower e2 */
  if(energy <= band.eek[band.tet[1][itet]])
    {
      /*cal. vectors (relative to origin of tetrahedron) pointing to the
        three nodes of the triangle */
      
      rfac=(energy-band.eek[band.tet[0][itet]]) / (band.eek[band.tet[1][itet]]-band.eek[band.tet[0][itet]]);
      ax=rfac*(band.xkk[band.tet[1][itet]] - band.xkk[band.tet[0][itet]]);
      ay=rfac*(band.ykk[band.tet[1][itet]] - band.ykk[band.tet[0][itet]]);
      az=rfac*(band.zkk[band.tet[1][itet]] - band.zkk[band.tet[0][itet]]);
      if(rfac<0.0||rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      rfac=(energy-band.eek[band.tet[0][itet]])/(band.eek[band.tet[2][itet]]-band.eek[band.tet[0][itet]]);
      bx=rfac*(band.xkk[band.tet[2][itet]]-band.xkk[band.tet[0][itet]]);
      by=rfac*(band.ykk[band.tet[2][itet]]-band.ykk[band.tet[0][itet]]);
      bz=rfac*(band.zkk[band.tet[2][itet]]-band.zkk[band.tet[0][itet]]);
      if(rfac<0.0||rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }

      rfac=(energy-band.eek[band.tet[0][itet]])/(band.eek[band.tet[3][itet]]-band.eek[band.tet[0][itet]]);
      cx=rfac*(band.xkk[band.tet[3][itet]]-band.xkk[band.tet[0][itet]]);
      cy=rfac*(band.ykk[band.tet[3][itet]]-band.ykk[band.tet[0][itet]]);
      cz=rfac*(band.zkk[band.tet[3][itet]]-band.zkk[band.tet[0][itet]]);
      if(rfac<0.0||rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      /*ok, choose a point on the equi-energy surface */
      xkl = rv1 * (bx-ax) + rv2*(cx-ax) + ax + band.xkk[band.tet[0][itet]];
      ykl = rv1 * (by-ay) + rv2*(cy-ay) + ay + band.ykk[band.tet[0][itet]];
      zkl = rv1 * (bz-az) + rv2*(cz-az) + az + band.zkk[band.tet[0][itet]];
    }
  /*Case 2 : equi-energy surface is  quadrangular with energy between e2 and e3 */
  else if(energy<band.eek[band.tet[2][itet]])
    {
      /*cal. vectors (relative to origin of kcell-space) pointing to the
        four nodes of the quadrangle*/
      rfac=(energy-band.eek[band.tet[0][itet]])/(band.eek[band.tet[2][itet]]-band.eek[band.tet[0][itet]]);
      ax=rfac*(band.xkk[band.tet[2][itet]]-band.xkk[band.tet[0][itet]])+band.xkk[band.tet[0][itet]];
      ay=rfac*(band.ykk[band.tet[2][itet]]-band.ykk[band.tet[0][itet]])+band.ykk[band.tet[0][itet]];
      az=rfac*(band.zkk[band.tet[2][itet]]-band.zkk[band.tet[0][itet]])+band.zkk[band.tet[0][itet]];
      if(rfac<0.0||rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      rfac=(energy-band.eek[band.tet[1][itet]])/(band.eek[band.tet[2][itet]]-band.eek[band.tet[1][itet]]);
      bx=rfac*(band.xkk[band.tet[2][itet]]-band.xkk[band.tet[1][itet]])+band.xkk[band.tet[1][itet]];
      by=rfac*(band.ykk[band.tet[2][itet]]-band.ykk[band.tet[1][itet]])+band.ykk[band.tet[1][itet]];
      bz=rfac*(band.zkk[band.tet[2][itet]]-band.zkk[band.tet[1][itet]])+band.zkk[band.tet[1][itet]];
      if(rfac<0.0||rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      rfac=(energy-band.eek[band.tet[0][itet]])/(band.eek[band.tet[3][itet]]-band.eek[band.tet[0][itet]]);
      cx=rfac*(band.xkk[band.tet[3][itet]]-band.xkk[band.tet[0][itet]])+band.xkk[band.tet[0][itet]];
      cy=rfac*(band.ykk[band.tet[3][itet]]-band.ykk[band.tet[0][itet]])+band.ykk[band.tet[0][itet]];
      cz=rfac*(band.zkk[band.tet[3][itet]]-band.zkk[band.tet[0][itet]])+band.zkk[band.tet[0][itet]];
      if(rfac<0.0||rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      rfac=(energy-band.eek[band.tet[1][itet]])/(band.eek[band.tet[3][itet]]-band.eek[band.tet[1][itet]]);
      dx=rfac*(band.xkk[band.tet[3][itet]]-band.xkk[band.tet[1][itet]])+band.xkk[band.tet[1][itet]];
      dy=rfac*(band.ykk[band.tet[3][itet]]-band.ykk[band.tet[1][itet]])+band.ykk[band.tet[1][itet]];
      dz=rfac*(band.zkk[band.tet[3][itet]]-band.zkk[band.tet[1][itet]])+band.zkk[band.tet[1][itet]];
      if(rfac<0.0||rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      /*cal. area of the four possible triangles made from the four nodes
        of the quadrangle
        TODO: why 4 triangles?
      */
      a1=band.FII(ax,ay,az,bx,by,bz,cx,cy,cz);
      a2=band.FII(dx,dy,dz,bx,by,bz,cx,cy,cz);
      a3=band.FII(ax,ay,az,dx,dy,dz,cx,cy,cz);
      a4=band.FII(ax,ay,az,bx,by,bz,dx,dy,dz);
      ar=Random()*(a1+a2+a3+a4);

      /*select one of the four triangles and calculate kcell-vector */
      if(ar<a1)
        {
          xkl=rv1*(bx-ax)+rv2*(cx-ax)+ax;
          ykl=rv1*(by-ay)+rv2*(cy-ay)+ay;
          zkl=rv1*(bz-az)+rv2*(cz-az)+az;
        }
      else if(ar<(a1+a2))
        {
          xkl=rv1*(bx-dx)+rv2*(cx-dx)+dx;
          ykl=rv1*(by-dy)+rv2*(cy-dy)+dy;
          zkl=rv1*(bz-dz)+rv2*(cz-dz)+dz;
        }
      else if(ar<(a1+a2+a3))
        {
          xkl=rv1*(dx-ax)+rv2*(cx-ax)+ax;
          ykl=rv1*(dy-ay)+rv2*(cy-ay)+ay;
          zkl=rv1*(dz-az)+rv2*(cz-az)+az;
        }
      else
        {
          xkl=rv1*(bx-ax)+rv2*(dx-ax)+ax;
          ykl=rv1*(by-ay)+rv2*(dy-ay)+ay;
          zkl=rv1*(bz-az)+rv2*(dz-az)+az;
        }
    }
  else
    /* Case 3 : equi-energy surface is a  triangular  with energy higher e3 */
    {
      /*cal. vectors (relative to origin of tetrahedron) pointing to the
        three nodes of the triangle*/
      rfac=(energy-band.eek[band.tet[3][itet]])/(band.eek[band.tet[0][itet]]-band.eek[band.tet[3][itet]]);
      ax=rfac*(band.xkk[band.tet[0][itet]]-band.xkk[band.tet[3][itet]]);
      ay=rfac*(band.ykk[band.tet[0][itet]]-band.ykk[band.tet[3][itet]]);
      az=rfac*(band.zkk[band.tet[0][itet]]-band.zkk[band.tet[3][itet]]);
      if(rfac<0.0&&rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      rfac=(energy-band.eek[band.tet[3][itet]])/(band.eek[band.tet[1][itet]]-band.eek[band.tet[3][itet]]);
      bx=rfac*(band.xkk[band.tet[1][itet]]-band.xkk[band.tet[3][itet]]);
      by=rfac*(band.ykk[band.tet[1][itet]]-band.ykk[band.tet[3][itet]]);
      bz=rfac*(band.zkk[band.tet[1][itet]]-band.zkk[band.tet[3][itet]]);
      
      if(rfac<0.0&&rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      rfac=(energy - band.eek[band.tet[3][itet]])/(band.eek[band.tet[2][itet]]-band.eek[band.tet[3][itet]]);
      cx=rfac*(band.xkk[band.tet[2][itet]]-band.xkk[band.tet[3][itet]]);
      cy=rfac*(band.ykk[band.tet[2][itet]]-band.ykk[band.tet[3][itet]]);
      cz=rfac*(band.zkk[band.tet[2][itet]]-band.zkk[band.tet[3][itet]]);
      if(rfac<0.0&&rfac>1.0)
        {
          cout<<"GetK  rfac<0.0||rfac>1.0";
        }
      xkl=rv1*(bx-ax)+rv2*(cx-ax)+ax+band.xkk[band.tet[3][itet]];
      ykl=rv1*(by-ay)+rv2*(cy-ay)+ay+band.ykk[band.tet[3][itet]];
      zkl=rv1*(bz-az)+rv2*(cz-az)+az+band.zkk[band.tet[3][itet]];
    }
  /*end of finding the final state*/

  /*apply the symmetry transformation according to isym */
 out_wedge(kx, ky, kz, xkl, ykl, zkl);
}

void MeshQuantities::out_wedge(double &xout,double &yout,double &zout,
                               double xin,double yin,double zin)
{

  band.OutWedge(xout, yout, zout, xin, yin, zin, isym);
}

void MeshQuantities::in_wedge(double &xout,double &yout,double &zout,
                             double xin,double yin,double zin){
  
  double vsort[3];

  vsort[band.matsym[0][isym]-1] = xin * band.matsym[3][isym];
  vsort[band.matsym[1][isym]-1] = yin * band.matsym[4][isym];
  vsort[band.matsym[2][isym]-1] = zin * band.matsym[5][isym];

  xout = vsort[0];
  yout = vsort[1];
  zout = vsort[2];
}
                             
void MeshQuantities::GetV()
{

  /*mc_band.itet is within irreduzible wedge*/
  
  out_wedge(vx, vy, vz, band.vt[4 * itet], band.vt[4 * itet + 1], band.vt[4 * itet + 2]);
}

void MeshQuantities::GetEnergy()
{
  double dx,dy,dz;

  /*transform kcell-vector into irreducible wedge */
  in_wedge(dx,dy,dz,kx,ky,kz);
  
  /*momentum relative to origin of tetrahedron*/
  dx -= band.xkk[band.tet[0][itet]];
  dy -= band.ykk[band.tet[0][itet]];
  dz -= band.zkk[band.tet[0][itet]];
  
  /* new energy */
  energy = band.vt[4*itet+3] + band.vt[4*itet] * dx + band.vt[4*itet+1] * dy + band.vt[4 * itet + 2] * dz;
}

void MeshQuantities::GetisymBZ()
{
 
  int it,itl,ie,ibase;
  double d1,d2,d3,d4;
  double xkl,ykl,zkl;
  bool Flag_Find;

  /*transform into the irreduzibile wedge*/
  in_wedge(xkl,ykl,zkl,kx,ky,kz);
  Flag_Find=false;
  //check original tetrahedron
  if(iband ==band.ibt[itet])
    {
      ibase=itet;
      d1=band.datantlin[0][0][ibase]*xkl +band.datantlin[1][0][ibase]*ykl
        +band.datantlin[2][0][ibase]*zkl-band.datantlin[3][0][ibase];
      d2=band.datantlin[0][1][ibase]*xkl +band.datantlin[1][1][ibase]*ykl
        +band.datantlin[2][1][ibase]*zkl-band.datantlin[3][1][ibase];
      d3=band.datantlin[0][2][ibase]*xkl +band.datantlin[1][2][ibase]*ykl
        +band.datantlin[2][2][ibase]*zkl-band.datantlin[3][2][ibase];
      d4=band.datantlin[0][3][ibase]*xkl +band.datantlin[1][3][ibase]*ykl
        +band.datantlin[2][3][ibase]*zkl-band.datantlin[3][3][ibase];
      if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))Flag_Find=true;
    }
  /* search only for tetraheda containing the given energy Energy and
       which are part of the surface of the BZ
  */
  if(!Flag_Find)
    {
      ie=(int)(energy*band.dlist);
      if(ie<0)ie=0;
      if(ie>=MWLE)ie=MWLE-1;
      for(itl=band.pslist[ie][iband];
          itl<=band.pslist[ie][iband]+band.nslist[ie][iband]-1;itl++)
        {
          it = band.slist[itl];
          if(iband == band.ibt[it])
            {
              ibase=it;
              d1=band.datantlin[0][0][ibase]*xkl + band.datantlin[1][0][ibase]*ykl
                +band.datantlin[2][0][ibase]*zkl - band.datantlin[3][0][ibase];
              d2=band.datantlin[0][1][ibase]*xkl + band.datantlin[1][1][ibase]*ykl
                +band.datantlin[2][1][ibase]*zkl - band.datantlin[3][1][ibase];
              d3=band.datantlin[0][2][ibase]*xkl + band.datantlin[1][2][ibase]*ykl
                +band.datantlin[2][2][ibase]*zkl - band.datantlin[3][2][ibase];
              d4=band.datantlin[0][3][ibase]*xkl + band.datantlin[1][3][ibase] * ykl
                +band.datantlin[2][3][ibase]*zkl-band.datantlin[3][3][ibase];
              if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))
                {
                  Flag_Find=true;
                  itet=it;
                  break;
                }
            }
        }
    }

  /*_____this should not happen, but ...*/
  if(!Flag_Find)
    Getisym();
}


void MeshQuantities::Getisym()
{
  //purpose:   update symmetry operation and band.tet index

  int pos[3]={1,1,1},it,itl,ie,ibase,MinK,MaxK;
  int icx,icy,icz,ic;
  double d1,d2,d3,d4;
  double xkl,ykl,zkl;
  double TempMin,TempMax;
  bool Flag_Find;

  if(fabs(kx)>fabs(ky))
    {
      TempMin=fabs(ky);
      TempMax=fabs(kx);
      MinK=1;
      MaxK=0;
    }
  else
    {
      TempMin=fabs(kx);
      TempMax=fabs(ky);
      MinK=0;
      MaxK=1;
    }
  if(fabs(kz)<TempMin) MinK=2;
  else if(fabs(kz)>TempMax) MaxK=2;
  pos[MinK]=2;
  pos[MaxK]=0;              

  isym=band.indmat[sign(kx)+1][sign(ky)+1][sign(kz)+1][pos[0]][pos[1]][pos[2]];

  //find tetrahedron index
  in_wedge(xkl,ykl,zkl,kx,ky,kz);
  Flag_Find=false;
  //check original tetrahedron
  if(iband==band.ibt[itet]) 
    {
      ibase=itet;
      d1=band.datantlin[0][0][ibase]*xkl +band.datantlin[1][0][ibase]*ykl
        +band.datantlin[2][0][ibase]*zkl-band.datantlin[3][0][ibase];
      d2=band.datantlin[0][1][ibase]*xkl +band.datantlin[1][1][ibase]*ykl
        +band.datantlin[2][1][ibase]*zkl-band.datantlin[3][1][ibase];
      d3=band.datantlin[0][2][ibase]*xkl +band.datantlin[1][2][ibase]*ykl
        +band.datantlin[2][2][ibase]*zkl-band.datantlin[3][2][ibase];
      d4=band.datantlin[0][3][ibase]*xkl +band.datantlin[1][3][ibase]*ykl
        +band.datantlin[2][3][ibase]*zkl-band.datantlin[3][3][ibase];
      if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))
        Flag_Find=true;
    }

  //search only for tetraheda contained in the given cube
  //search fine list
  if((!Flag_Find)&&(iband==band.bandof[PELEC]))
    {
      icx=(int)((xkl/a0pi-0.7)/0.3 * MNCubeX);
      icy=(int)(ykl/(a0pi*0.075) * MNCubeY);
      icz=(int)(zkl/(a0pi*0.075) * MNCubeZ);
      if (icx>=0 && icx<MNCubeX && icy>=0 && icy<MNCubeY && icz>=0 && icz<MNCubeZ)
        {                                
          ic=icz+ MNCubeZ*(icy+MNCubeY*icx);
          for(itl=band.pclist[ic];itl<=band.pclist[ic]+band.nclist[ic]-1;itl++)
            {
              it=band.clist[itl];
              ibase=it;
              d1=band.datantlin[0][0][ibase]*xkl +band.datantlin[1][0][ibase]*ykl
                +band.datantlin[2][0][ibase]*zkl-band.datantlin[3][0][ibase];
              d2=band.datantlin[0][1][ibase]*xkl +band.datantlin[1][1][ibase]*ykl
                +band.datantlin[2][1][ibase]*zkl-band.datantlin[3][1][ibase];
              d3=band.datantlin[0][2][ibase]*xkl +band.datantlin[1][2][ibase]*ykl
                +band.datantlin[2][2][ibase]*zkl-band.datantlin[3][2][ibase];
              d4=band.datantlin[0][3][ibase]*xkl +band.datantlin[1][3][ibase]*ykl
                +band.datantlin[2][3][ibase]*zkl-band.datantlin[3][3][ibase];
              if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))
                {
                  Flag_Find=true;
                  itet=it;
                  break;
                }
            }
        }
    }

  //search only for tetraheda containing the given energy Energy
  //search fine list
  if(!Flag_Find)
    {
      ie=(int)(energy*band.dlists);
      if(ie<0)ie=0;
      if(ie<MWLES)
        for(itl=band.ptlists[ie][iband];itl<=band.ptlists[ie][iband]+band.ntlists[ie][iband]-1;itl++)
          {
            it=band.tlists[itl];

            if(iband==band.ibt[it])
              {
                ibase=it;
                d1=band.datantlin[0][0][ibase]*xkl +band.datantlin[1][0][ibase]*ykl
                  +band.datantlin[2][0][ibase]*zkl -band.datantlin[3][0][ibase];
                d2=band.datantlin[0][1][ibase]*xkl +band.datantlin[1][1][ibase]*ykl
                  +band.datantlin[2][1][ibase]*zkl -band.datantlin[3][1][ibase];
                d3=band.datantlin[0][2][ibase]*xkl +band.datantlin[1][2][ibase]*ykl
                  +band.datantlin[2][2][ibase]*zkl -band.datantlin[3][2][ibase];
                d4=band.datantlin[0][3][ibase]*xkl +band.datantlin[1][3][ibase]*ykl
                  +band.datantlin[2][3][ibase]*zkl -band.datantlin[3][3][ibase];
                if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))
                  {
                    Flag_Find=true;
                    itet=it;
                    break;
                  }
              }
          }
    }

  //search only for tetraheda containing the given energy Energy
  if(!Flag_Find)
    {
      ie=(int)(energy*band.dlist);
      if(ie<0)ie=0;
      if(ie>=MWLE)ie=MWLE-1;
      for(itl=band.ptlist[ie][iband];
          itl<=band.ptlist[ie][iband]+band.ntlist[ie][iband]-1;itl++)
        {
          it=band.tlist[itl];
          if(iband==band.ibt[it])
            {
              ibase=it;
              d1=band.datantlin[0][0][ibase]*xkl +band.datantlin[1][0][ibase]*ykl
                +band.datantlin[2][0][ibase]*zkl -band.datantlin[3][0][ibase];
              d2=band.datantlin[0][1][ibase]*xkl +band.datantlin[1][1][ibase]*ykl
                +band.datantlin[2][1][ibase]*zkl -band.datantlin[3][1][ibase];
              d3=band.datantlin[0][2][ibase]*xkl +band.datantlin[1][2][ibase]*ykl
                +band.datantlin[2][2][ibase]*zkl -band.datantlin[3][2][ibase];
              d4=band.datantlin[0][3][ibase]*xkl +band.datantlin[1][3][ibase]*ykl
                +band.datantlin[2][3][ibase]*zkl -band.datantlin[3][3][ibase];
              if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))
                {
                  Flag_Find=true;
                  itet=it;
                  break;
                }
            }
        }
    }

  //search only for tetraheda containing the given energy Energy + dee
  ie=ie+1;
  if((!Flag_Find)&&(ie<MWLE))
    for(itl=band.ptlist[ie][iband];
        itl<=band.ptlist[ie][iband]+band.ntlist[ie][iband]-1;itl++)
      {
        it=band.tlist[itl];
        if(iband==band.ibt[it])
          {
            ibase=it;
            d1=band.datantlin[0][0][ibase]*xkl +band.datantlin[1][0][ibase]*ykl
              +band.datantlin[2][0][ibase]*zkl -band.datantlin[3][0][ibase];
            d2=band.datantlin[0][1][ibase]*xkl +band.datantlin[1][1][ibase]*ykl
              +band.datantlin[2][1][ibase]*zkl -band.datantlin[3][1][ibase];
            d3=band.datantlin[0][2][ibase]*xkl +band.datantlin[1][2][ibase]*ykl
              +band.datantlin[2][2][ibase]*zkl -band.datantlin[3][2][ibase];
            d4=band.datantlin[0][3][ibase]*xkl +band.datantlin[1][3][ibase]*ykl
              +band.datantlin[2][3][ibase]*zkl -band.datantlin[3][3][ibase];
            if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))
              {
                Flag_Find=true;
                itet=it;
                break;
              }
          }
      }

  //search only for tetraheda containing the given energy Energy - dee
  ie=ie-2;
  if((!Flag_Find)&&(ie>=0))
    for(itl=band.ptlist[ie][iband];
        itl<=band.ptlist[ie][iband]+band.ntlist[ie][iband]-1;itl++)
      {
        it=band.tlist[itl];
        if(iband==band.ibt[it])
          {
            ibase=it;
            d1=band.datantlin[0][0][ibase]*xkl +band.datantlin[1][0][ibase]*ykl
              +band.datantlin[2][0][ibase]*zkl -band.datantlin[3][0][ibase];
            d2=band.datantlin[0][1][ibase]*xkl +band.datantlin[1][1][ibase]*ykl
              +band.datantlin[2][1][ibase]*zkl -band.datantlin[3][1][ibase];
            d3=band.datantlin[0][2][ibase]*xkl +band.datantlin[1][2][ibase]*ykl
              +band.datantlin[2][2][ibase]*zkl -band.datantlin[3][2][ibase];
            d4=band.datantlin[0][3][ibase]*xkl +band.datantlin[1][3][ibase]*ykl
              +band.datantlin[2][3][ibase]*zkl -band.datantlin[3][3][ibase];
            if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))
              {
                Flag_Find=true;
                itet=it;
                break;
              }
          }
      }

  //search all tetraheda, if not found above
  if(!Flag_Find)
    for(it=0;it<band.nt;it++)
      if(iband==band.ibt[it])
        {
          ibase=it;
          d1=band.datantlin[0][0][ibase]*xkl +band.datantlin[1][0][ibase]*ykl
            +band.datantlin[2][0][ibase]*zkl -band.datantlin[3][0][ibase];
          d2=band.datantlin[0][1][ibase]*xkl +band.datantlin[1][1][ibase]*ykl
            +band.datantlin[2][1][ibase]*zkl -band.datantlin[3][1][ibase];
          d3=band.datantlin[0][2][ibase]*xkl +band.datantlin[1][2][ibase]*ykl
            +band.datantlin[2][2][ibase]*zkl -band.datantlin[3][2][ibase];
          d4=band.datantlin[0][3][ibase]*xkl +band.datantlin[1][3][ibase]*ykl
            +band.datantlin[2][3][ibase]*zkl -band.datantlin[3][3][ibase];
          if((d1>=-MY_ZERO)&&(d2>=-MY_ZERO)&&(d3>=-MY_ZERO)&&(d4>=-MY_ZERO))
            {
              Flag_Find=true;
              itet=it;
              break;
            }
        }

  if(!Flag_Find)
    {
      cout  <<"----------------------------------------------------------"<<endl;
      cout <<"Getisym: no itet found"<<endl;
      cout << "Getisym: no itet found"<<endl;
      Flag_Catch = true;
//      exit(0);
    }
}
void MeshQuantities::TALIND(int ityp)
{

  int talend[5],ital;

  //number of valley
  if((fabs(kx)>=fabs(ky))&&(fabs(kx)>=fabs(kz)))
    {
      if(kx>=0.0)
        {
          //location valley  1
          if(ityp==1)
            {
              ital=1;
            }
          if(ityp==2)
            {
              ital=2;
            }
          if(ityp==3)
            {
              talend[1]=3;
              talend[2]=4;
              talend[3]=5;
              talend[4]=6;
              ital=talend[(int)(Random()*4.0+1.0)];
            }
        }
      if(kx<0.0)
        {
          //location valley  2
          if(ityp==1)
            {
              ital=2;
            }
          if(ityp==2)
            {
              ital=1;
            }
          if(ityp==3)
            {
              talend[1]=3;
              talend[2]=4;
              talend[3]=5;
              talend[4]=6;
              ital=talend[(int)(Random()*4.0+1.0)];
            }
        }
    }
  if((fabs(ky)>=fabs(kz))&&(fabs(ky)>=fabs(kx)))
    {
      if(ky>=0.0)
        {
          //location valley  3
          if(ityp==1)
            ital=3;
          if(ityp==2)
            ital=4;
          if(ityp==3)
            {
              talend[1]=1;
              talend[2]=2;
              talend[3]=5;
              talend[4]=6;
              ital=talend[(int)(Random()*4.0+1.0)];
            }
        }
      if(ky<0.0)
        {
          //location valley  4
          if(ityp==1)
            ital=4;
          if(ityp==2)
            ital=3;
          if(ityp==3)
            {
              talend[1]=1;
              talend[2]=2;
              talend[3]=5;
              talend[4]=6;
              ital=talend[(int)(Random()*4.0+1.0)];
            }
        }
    }

  if((fabs(kz)>fabs(ky))&&(fabs(kz)>fabs(kx)))
    {
      if(kz>=0.0)
        {
          //location valley  5
          if(ityp==1)
            ital=5;
          if(ityp==2)
            ital=6;
          if(ityp==3)
            {
              talend[1]=1;
              talend[2]=2;
              talend[3]=3;
              talend[4]=4;
              ital=talend[(int)(Random()*4.0+1.0)];
            }
        }
      if(kz<0.0)
        {
          //location valley  6
          if(ityp==1)
            ital=6;
          if(ityp==2)
            ital=5;
          if(ityp==3)
            {
              talend[1]=1;
              talend[2]=2;
              talend[3]=3;
              talend[4]=4;
              ital=talend[(int)(Random()*4.0+1.0)];
            }
        }
    }

  isym=(int)(Random()*8.0)+8*(ital-1);
}

void MeshQuantities::init_deep() {
  	int it,ik;
    int j;
    int size;
    ifstream ftpkloem;

    /*
     * k: const or k-space ?
     * lo, to, la, ta: 光学波和声学波
     * ab: absorption
     * em: emission
     */

    string filename;
    filename = bs_path + "/kloem.txt";
    ftpkloem.open(filename.c_str());

    for(it = 0; it < 120; it++)
      for(ik = 0; ik < 201; ik++)
        ftpkloem >> kloem[ik][it];
    
    ftpkloem.close();
      
    //c	   kloab
    ifstream ftpkloab;
    filename = bs_path + "/kloab.txt";
    ftpkloab.open(filename.c_str());
    
    for(it = 0; it < 110; it++)
      for(ik = 0;ik < 201; ik++)
        ftpkloab >> kloab[ik][it];
    ftpkloab.close();

    //c	   ktoem
    ifstream ftpktoem;
    filename = bs_path + "/ktoem.txt" ;
    ftpktoem.open(filename.c_str());
    for(it = 0; it < 160; it++)
      for(ik = 0;ik < 51; ik++)
        ftpktoem >> ktoem[ik][it];
    ftpktoem.close();

    //c	   ktoab
    ifstream ftpktoab;
    filename = bs_path + "/ktoab.txt";
    ftpktoab.open(filename.c_str());
    for(it = 0; it < 150; it++)
      for(ik = 0; ik < 51; ik++)
        ftpktoab >> ktoab[ik][it];
    ftpktoab.close();


    //c	   klaem
    ifstream ftpklaem;
    filename = bs_path + "/klaem.txt";
    ftpklaem.open(filename.c_str());
    for(it = 0; it < 160; it++)
      for(ik  =0; ik < 51; ik++)
        ftpklaem >> klaem[ik][it];
    ftpklaem.close();
      
    //c	   klaab
    ifstream ftpklaab;
    filename = bs_path + "/klaab.txt";
    ftpklaab.open(filename.c_str());
    for(it = 0; it < 150; it++)
      for(ik = 0; ik < 51; ik++)
        ftpklaab >> klaab[ik][it];
    ftpklaab.close();
}

