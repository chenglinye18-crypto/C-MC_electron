#include "Band.h"

void Band::VELOTET(void)
{
    //____Purpose : Calculate velocities within each tetrahedron
    //利用布里渊区四面体结构计算 速度
    //____local variables
    int it,ibase;
    double m[4][4],ar[4],vv,mi[4][4],det;  //matrix,mi=inverse matrix,

    //____loop over all tetrahedra   所有的四面体都会走一遍
    
    for(it=0;it<nt;it++)  // it=四面体[1]
    {

        //      ..build matrix   xkk ykk zkk k空间中的向量格点   （以每一个面的中心为格点的值？）
        m[1][1]=xkk[tet[1][it]]-xkk[tet[0][it]];
        m[1][2]=ykk[tet[1][it]]-ykk[tet[0][it]];
        m[1][3]=zkk[tet[1][it]]-zkk[tet[0][it]];
        m[2][1]=xkk[tet[2][it]]-xkk[tet[0][it]];
        m[2][2]=ykk[tet[2][it]]-ykk[tet[0][it]];
        m[2][3]=zkk[tet[2][it]]-zkk[tet[0][it]];
        m[3][1]=xkk[tet[3][it]]-xkk[tet[0][it]];
        m[3][2]=ykk[tet[3][it]]-ykk[tet[0][it]];
        m[3][3]=zkk[tet[3][it]]-zkk[tet[0][it]];
        //      ..determinat  行列式
        //   [ x1,y1,z1]
        //   [ x2,y2,z2]
        //   [ x3,y3,z3]   det 进行行列式求值
        det=m[1][1]*(m[2][2]*m[3][3]-m[3][2]*m[2][3])
            -m[1][2]*(m[2][1]*m[3][3]-m[3][1]*m[2][3])
            +m[1][3]*(m[2][1]*m[3][2]-m[3][1]*m[2][2]);

        //      ..inverse matrix
        mi[1][1]=(m[2][2]*m[3][3]-m[3][2]*m[2][3])/det;
        mi[2][1] =-(m[2][1]*m[3][3]-m[3][1]*m[2][3])/det;
        mi[3][1]=(m[2][1]*m[3][2]-m[3][1]*m[2][2])/det;
        mi[1][2] =-(m[1][2]*m[3][3]-m[3][2]*m[1][3])/det;
        mi[2][2]=(m[1][1]*m[3][3]-m[3][1]*m[1][3])/det;
        mi[3][2] =-(m[1][1]*m[3][2]-m[3][1]*m[1][2])/det;
        mi[1][3]=(m[1][2]*m[2][3]-m[2][2]*m[1][3])/det;
        mi[2][3] =-(m[1][1]*m[2][3]-m[2][1]*m[1][3])/det;
        mi[3][3]=(m[1][1]*m[2][2]-m[2][1]*m[1][2])/det;
        
        //k空间中的格点的能量值
        
        ar[1]=eek[tet[1][it]]-eek[tet[0][it]];  
        ar[2]=eek[tet[2][it]]-eek[tet[0][it]];
        ar[3]=eek[tet[3][it]]-eek[tet[0][it]];
        // ..group velocity within tetrahedron
        ibase=4*it;              // 4*(it-1)  it changed
        vt[ibase+0]=mi[1][1]*ar[1]+mi[1][2]*ar[2]+mi[1][3]*ar[3];  //这里的mi[x][x]已经是inverse之后的了   群速度的定义 蒙特卡洛模拟资料
        vt[ibase+1]=mi[2][1]*ar[1]+mi[2][2]*ar[2]+mi[2][3]*ar[3];
        vt[ibase+2]=mi[3][1]*ar[1]+mi[3][2]*ar[2]+mi[3][3]*ar[3];
        vt[ibase+3]=eek[tet[0][it]];


        //      ..absolute value of group velocity
        vv=sqrt(vt[ibase+0]*vt[ibase+0]+vt[ibase+1]*vt[ibase+1]+vt[ibase+2]*vt[ibase+2]);
        vgt[it]=vv;

        //      .. check velocity

    }
    //____End of VELOTET
    return;
}


void Band::ZUDI(void)
{
    //________ Purpose   : calculate the DOS-table and maximum DOS in tet.
    //____local variables
    int itab,it,ib,ie,itl,is;
    double eps,area,epsu,epsm,epsl;
    double  areal,areau;
    //____set all dos-values
    for(itab=0;itab<=MTAB;itab++)
    {
        //DO itab=0,MTAB
        for(ib=0;ib<MNB;ib++)
        {
            //DO ib=1,MNB
            dos[ib][itab]=0.0;
            dostetmax[ib][itab]=0.0;
        }
    }
    //____first: claculate DOS
    //____loop over all energy-tab-steps
    for(itab=0;itab<=MTAB;itab++)
    {
        //DO itab=0,MTAB
        eps=energy[itab];
        //_______energy index of tetraheda list
        ie=int(eps*dlist);       //+ 1//minim=0; dlist=网格的间距
        if(ie<0)ie=0;            //1;
        if(ie>=MWLE)ie=MWLE-1;   //maxim=
        //_______cal. DOS for each tetrahedron
        //      ..loop over all tetraheda,which may contain the energy eps
        for(ib=0;ib<nbt;ib++)   //nbt=总的能带数
        {
            //DO ib=1,nbt
            for(itl=ptlist[ie][ib];itl<=ptlist[ie][ib]+ntlist[ie][ib]-1;itl++)
            {
                // DO itl=
                it=tlist[itl];
                //         ..check if the tetrahedron contains the energy eps
                if(eps>eek[tet[0][it]]&&eps<eek[tet[3][it]]) //这里的条件为什么这样写？
                {
                    //             ..calculate area of intersection
                    area=FASTSURF(eps,it);
                    //             ..cal. density of states
                                 //**3
                    dos[ibt[it]][itab]+=area/vgt[it]*48.0/(8.0*PI*PI*PI);  //ibt指的是每个能带的编号 ，小尺寸蒙特卡洛模拟447页有解释
                }  //最终求出了对应的态密度
            }
        }
        //____end of energy-loop
    }
    //____write dos table to file
    
    //____cal. maximum dos in tetrahedron
    //____loop over all energy-tab-steps

    for(itab=1;itab<=MTAB;itab++)///itab=0 is set to 0;
    {
        //DO itab=1,MTAB
        //_______cal. maximum area for each tetrahedron between the energies
        //       epsl and epsu
        //      .. loop over all tetraheda
        for(it=0;it<nt;it++)
        {
            //DO it=1,nt
            epsl=energy[itab-1];
            epsu=energy[itab];
            //         ..check,if the tetrahedron overlapps with the given energy range
            if((epsu>eek[tet[0][it]])&&(epsl<eek[tet[3][it]]))
            {
                if(epsu<=eek[tet[1][it]])
                    //            ..the area grows quadratically between the energy of
                    //              node 1 and 2
                    area=FASTSURF(epsu,it);
                else if(epsl>=eek[tet[2][it]])
                    //            ..the area decreases quadratically between the energy of
                    //              node 3 and 4
                    area=FASTSURF(epsl,it);
                else
                {
                    //               ..between the energy of node 2 and 3 the area has a maximum
                    //                 the maximum in the interval between epsl and epsu is found
                    //                 with bisection (it has been assumed,that 64 iterations are
                    //                 sufficient).
                    //               ..limit the energy interval to the range between the energy
                    //                 of node 2 and 3
                    epsl=Max(epsl,eek[tet[1][it]]);
                    epsu=Min(epsu,eek[tet[2][it]]);
                    areal=FASTSURF(epsl,it);
                    areau=FASTSURF(epsu,it);
                    for(is=0;is<64;is++)  //迭代比较64次才能比较出来最终结果
                    {
                        //DO is=1,64
                        epsm=0.50*(epsl+epsu);
                        if(areau>areal)
                        {
                            epsl=epsm;
                            areal=FASTSURF(epsl,it);
                        }
                        else
                        {
                            epsu=epsm;
                            areau=FASTSURF(epsu,it);
                        }
                    }
                    area=Max(areal,areau);
                }
                //            ..take maximum
                dostetmax[ibt[it]][itab]=Max(area/vgt[it],
                    dostetmax[ibt[it]][itab]);
            }
        }
        //____end of energy-loop
    }
    //____set maximum area for zero energy
    for(ib=0;ib<MNB;ib++)
        // DO ib=1,MNB
        dostetmax[ib][0]=dostetmax[ib][1];

    //____End of ZUDI
    return;
}


//====
double Band::SURF(double eps,int it)
{
    //     double eps
    //    int it
    //    : Calculate the area of the intersection of an equi energy plane with
    //      a tetrahedron
    //
    //    : between the energies of the node 1 and 2 the intersection is always
    //      of triangular shape
    //      between node 2 and 3 always a quadrangle and between 3 and 4 again
    //      a triangle
    //____local variables
    double rfac,ax,ay,az,bx,by,bz,cx,cy,cz;
    double dx,dy,dz,a1,a2,a3,a4; //,FII;
    double rtnv;
    if(eps>eek[tet[0][it]]&&eps<eek[tet[3][it]])
    {
        //_______triangular intersection with energy lower e2
        if (eps<=eek[tet[1][it]])
        {
            //         ..cal. vectors (relative to origin of tetrahedron) pointing to the
            //           three nodes of the triangle
            rfac= (eps-eek[tet[0][it]])/(eek[tet[1][it]]-eek[tet[0][it]]);
            ax=rfac*(xkk[tet[1][it]]-xkk[tet[0][it]]);
            ay=rfac*(ykk[tet[1][it]]-ykk[tet[0][it]]);
            az=rfac*(zkk[tet[1][it]]-zkk[tet[0][it]]);
            if (rfac<0.0||rfac>1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);;
            }
            rfac=(eps-eek[tet[0][it]])/(eek[tet[2][it]]-eek[tet[0][it]]);
            bx=rfac*(xkk[tet[2][it]]-xkk[tet[0][it]]);
            by=rfac*(ykk[tet[2][it]]-ykk[tet[0][it]]);
            bz=rfac*(zkk[tet[2][it]]-zkk[tet[0][it]]);
            if (rfac<0.0||rfac>1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            rfac=(eps-eek[tet[0][it]])/(eek[tet[3][it]]-eek[tet[0][it]]);
            cx=rfac*(xkk[tet[3][it]]-xkk[tet[0][it]]);
            cy=rfac*(ykk[tet[3][it]]-ykk[tet[0][it]]);
            cz=rfac*(zkk[tet[3][it]]-zkk[tet[0][it]]);
            if (rfac<0.0||rfac>1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            //         ..cal. area
            rtnv=FII(ax,ay,az,bx,by,bz,cx,cy,cz);
        }

        //_______quadrangular intersection with energy between e2 and e3
        else if(eps>eek[tet[1][it]]&&eps<eek[tet[2][it]])
        {
            //         ..cal. vectors (relative to origin of k-space) pointing to the
            //           four nodes of the quadrangle
            rfac=(eps-eek[tet[0][it]])/(eek[tet[2][it]]-eek[tet[0][it]]);
            ax=rfac*(xkk[tet[2][it]]-xkk[tet[0][it]])+xkk[tet[0][it]];
            ay=rfac*(ykk[tet[2][it]]-ykk[tet[0][it]])+ykk[tet[0][it]];
            az=rfac*(zkk[tet[2][it]]-zkk[tet[0][it]])+zkk[tet[0][it]];
            if (rfac< 0.0||rfac> 1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            rfac=(eps-eek[tet[1][it]])
                /(eek[tet[2][it]]-eek[tet[1][it]]);
            bx=rfac*(xkk[tet[2][it]]-xkk[tet[1][it]])+xkk[tet[1][it]];
            by=rfac*(ykk[tet[2][it]]-ykk[tet[1][it]])+ykk[tet[1][it]];
            bz=rfac*(zkk[tet[2][it]]-zkk[tet[1][it]])+zkk[tet[1][it]];
            if (rfac< 0.0||rfac> 1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            rfac= (eps-eek[tet[0][it]])
                /(eek[tet[3][it]]-eek[tet[0][it]]);
            cx=rfac*(xkk[tet[3][it]]-xkk[tet[0][it]])+xkk[tet[0][it]];
            cy=rfac*(ykk[tet[3][it]]-ykk[tet[0][it]])+ykk[tet[0][it]];
            cz=rfac*(zkk[tet[3][it]]-zkk[tet[0][it]])+zkk[tet[0][it]];
            if (rfac< 0.0||rfac> 1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            rfac= (eps-eek[tet[1][it]])
                /(eek[tet[3][it]]-eek[tet[1][it]]);
            dx=rfac*(xkk[tet[3][it]]-xkk[tet[1][it]])+xkk[tet[1][it]];
            dy=rfac*(ykk[tet[3][it]]-ykk[tet[1][it]])+ykk[tet[1][it]];
            dz=rfac*(zkk[tet[3][it]]-zkk[tet[1][it]])+zkk[tet[1][it]];
            if (rfac< 0.0||rfac> 1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            //      ..cal. area of the four possible triangles made from the four nodes
            //        of the quadrangle
            a1=FII(ax,ay,az,bx,by,bz,cx,cy,cz);
            a2=FII(dx,dy,dz,bx,by,bz,cx,cy,cz);
            a3=FII(ax,ay,az,dx,dy,dz,cx,cy,cz);
            a4=FII(ax,ay,az,bx,by,bz,dx,dy,dz);
            rtnv=0.50*(a1+a2+a3+a4);

        }
        //_______triangular intersection with energy higher e3
        else if(eps>eek[tet[2][it]]&&eps<eek[tet[3][it]])
        {
            //         ..cal. vectors (relative to origin of tetrahedron) pointing to the
            //           three nodes of the triangle
            rfac= (eps-eek[tet[3][it]])/(eek[tet[0][it]]-eek[tet[3][it]]);
            ax=rfac*(xkk[tet[0][it]]-xkk[tet[3][it]]);
            ay=rfac*(ykk[tet[0][it]]-ykk[tet[3][it]]);
            az=rfac*(zkk[tet[0][it]]-zkk[tet[3][it]]);
            if (rfac< 0.0||rfac> 1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            rfac= (eps-eek[tet[3][it]])/(eek[tet[1][it]]-eek[tet[3][it]]);
            bx=rfac*(xkk[tet[1][it]]-xkk[tet[3][it]]);
            by=rfac*(ykk[tet[1][it]]-ykk[tet[3][it]]);
            bz=rfac*(zkk[tet[1][it]]-zkk[tet[3][it]]);
            if (rfac< 0.0||rfac> 1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            rfac= (eps-eek[tet[3][it]])/(eek[tet[2][it]]-eek[tet[3][it]]);
            cx=rfac*(xkk[tet[2][it]]-xkk[tet[3][it]]);
            cy=rfac*(ykk[tet[2][it]]-ykk[tet[3][it]]);
            cz=rfac*(zkk[tet[2][it]]-zkk[tet[3][it]]);
            if (rfac< 0.0||rfac> 1.0)
            {
                cout<<"error SURF: rfac< 0 or>1";exit(0);
            }
            rtnv=FII(ax,ay,az,bx,by,bz,cx,cy,cz);
        }
        else
            rtnv=0.0;
    }
    else
        rtnv=0.0;
    //____End of SURF
    return rtnv;
}


////====
double Band:: FII(double x1,double y1,double z1,
				  double x2,double y2,double z2,
				  double x3,double y3,double z3)
{

    double sbx,sby,sbz,rtnv;
    //____Purpose : calculate the square in a triangle
    //    Parameter : values of the three vectors pointing to the nodes
    //    : The area is given by the absolute vale of the cross product
    //      of the two vectors along the edges between the nodes 2/1 and 3/1

    sbx=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
    sby=(x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
    sbz=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
    rtnv=0.50*sqrt(sbx*sbx+sby*sby+sbz*sbz);
    //____End of FII
    return rtnv;
}


//====
void Band::HESSTET(void)
{
    //____Purpose : cal. Hessian form for the four surfaces of the tetraheda,
    //              the normal vectors are not normalized,
    //              they allways point into the tetrahedron.
    //____local variables
    double ax,ay,az,bx,by,bz,cx,cy,cz;
    double dx,dy,dz,ex,ey,ez,fx,fy,fz;
    int it;
	//____loop over all tetraheda
    for(it=0;it<nt;it++)
    {
        //DO it=1,nt
        //      ..cal. the six vectors along the edges of the tetrahedron
        ax=xkk[tet[1][it]]-xkk[tet[0][it]];
        ay=ykk[tet[1][it]]-ykk[tet[0][it]];
        az=zkk[tet[1][it]]-zkk[tet[0][it]];
        bx=xkk[tet[2][it]]-xkk[tet[1][it]];
        by=ykk[tet[2][it]]-ykk[tet[1][it]];
        bz=zkk[tet[2][it]]-zkk[tet[1][it]];
        cx=xkk[tet[0][it]]-xkk[tet[2][it]];
        cy=ykk[tet[0][it]]-ykk[tet[2][it]];
        cz=zkk[tet[0][it]]-zkk[tet[2][it]];
        dx=xkk[tet[3][it]]-xkk[tet[0][it]];
        dy=ykk[tet[3][it]]-ykk[tet[0][it]];
        dz=zkk[tet[3][it]]-zkk[tet[0][it]];
        ex=xkk[tet[3][it]]-xkk[tet[1][it]];
        ey=ykk[tet[3][it]]-ykk[tet[1][it]];
        ez=zkk[tet[3][it]]-zkk[tet[1][it]];
        fx=xkk[tet[3][it]]-xkk[tet[2][it]];
        fy=ykk[tet[3][it]]-ykk[tet[2][it]];
        fz=zkk[tet[3][it]]-zkk[tet[2][it]];
        //      .. surface 1
        //      .. cal. vector perpendicular to the surface (cross product)
        datantlin[0][0][it]=ey*fz-ez*fy;
        datantlin[1][0][it]=ez*fx-ex*fz;
        datantlin[2][0][it]=ex*fy-ey*fx;
        //      .. cal. distance of surface to origin of tetrahedron
        datantlin[3][0][it]=datantlin[0][0][it]*xkk[tet[1][it]]
            +datantlin[1][0][it]*ykk[tet[1][it]]
            +datantlin[2][0][it]*zkk[tet[1][it]];
        //      .. force direction of perpendicular vector to point into the tetrahedron
        if ( datantlin[0][0][it]*xkk[tet[0][it]]
            +datantlin[1][0][it]*ykk[tet[0][it]]
            +datantlin[2][0][it]*zkk[tet[0][it]]<datantlin[3][0][it])
        {
            datantlin[0][0][it] =-datantlin[0][0][it];
            datantlin[1][0][it] =-datantlin[1][0][it];
            datantlin[2][0][it] =-datantlin[2][0][it];
            datantlin[3][0][it] =-datantlin[3][0][it];
        }

        //      .. surface 2
        datantlin[0][1][it]=dy*fz-dz*fy;
        datantlin[1][1][it]=dz*fx-dx*fz;
        datantlin[2][1][it]=dx*fy-dy*fx;
        datantlin[3][1][it]=datantlin[0][1][it]*xkk[tet[2][it]]
            +datantlin[1][1][it]*ykk[tet[2][it]]
            +datantlin[2][1][it]*zkk[tet[2][it]];
        if ( datantlin[0][1][it]*xkk[tet[1][it]]
            +datantlin[1][1][it]*ykk[tet[1][it]]
            +datantlin[2][1][it]*zkk[tet[1][it]]<datantlin[3][1][it])
        {
            datantlin[0][1][it] =-datantlin[0][1][it];
            datantlin[1][1][it] =-datantlin[1][1][it];
            datantlin[2][1][it] =-datantlin[2][1][it];
            datantlin[3][1][it] =-datantlin[3][1][it];
        }
        //      .. surface 3
        datantlin[0][2][it]=ey*dz-ez*dy;
        datantlin[1][2][it]=ez*dx-ex*dz;
        datantlin[2][2][it]=ex*dy-ey*dx;
        datantlin[3][2][it]=datantlin[0][2][it]*xkk[tet[3][it]]
            +datantlin[1][2][it]*ykk[tet[3][it]]
            +datantlin[2][2][it]*zkk[tet[3][it]];
        if ( datantlin[0][2][it]*xkk[tet[2][it]]
            +datantlin[1][2][it]*ykk[tet[2][it]]
            +datantlin[2][2][it]*zkk[tet[2][it]]<datantlin[3][2][it])
        {
            datantlin[0][2][it] =-datantlin[0][2][it];
            datantlin[1][2][it] =-datantlin[1][2][it];
            datantlin[2][2][it] =-datantlin[2][2][it];
            datantlin[3][2][it] =-datantlin[3][2][it];
        }
        //      .. surface 4
        datantlin[0][3][it]=ay*bz-az*by;
        datantlin[1][3][it]=az*bx-ax*bz;
        datantlin[2][3][it]=ax*by-ay*bx;
        datantlin[3][3][it]=datantlin[0][3][it]*xkk[tet[0][it]]
            +datantlin[1][3][it]*ykk[tet[0][it]]
            +datantlin[2][3][it]*zkk[tet[0][it]] ;
        if ( datantlin[0][3][it]*xkk[tet[3][it]]
            +datantlin[1][3][it]*ykk[tet[3][it]]
            +datantlin[2][3][it]*zkk[tet[3][it]]<datantlin[3][3][it])
        {
            datantlin[0][3][it] =-datantlin[0][3][it];
            datantlin[1][3][it] =-datantlin[1][3][it];
            datantlin[2][3][it] =-datantlin[2][3][it];
            datantlin[3][3][it] =-datantlin[3][3][it];
        }
    }
	

    //____End of HESSTET
    return ;
}


//=====

void Band::SETNEIBSYM(void)
{
    //Purpose : cal. the five neighbouring wedges for each of the 48 wedges

    int jsym,iv,isym;
    int pos[3];
    double ax,ay,az,xkl,ykl,zkl;
    double xtrans[5],ytrans[5],ztrans[5];
    double xk,yk,zk;
    double TempMin,TempMax;
    int MinK,MaxK;
    //five vectors which are outside of the irreducible wedge
    //    and within the five neighbour wedges

    //Kz=0
    xtrans[0]= 0.5*a0pi;
    ytrans[0]= 0.3*a0pi;
    ztrans[0]=-0.01*a0pi;
    //Kx-Ky=0
    xtrans[1]= 0.49*a0pi;
    ytrans[1]= 0.51*a0pi;
    ztrans[1]= 0.1*a0pi;
    //Ky-Kz=0
    xtrans[2]= 0.5*a0pi;
    ytrans[2]= 0.24*a0pi;
    ztrans[2]= 0.26*a0pi;
    //Kx+Ky+Kz=1.5
    xtrans[3]= 0.71*a0pi;
    ytrans[3]= 0.51*a0pi;
    ztrans[3]= 0.31*a0pi;
    //Kx=1
    xtrans[4]= 1.01*a0pi;
    ytrans[4]= 0.1*a0pi;
    ztrans[4]= 0.05*a0pi;
    // 

    //loop over all 48 wedges of the first BZ
    for(isym=0;isym<48;isym++)
    {
        //loop over all five neighbours
        for(iv=0;iv<=4;iv++)     /////////becareful
        {
            ax=xtrans[iv];
            ay=ytrans[iv];
            az=ztrans[iv];
            OutWedge(xk,yk,zk,ax,ay,az,isym);
            xkl=xk;
            ykl=yk;
            zkl=zk;

            //in the case that the vector is outside of the first BZ
            //reduce it to the first BZ
            CONFBZ(xk,yk,zk);

            //save k-vector displacement
            xbz[iv][isym]=xk-xkl;
            ybz[iv][isym]=yk-ykl;
            zbz[iv][isym]=zk-zkl;
            if(fabs(xbz[iv][isym])<1e-10) xbz[iv][isym]=0.0;
            if(fabs(ybz[iv][isym])<1e-10) ybz[iv][isym]=0.0;
            if(fabs(zbz[iv][isym])<1e-10) zbz[iv][isym]=0.0;
            //cal. index of the new wedge
            pos[0]=pos[1]=pos[2]=1;
			if(fabs(xk)>fabs(yk))
            {
                TempMin=fabs(yk);
                TempMax=fabs(xk);
                MinK=1;
                MaxK=0;
            }
            else
            {
                TempMin=fabs(xk);
                TempMax=fabs(yk);
                MinK=0;
                MaxK=1;
            }
            if(fabs(zk)<TempMin) MinK=2;
            else if(fabs(zk)>TempMax) MaxK=2;
            pos[MinK]=2;
            pos[MaxK]=0;         //pos[a]=bakb   a=0~2 b=0~2
            jsym=indmat[sign(xk)+1][sign(yk)+1][sign(zk)+1][pos[0]][pos[1]][pos[2]];
            newisym[iv][isym]=jsym;
        }
	}

}


void Band::OutWedge(double &xout,double &yout,double &zout,
double xin,double yin,double zin,int isym)
{
    //transform k-vector from the irreducible wedge into the BZ

    double  vsort[3];

    vsort[0]=xin;
    vsort[1]=yin;
    vsort[2]=zin;

    xout=vsort[matsym[0][isym]-1]*matsym[3][isym];
    yout=vsort[matsym[1][isym]-1]*matsym[4][isym];
    zout=vsort[matsym[2][isym]-1]*matsym[5][isym];

}


void Band::CONFBZ (double &xkl,double &ykl,double &zkl)
{
    //     Purpose: confine k-vector to first BZ
    //     -------
    //____local variables
    bool flag,ndfl;

    //____find the equivalent k-vektor in the
    //     first brillouin zone
    ndfl=false;
    while(!ndfl)
    {
        flag=false;
        if(fabs(xkl)>a0pi)
            flag=true;
        else if (fabs(ykl)>a0pi)
            flag=true;
        else if (fabs(zkl)>a0pi)
            flag=true;
        else if (fabs(xkl)+fabs(ykl)+fabs(zkl)>1.50*a0pi)
            flag=true;

        if(flag)
        {
            if (xkl>0.0)
                xkl=xkl-a0pi;
            else
                xkl=xkl+a0pi;
            if (ykl>0.0)
                ykl=ykl-a0pi;
            else
                ykl=ykl+a0pi;
            if (zkl>0.0)
                zkl=zkl-a0pi;
            else
                zkl=zkl+a0pi;
        }
        else
            ndfl=true;
    }                            //GOTO 10

    //____end of CONFBZ
    return;
}


void Band::CALCFAC(void)
{
    //____Purpose : Calculate factors for area calculation within each tetrahedron
    //              used in FASTSURF
    //
    //    : In the energy intervals between the energy of the nodes the
    //      dependence of the area on energy is given by a quadrati//relation
    //____local variables
    int it;
    double eps;
    //____loop over all tetrahedra
    for(it=0;it<nt;it++)
    {
        //DO it=1,nt
        //      .. between e1 and e2 the intersection is triangular
        if (eek[tet[1][it]]>eek[tet[0][it]])
        {
            eps=0.50*(eek[tet[0][it]]+eek[tet[1][it]]);
            faclow[it]=SURF(eps,it)/((eps-eek[tet[0][it]])*(eps-eek[tet[0][it]]));
        }
        else
            faclow[it]=0.0;
        //      .. between e3 and e4 the intersection is triangular
        if (eek[tet[3][it]]>eek[tet[2][it]])
        {
            eps=0.50*(eek[tet[2][it]]+eek[tet[3][it]]);
            fachigh[it]=SURF(eps,it)/((eps-eek[tet[3][it]])*(eps-eek[tet[3][it]]));
        }
        else
            fachigh[it]=0.0;
        //      .. between e2 and e3 the intersection is a quadrangle
        if (eek[tet[2][it]]>eek[tet[1][it]])
        {
            eps=0.5*(eek[tet[2][it]]+eek[tet[1][it]]);
            facmedium[it]=(faclow[it]*(eps-eek[tet[0][it]])*(eps-eek[tet[0][it]])
                -SURF(eps,it))/((eps-eek[tet[1][it]])*(eps-eek[tet[1][it]]));
        }
        else
            facmedium[it]=0.0;
    }
    //____End of CALCFAC
    return ;
}


//====
double Band::FASTSURF(double eps,int it)
{
    double rtnv;
    //int it;
    //    : Calculate the area of the intersection of a equi energy plane with
    //      a tetrahedron
    if (eps>eek[tet[0][it]]&&eps<eek[tet[3][it]])
    {
        //_______triangular intersection with energy lower than e2
        if (eps<=eek[tet[1][it]])
            rtnv=faclow[it]*(eps-eek[tet[0][it]])*(eps-eek[tet[0][it]]);
        //_______quadrangular intersection with energy between e2 and e3
        else if (eps<eek[tet[2][it]])
        {

            if (eek[tet[1][it]]-eek[tet[0][it]]>1e-10)
                rtnv=faclow[it]*(eps-eek[tet[0][it]])*(eps-eek[tet[0][it]])
                    -facmedium[it]*(eps-eek[tet[1][it]])*(eps-eek[tet[1][it]]);
            else
                rtnv=SURF(eps,it);
        }
        //_______triangular intersection with energy higher than e3
        else
            rtnv=fachigh[it]*(eps-eek[tet[3][it]])*(eps-eek[tet[3][it]]);
    }
    else
        rtnv=0.0;
    if (rtnv<0.0)
    {
        /*
                 WRITE (LUTTYO,*) 'FASTSURF: Area<0'
                 WRITE (LUOUT ,*) ''
                 STOP
                */
        cout<<"error FASTSURF: Area<0";exit(0);
    }
    //____End of FASTSURF
    return rtnv;
}


//====
void Band::GETMAXAOVTET(void)
{
    //________ Purpose   : calculate the maximum AOV in tet.
    //                     this quantity is area over velocity
    //____local variables
    int it,is;
    double  area,epsu,epsm,epsl;
    double  areal,areau;
    //____loop over all tetrahedra
    for(it=0;it<nt;it++)
    {
        //DO it=1,nt
        //_______cal. maximum area for each tetrahedron between the energies e2 and e3
        //      ..between the energy of node 2 and 3 the area has a maximum
        //        the maximum in the interval between epsl and epsu is found
        //        with bisection (it has been assumed,that 20 iterations are
        //        sufficient).
        epsl=eek[tet[1][it]];
        epsu=eek[tet[2][it]];
        areal=FASTSURF(epsl,it);
        areau=FASTSURF(epsu,it);
        for(is=0;is<20;is++)
        {
            //DO is=1,20
            epsm=0.50*(epsl+epsu);
            if (areau>areal)
            {
                epsl=epsm;
                areal=FASTSURF(epsl,it);
            }
            else
            {
                epsu=epsm;
                areau=FASTSURF(epsu,it);
            }
        }
        area=Max(areal,areau);
        //   ..take maximum
        maxaovtet[it]=area/vgt[it];
    }
    //____End of GETMAXAOVTET
    return;
}


void  Band::READBS(void)
{
    //READ band structure data

    int  it,ik,iptype;

    //FILE*lutmp;
    ifstream ftp;
    string filename;
        //read tetrahedra data (not normalized), the nodes of the tetrahedra
        //must be sorted in such a way, that their energy is in ascending order
        //with the node number (e1 < e2 < e3 < e4)!

    if(sifl)
        {
          filename = pathname + "/bs.si.asc";
            ftp.open(filename.c_str());
            //assert(ftp);
        }
		else if (gefl)
		{
          filename = pathname + "/bs.ge.trs";
          ftp.open(filename.c_str());
		}
        else if(gaasfl)
        {
          filename = pathname + "/bs.gaas.asc";
            ftp.open(filename.c_str());
        }

        //Number of energy bands
        //READ(lutmp) nband
   
        for(int i=0;i<NumParType;i++)
          ftp>>nband[i];
             
        if(nband[PELEC]>NBE)
        {

            cout<<"error ETABF: nbe > MNBE";exit(0);
            //ENDIF
        }
        if(nband[PHOLE]>NBH)
        {
            cout<<"error ETABF: nbh > MNBH";exit(0);
            //ENDIF
        }
        if(nband[POXEL]>NBOE)
        {

            cout<<"error ETABF: nboe > MNBOE";exit(0);
            //ENDIF
        }
        nbt = nband[PELEC] + nband[PHOLE] + nband[POXEL];

        //set band pointers
        bandof[PELEC]=0;
        bandof[PHOLE]=nband[PELEC];
        bandof[POXEL]=nband[PELEC]+nband[PHOLE];
        for(iptype=0;iptype<NumParType;iptype++)
        {
            //DO iptype=1, NumParType
            for(int iband=bandof[iptype];iband<bandof[iptype]+nband[iptype];iband++)
            {
                //DO iband=bandof[iptype]+1, bandof[iptype]+nband[iptype]
                partyp[iband]=iptype;
                //ENDDO
            }
            //ENDDO
        }

        //Number of k-space points
        //READ(lutmp) nk
        ftp>>nk;

        if(nk>MNK)
        {

        }
        //read k-vector and energy

        for(ik=0;ik<nk;ik++)ftp>>xkk[ik];		
		for(ik=0;ik<nk;ik++)ftp>>ykk[ik];
		for(ik=0;ik<nk;ik++)ftp>>zkk[ik];
  
		for(ik=0;ik<nk;ik++)
          ftp>>eek[ik];
  
		for(ik=0;ik<nk;ik++)
        {
            //DO ik=1, nk
            if(eek[ik]<0.0)
            {

                cout<<"error ETABF: eek < 0";exit(0);
                //ENDIF
            }
        }
        //ENDDO

        //read number of tetrahedra

        ftp>>nt;
        if(nt>MNTet)
        {

            cout<<"error ETABF: nt > MNTet";exit(0);
            //ENDIF
        }
        //read grid point indices for each of the four corners of a tetrahedron
   
        int iitemp;
        
        for(it=0;it<nt;it++){ftp>>iitemp;tet[0][it]=iitemp-1;}
        for(it=0;it<nt;it++){ftp>>iitemp;tet[1][it]=iitemp-1;}
        for(it=0;it<nt;it++){ftp>>iitemp;tet[2][it]=iitemp-1;}
        for(it=0;it<nt;it++){ftp>>iitemp;tet[3][it]=iitemp-1;}

        for(it=0;it<nt;it++){ftp>>iitemp;NeibrTet[0][it]=iitemp-1;}
        for(it=0;it<nt;it++){ftp>>iitemp;NeibrTet[1][it]=iitemp-1;}
        for(it=0;it<nt;it++){ftp>>iitemp;NeibrTet[2][it]=iitemp-1;}
        for(it=0;it<nt;it++){ftp>>iitemp;NeibrTet[3][it]=iitemp-1;}
        //read band index of tetrahedron

        for(it=0;it<nt;it++){ftp>>iitemp;ibt[it]=iitemp-1;}
        ftp.close();             //CLOSE(lutmp)

		//normalize k-space grid point data

        for(ik=0;ik<nk;ik++)
        {
            //DO ik=1, nk
            xkk[ik]=xkk[ik]*a0pi;
            ykk[ik]=ykk[ik]*a0pi;
            zkk[ik]=zkk[ik]*a0pi;
            eek[ik]=eek[ik]/eV0;
            //ENDDO
        }
        //check order of the energy of at the corners of each tetrahedron
        for(it=0;it<nt;it++)
        {
            //DO it=1, nt
            if(!((eek[tet[0][it]]<=eek[tet[1][it]])||
                (eek[tet[1][it]]<=eek[tet[2][it]])||
                (eek[tet[2][it]]<=eek[tet[3][it]])))
            {

                cout<<"error ETABF: energy order wrong";exit(0);
                //ENDIF
            }
            //ENDDO
        }

    //cal. center of tetrahedron
    for(it=0;it<nt;it++)
    {
        //DO it=1, nt
        xkct[it]=0.250*(xkk[tet[0][it]]+xkk[tet[1][it]]
            +xkk[tet[2][it]]+xkk[tet[3][it]]);
        ykct[it]=0.250*(ykk[tet[0][it]]+ykk[tet[1][it]]+ykk[tet[2][it]]
            +ykk[tet[3][it]]);
        zkct[it]=0.250*(zkk[tet[0][it]]+zkk[tet[1][it]]+zkk[tet[2][it]]
            +zkk[tet[3][it]]);
        //ENDDO
    }
    //cal. velocities within each tetrahedron
    VELOTET();
    //get vectors normal to the four surfaces of each tetrahedron
    //(Hessian form)
    HESSTET();
    
    //matrix operations (the 48 symmetry transformations of the BZ)
    SETMATSYM();
    
    //cal. the five neighbouring wedges for each of the 48 wedges
    SETNEIBSYM();
    
    //cal. factors for the determination of area within a tetrahedron
    //used in FASTSURF
    CALCFAC();
    
    //calculate maximum area over velocity for each tetrahedron
    GETMAXAOVTET();
    
    //check tetrahedra

    //End of READBS
    return;
    //END
}


//subroutine
//corrected by dugang DEC.15
void Band::SETMATSYM(void)
{
    //Purpose:   initialize the 48 symmetry transformations of silicon
    //--------
    //local variables
    int isym;

    //matrix operations (the 48 symmetry transformations of the BZ)

	  matsym[3][0]=1;
      matsym[4][0]=1;
      matsym[5][0]=1;
      matsym[0][0]=1;
      matsym[1][0]=2;
      matsym[2][0]=3;
      matsym[3][1]=1;
      matsym[4][1]=1;
      matsym[5][1]=-1;
      matsym[0][1]=1;
      matsym[1][1]=2;
      matsym[2][1]=3;
      matsym[3][2]=1;
      matsym[4][2]=-1;
      matsym[5][2]=1;
      matsym[0][2]=1;
      matsym[1][2]=2;
      matsym[2][2]=3;
      matsym[3][3]=1;
      matsym[4][3]=-1;
      matsym[5][3]=-1;
      matsym[0][3]=1;
      matsym[1][3]=2;
      matsym[2][3]=3;
      matsym[3][4]=1;
      matsym[4][4]=1;
      matsym[5][4]=1;
      matsym[0][4]=1;
      matsym[1][4]=3;
      matsym[2][4]=2;
      matsym[3][5]=1;
      matsym[4][5]=1;
      matsym[5][5]=-1;
      matsym[0][5]=1;
      matsym[1][5]=3;
      matsym[2][5]=2;
      matsym[3][6]=1;
      matsym[4][6]=-1;
      matsym[5][6]=1;
      matsym[0][6]=1;
      matsym[1][6]=3;
      matsym[2][6]=2;
      matsym[3][7]=1;
      matsym[4][7]=-1;
      matsym[5][7]=-1;
      matsym[0][7]=1;
      matsym[1][7]=3;
      matsym[2][7]=2;
      matsym[3][8]=-1;
      matsym[4][8]=1;
      matsym[5][8]=1;
      matsym[0][8]=1;
      matsym[1][8]=2;
      matsym[2][8]=3;
      matsym[3][9]=-1;
      matsym[4][9]=1;
      matsym[5][9]=-1;
      matsym[0][9]=1;
      matsym[1][9]=2;
      matsym[2][9]=3;
      matsym[3][10]=-1;
      matsym[4][10]=-1;
      matsym[5][10]=1;
      matsym[0][10]=1;
      matsym[1][10]=2;
      matsym[2][10]=3;
      matsym[3][11]=-1;
      matsym[4][11]=-1;
      matsym[5][11]=-1;
      matsym[0][11]=1;
      matsym[1][11]=2;
      matsym[2][11]=3;
      matsym[3][12]=-1;
      matsym[4][12]=1;
      matsym[5][12]=1;
      matsym[0][12]=1;
      matsym[1][12]=3;
      matsym[2][12]=2;
      matsym[3][13]=-1;
      matsym[4][13]=1;
      matsym[5][13]=-1;
      matsym[0][13]=1;
      matsym[1][13]=3;
      matsym[2][13]=2;
      matsym[3][14]=-1;
      matsym[4][14]=-1;
      matsym[5][14]=1;
      matsym[0][14]=1;
      matsym[1][14]=3;
      matsym[2][14]=2;
      matsym[3][15]=-1;
      matsym[4][15]=-1;
      matsym[5][15]=-1;
      matsym[0][15]= 1;
      matsym[1][15]= 3;
      matsym[2][15]= 2;
      matsym[3][16]=1;
      matsym[4][16]=1;
      matsym[5][16]=1;
      matsym[0][16]=2;
      matsym[1][16]=1;
      matsym[2][16]=3;
      matsym[3][17]=1;
      matsym[4][17]=1;
      matsym[5][17]=-1;
      matsym[0][17]=2;
      matsym[1][17]=1;
      matsym[2][17]=3;
      matsym[3][18]=-1;
      matsym[4][18]= 1;
      matsym[5][18]=-1;
      matsym[0][18]= 2;
      matsym[1][18]= 1;
      matsym[2][18]= 3;
      matsym[3][19]=-1;
      matsym[4][19]= 1;
      matsym[5][19]= 1;
      matsym[0][19]= 2;
      matsym[1][19]= 1;
      matsym[2][19]= 3;
      matsym[3][20]=1;
      matsym[4][20]=1;
      matsym[5][20]=1;
      matsym[0][20]=3;
      matsym[1][20]=1;
      matsym[2][20]=2;
      matsym[3][21]=1;
      matsym[4][21]=1;
      matsym[5][21]=-1;
      matsym[0][21]=3;
      matsym[1][21]=1;
      matsym[2][21]=2;
      matsym[3][22]=-1;
      matsym[4][22]= 1;
      matsym[5][22]=-1;
      matsym[0][22]= 3;
      matsym[1][22]= 1;
      matsym[2][22]= 2;
      matsym[3][23]=-1;
      matsym[4][23]= 1;
      matsym[5][23]= 1;
      matsym[0][23]= 3;
      matsym[1][23]= 1;
      matsym[2][23]= 2;
      matsym[3][24]=1;
      matsym[4][24]=-1;
      matsym[5][24]=1;
      matsym[0][24]=2;
      matsym[1][24]=1;
      matsym[2][24]=3;
      matsym[3][25]=1;
      matsym[4][25]=-1;
      matsym[5][25]=-1;
      matsym[0][25]=2;
      matsym[1][25]=1;
      matsym[2][25]=3;
      matsym[3][26]=-1;
      matsym[4][26]=-1;
      matsym[5][26]=-1;
      matsym[0][26]= 2;
      matsym[1][26]= 1;
      matsym[2][26]= 3;
      matsym[3][27]=-1;
      matsym[4][27]=-1;
      matsym[5][27]= 1;
      matsym[0][27]= 2;
      matsym[1][27]= 1;
      matsym[2][27]= 3;
      matsym[3][28]=1;
      matsym[4][28]=-1;
      matsym[5][28]=1;
      matsym[0][28]=3;
      matsym[1][28]=1;
      matsym[2][28]=2;
      matsym[3][29]=1;
      matsym[4][29]=-1;
      matsym[5][29]=-1;
      matsym[0][29]=3;
      matsym[1][29]=1;
      matsym[2][29]=2;
      matsym[3][30]=-1;
      matsym[4][30]=-1;
      matsym[5][30]=-1;
      matsym[0][30]= 3;
      matsym[1][30]= 1;
      matsym[2][30]= 2;
      matsym[3][31]=-1;
      matsym[4][31]=-1;
      matsym[5][31]= 1;
      matsym[0][31]= 3;
      matsym[1][31]= 1;
      matsym[2][31]= 2;
      matsym[3][32]=1;
      matsym[4][32]=1;
      matsym[5][32]=1;
      matsym[0][32]=2;
      matsym[1][32]=3;
      matsym[2][32]=1;
      matsym[3][33]=1;
      matsym[4][33]=-1;
      matsym[5][33]=1;
      matsym[0][33]=2;
      matsym[1][33]=3;
      matsym[2][33]=1;
      matsym[3][34]=-1;
      matsym[4][34]=-1;
      matsym[5][34]= 1;
      matsym[0][34]= 2;
      matsym[1][34]= 3;
      matsym[2][34]= 1;
      matsym[3][35]=-1;
      matsym[4][35]= 1;
      matsym[5][35]= 1;
      matsym[0][35]= 2;
      matsym[1][35]= 3;
      matsym[2][35]= 1;
      matsym[3][36]=1;
      matsym[4][36]=1;
      matsym[5][36]=1;
      matsym[0][36]=3;
      matsym[1][36]=2;
      matsym[2][36]=1;
      matsym[3][37]=1;
      matsym[4][37]=-1;
      matsym[5][37]=1;
      matsym[0][37]=3;
      matsym[1][37]=2;
      matsym[2][37]=1;
      matsym[3][38]=-1;
      matsym[4][38]=-1;
      matsym[5][38]= 1;
      matsym[0][38]= 3;
      matsym[1][38]= 2;
      matsym[2][38]= 1;
      matsym[3][39]=-1;
      matsym[4][39]= 1;
      matsym[5][39]= 1;
      matsym[0][39]= 3;
      matsym[1][39]= 2;
      matsym[2][39]= 1;
      matsym[3][40]=1;
      matsym[4][40]=1;
      matsym[5][40]=-1;
      matsym[0][40]=2;
      matsym[1][40]=3;
      matsym[2][40]=1;
      matsym[3][41]=1;
      matsym[4][41]=-1;
      matsym[5][41]=-1;
      matsym[0][41]=2;
      matsym[1][41]=3;
      matsym[2][41]=1;
      matsym[3][42]=-1;
      matsym[4][42]=-1;
      matsym[5][42]=-1;
      matsym[0][42]= 2;
      matsym[1][42]= 3;
      matsym[2][42]= 1;
      matsym[3][43]=-1;
      matsym[4][43]= 1;
      matsym[5][43]=-1;
      matsym[0][43]= 2;
      matsym[1][43]= 3;
      matsym[2][43]= 1;
      matsym[3][44]=1;
      matsym[4][44]=1;
      matsym[5][44]=-1;
      matsym[0][44]=3;
      matsym[1][44]=2;
      matsym[2][44]=1;
      matsym[3][45]=1;
      matsym[4][45]=-1;
      matsym[5][45]=-1;
      matsym[0][45]=3;
      matsym[1][45]=2;
      matsym[2][45]=1;
      matsym[3][46]=-1;
      matsym[4][46]=-1;
      matsym[5][46]=-1;
      matsym[0][46]= 3;
      matsym[1][46]= 2;
      matsym[2][46]= 1;
      matsym[3][47]=-1;
      matsym[4][47]= 1;
      matsym[5][47]=-1;
      matsym[0][47]= 3;
      matsym[1][47]= 2;
      matsym[2][47]= 1;

    for(isym=0;isym<48;isym++)
    {
        //DO isym=1,48
        indmat[matsym[3][isym]+1][matsym[4][isym]+1][matsym[5][isym]+1]
            [matsym[0][isym]-1][matsym[1][isym]-1][matsym[2][isym]-1]
            =isym;
    }

    //end of SETMATSYM
    return;
    //END
}


//rewrite by dugang DEC.14
void Band::BUILDPHSCATT(void)    //(cscat)
{
    //____calculate phonon scattering rates

    //CHARACTER*(CPCVL) cscat
    int  itab,iscat,it;
    int  itabl,itabu,itabfp;
    int  iinital,ifinal;
    double k1em,k2em,k3em,k4em,k5em,k6em ;
    double k1ab,k2ab,k3ab,k4ab,k5ab,k6ab,
        kinth,kemh,kabh,besh,
        kintoe,kemoe,kaboe,besoe,
        bes[6],dosd1,dosd2,kint,
        eeafter,sii,
        kopemlow,kopablow,kopemhigh,kopabhigh,
        kaptemlow,kaptablow,kaptemhigh,kaptabhigh,
        kaplemlow,kaplablow,kaplemhigh,kaplabhigh;
    //double CALDOS,CALDOSSUM,HIIRATE,
    double dosfac;
    double ee;
    int ipt;

    if(jacophfl)
    {
        //_______Bose-Einstein for optical phonons
        bes[0] =1.0/(exp(temptag)-1.0);
        bes[1] =1.0/(exp(templag)-1.0);
        bes[2] =1.0/(exp(templog)-1.0);
        bes[3] =1.0/(exp(temptaf)-1.0);
        bes[4] =1.0/(exp(templaf)-1.0);
        bes[5] =1.0/(exp(temptof)-1.0);

        //_______optical phonon constants for emission und absorbtion(electrons)
        //     (phonons+II)
        scpre=14;
        kint=2.0*PI*dfelast*dfelast/sirho/6e0/(siul*siul);

        k1em=(dftag)*(dftag)*PI/temptag/sirho/6e0;
        k2em=(dflag)*(dflag)*PI/templag/sirho/6e0;
        k3em=(dflog)*(dflog)*PI/templog/sirho/6e0;
        k4em=(dftaf)*(dftaf)*PI/temptaf/sirho*0.666e0;
        k5em=(dflaf)*(dflaf)*PI/templaf/sirho*0.666e0;
        k6em=(dftof)*(dftof)*PI/temptof/sirho*0.666e0;

        k1ab=bes[0]*k1em;
        k2ab=bes[1]*k2em;
        k3ab=bes[2]*k3em;
        k4ab=bes[3]*k4em;
        k5ab=bes[4]*k5em;
        k6ab=bes[5]*k6em;

        k1em=(bes[0]+1.0)*k1em;
        k2em=(bes[1]+1.0)*k2em;
        k3em=(bes[2]+1.0)*k3em;
        k4em=(bes[3]+1.0)*k4em;
        k5em=(bes[4]+1.0)*k5em;
        k6em=(bes[5]+1.0)*k6em;

        //_______calculate scattering-rates(electrons)
        ipt=PELEC;
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
        {
            //DO ifinal  =1,nband[ipt]
            for(itab=0;itab<=MTAB;itab++)
            {                    //DO itab=0,MTAB
                ee=energy[itab];

                //__________intra-valley elastic Phonon
                dose[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);

                //__________transversal acoustic g-Phonon(Absorption)
                dose[1][ifinal][itab]=CALDOS(ee+temptag,ifinal+bandof[ipt]);

                //__________transversal acoustic g-Phonon(Emission)
                dose[2][ifinal][itab]=CALDOS(ee-temptag,ifinal+bandof[ipt]);

                //__________longitudinal acoustic g-Phonon(Absorption)
                dose[3][ifinal][itab]=CALDOS(ee+templag,ifinal+bandof[ipt]);

                //__________longitudinal acoustic g-Phonon(Emission)
                dose[4][ifinal][itab]=CALDOS(ee-templag,ifinal+bandof[ipt]);

                //__________longitudinal optical g-Phonon(Absorption)
                dose[5][ifinal][itab]=CALDOS(ee+templog,ifinal+bandof[ipt]);

                //__________longitudinal optical g-Phonon(Emission)
                dose[6][ifinal][itab]=CALDOS(ee-templog,ifinal+bandof[ipt]);

                //__________transversal acoustic f-Phonon(Absorption)
                dose[7][ifinal][itab]=CALDOS(ee+temptaf,ifinal+bandof[ipt]);

                //__________transversal acoustic f-Phonon(Emission)
                dose[8][ifinal][itab]=CALDOS(ee-temptaf,ifinal+bandof[ipt]);

                //__________longitudinal acoustic f-Phonon(Absorption)
                dose[9][ifinal][itab]=CALDOS(ee+templaf,ifinal+bandof[ipt]);

                //__________longitudinal acoustic f-Phonon(Emission)
                dose[10][ifinal][itab]=CALDOS(ee-templaf,ifinal+bandof[ipt]);

                //__________transversal optical f-Phonon(Absorption)
                dose[11][ifinal][itab]=CALDOS(ee+temptof,ifinal+bandof[ipt]);

                //__________transversal optical f-Phonon(Emission)
                dose[12][ifinal][itab]=CALDOS(ee-temptof,ifinal+bandof[ipt]);

            }                    //ENDDO
        }                        //ENDDO

        for(iinital=0;iinital<nband[ipt];iinital++)
        {                        //DO iinital =1,nband[ipt]
            for(ifinal=0;ifinal<nband[ipt];ifinal++)
            {                    //DO ifinal  =1,nband[ipt]

                                 ///be caareful not 1 but 0
                if(iinital==0&&ifinal==0)
                    dosfac=1.0;
                else
                    dosfac=ephb;

                //__________intra-valley elastic Phonon
                scatte[0][ifinal][iinital]=kint*dosfac;

                //__________transversal acoustic g-Phonon(Absorption)
                scatte[1][ifinal][iinital]=k1ab*dosfac;

                //__________transversal acoustic g-Phonon(Emission)
                scatte[2][ifinal][iinital]=k1em*dosfac;

                //__________longitudinal acoustic g-Phonon(Absorption)
                scatte[3][ifinal][iinital]=k2ab*dosfac;

                //__________longitudinal acoustic g-Phonon(Emission)
                scatte[4][ifinal][iinital]=k2em*dosfac;

                //__________longitudinal optical g-Phonon(Absorption)
                scatte[5][ifinal][iinital]=k3ab*dosfac;

                //__________longitudinal optical g-Phonon(Emission)
                scatte[6][ifinal][iinital]=k3em*dosfac;

                //__________transversal acoustic f-Phonon(Absorption)
                scatte[7][ifinal][iinital]=k4ab*dosfac;

                //__________transversal acoustic f-Phonon(Emission)
                scatte[8][ifinal][iinital]=k4em*dosfac;

                //__________longitudinal acoustic f-Phonon(Absorption)
                scatte[9][ifinal][iinital]=k5ab*dosfac;

                //__________longitudinal acoustic f-Phonon(Emission)
                scatte[10][ifinal][iinital]=k5em*dosfac;

                //__________transversal optical f-Phonon(Absorption)
                scatte[11][ifinal][iinital]=k6ab*dosfac;

                //__________transversal optical f-Phonon(Emission)
                scatte[12][ifinal][iinital]=k6em*dosfac;

            }                    //ENDDO
        }                        //ENDDO

        //_______hole constants(phonons+II)
        scprh=4;
        besh =1.0/(exp(temphop)-1.0);
        kinth=2e0*PI*dfhelast*dfhelast/sirho/2.0/(siul*siul);
        kemh=dfhop*dfhop*PI/temphop/sirho;
        kabh=besh*kemh;
        kemh=(besh+1.0)*kemh;

        //_______calculate scattering-rates(holes)
        ipt=PHOLE;
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
        {                        //DO ifinal  =1,nband[ipt]
            for(itab=0;itab<=MTAB;itab++)
            {                    //DO itab=0,MTAB
                ee=energy[itab];

                //__________intra-valley elastic Phonon
                dosh[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);

                //__________(Absorption)
                dosh[1][ifinal][itab]=CALDOS(ee+temphop,ifinal+bandof[ipt]);

                //__________(Emission)
                dosh[2][ifinal][itab]=CALDOS(ee-temphop,ifinal+bandof[ipt]);

            }                    //ENDDO
        }                        //ENDDO

        for(iinital=0;iinital<nband[ipt];iinital++)
        {                        //DO iinital =1,nband[ipt]
            for(ifinal=0;ifinal<nband[ipt];ifinal++)
            {                    //DO ifinal  =1,nband[ipt]

                //__________intra-valley elastic Phonon
                scatth[0][ifinal][iinital]=kinth;

                //__________(Absorption)
                scatth[1][ifinal][iinital]=kabh;

                //__________(Emission)
                scatth[2][ifinal][iinital]=kemh;

            }                    //ENDDO
        }                        //ENDDO
    }                            //if(jcophfl)
    else if(fiscphfl)
    {

        //_______Bose-Einstein for optical phonons
        bes[0] =1.0/(exp(efopee)-1.0);

        //_______seven processes(phonons+II)
        scpre=7;

        //_______optical phonon constants for emission und absorbtion(electrons)
        kopemlow=efoplow*efoplow*PI/efopee/sirho;
        kopemhigh=efophigh*efophigh*PI/efopee/sirho;

        kopablow =bes[0]*kopemlow;
        kopabhigh=bes[0]*kopemhigh;
        kopemlow =(bes[0]+1.0)*kopemlow;
        kopemhigh=(bes[0]+1.0)*kopemhigh;

        //_______acoustic phonon constants for emission und absorbtion(electrons)
        //      including 1 TA and 2 LA modes(LA and TA modes are exchanged)

        kaptablow =35e0*efaplow*efaplow*PI/(efapeet*efapeet*sirho*sia0*sia0);
        kaptabhigh=35e0*(efaphigh*efaphigh)*PI/(efapeet*efapeet)/sirho/(sia0*sia0);

        kaptemlow =(efaplow*efaplow)*(PI*PI*PI)/efapeet/sirho/(sia0*sia0)
            *(1.0/(exp(efapeet)-1.0)+1.0)*6e0;
        kaptemhigh=(efaphigh*efaphigh)*(PI*PI*PI)/efapeet/sirho/(sia0*sia0)
            *(1.0/(exp(efapeet)-1.0)+1.0)*6e0;

        kaplablow =35e0*(efaplow*efaplow)*PI/(efapeel*efapeel)/sirho/(sia0*sia0)
            *2e0;
        kaplabhigh=35e0*(efaphigh*efaphigh)*PI/(efapeel*efapeel)/sirho/(sia0*sia0)
            *2e0;

        kaplemlow =(efaplow*efaplow)*(PI*PI*PI)/efapeel/sirho/(sia0*sia0)
            *(1.0/(exp(efapeel)-1.0)+1.0)*6e0
            *2e0;
        kaplemhigh=(efaphigh*efaphigh)*(PI*PI*PI)/efapeel/sirho/(sia0*sia0)
            *(1.0/(exp(efapeel)-1.0)+1.0)*6e0
            *2e0;

        //_______calculate scattering-rates(electrons)
        ipt=PELEC;
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
        {                        //DO ifinal  =1,nband[ipt]
            for(itab=0;itab<=MTAB;itab++)
            {                    //DO itab=0,MTAB
                ee=energy[itab];
                                 //+ 1 be careful!!!!!
                itabfp=(int)(ee*dlistfp);
                if(itabfp>MWLEFP) itabfp=MWLEFP;

                //__________optical Phonon(Absorption)
                dose[0][ifinal][itab]=CALDOS(ee+efopee,ifinal+bandof[ipt]);

                //__________optical Phonon(Emission)
                dose[1][ifinal][itab]=CALDOS(ee-efopee,ifinal+bandof[ipt]);

                //__________transversal acoustic  Phonon(Absorption)
                dose[2][ifinal][itab]=dosfpab[ifinal+bandof[ipt]][itabfp];

                //__________transversal acoustic  Phonon(Emission)
                dose[3][ifinal][itab]=dosfpem[ifinal+bandof[ipt]][itabfp];

                //__________longitudinal acoustic  Phonon(Absorption)
                dose[4][ifinal][itab]=dosfpab[ifinal+bandof[ipt]][itabfp];

                //__________longitudinal acoustic  Phonon(Emission)
                dose[5][ifinal][itab]=dosfpem[ifinal+bandof[ipt]][itabfp];

            }                    //ENDDO
        }                        //ENDDO

        for(iinital=0;iinital<nband[ipt];iinital++)
        {                        //DO iinital =1,nband[ipt]
            for(ifinal=0;ifinal<nband[ipt];ifinal++)
            {                    //DO ifinal  =1,nband[ipt]

                                 ///be caredul
                if(iinital==0&&ifinal==0)
                {
                    //_____________optical Phonon(Absorption)
                    scatte[0][ifinal][iinital]=kopablow;

                    //_____________optical Phonon(Emission)
                    scatte[1][ifinal][iinital]=kopemlow;

                    //_____________transversal acoustic Phonon(Absorption)
                    scatte[2][ifinal][iinital]=kaptablow;

                    //_____________transversal acoustic Phonon(Emission)
                    scatte[3][ifinal][iinital]=kaptemlow;

                    //_____________longitudinal acoustic Phonon(Absorption)
                    scatte[4][ifinal][iinital]=kaplablow;

                    //_____________longitudinal acoustic Phonon(Emission)
                    scatte[5][ifinal][iinital]=kaplemlow;

                }
                else
                {

                    //_____________transversal optical f-Phonon(Absorption)
                    scatte[0][ifinal][iinital]=kopabhigh;

                    //_____________transversal optical f-Phonon(Emission)
                    scatte[1][ifinal][iinital]=kopemhigh;

                    //_____________transversal acoustic Phonon(Absorption)
                    scatte[2][ifinal][iinital]=kaptabhigh;

                    //_____________transversal acoustic Phonon(Emission)
                    scatte[3][ifinal][iinital]=kaptemhigh;

                    //_____________longitudinal acoustic Phonon(Absorption)
                    scatte[4][ifinal][iinital]=kaplabhigh;

                    //_____________longitudinal acoustic Phonon(Emission)
                    scatte[5][ifinal][iinital]=kaplemhigh;

                }                //endif

            }                    //ENDDO
        }                        //ENDDO

        //_______hole constants(phonons+II)
        scprh=4;
        besh =1.0/(exp(temphop)-1.0);
        kinth=2e0*PI*(dfhelast*dfhelast)/sirho/2e0/(siul*siul);
        kemh=(dfhop*dfhop)*PI/temphop/sirho;
        kabh=besh*kemh;
        kemh=(besh+1.0)*kemh;

        //_______calculate scattering-rates(holes)
        ipt=PHOLE;
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
        {                        //DO ifinal  =1,nband[ipt]
            for(itab=0;itab<=MTAB;itab++)
            {                    //DO itab=0,MTAB
                ee=energy[itab];

                //__________intra-valley elastic Phonon
                dosh[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);

                //__________(Absorption)
                dosh[1][ifinal][itab]=CALDOS(ee+temphop,ifinal+bandof[ipt]);

                //__________(Emission)
                dosh[2][ifinal][itab]=CALDOS(ee-temphop,ifinal+bandof[ipt]);

            }                    //ENDDO
        }                        //ENDDO

        for(iinital=0;iinital<nband[ipt];iinital++)
        {                        //DO iinital =1,nband[ipt]
            for(ifinal=0;ifinal<nband[ipt];ifinal++)
            {                    //DO ifinal  =1,nband[ipt]

                //__________intra-valley elastic Phonon
                scatth[0][ifinal][iinital]=kinth;

                //__________(Absorption)
                scatth[1][ifinal][iinital]=kabh;

                //__________(Emission)
                scatth[2][ifinal][iinital]=kemh;

            }                    //ENDDO
        }                        //ENDDO

    }
    else
    {
        /*
                 WRITE(LUTTYO,*) 'BUILDPHSCAT: Specify phonon system'
                 WRITE(LUOUT ,*) ''
                 STOP
                */
        cout<<"error BUILDPHSCAT: Specify phonon system";exit(0);
    }                            //endif

    //____II-Rate
    //____calculate scattering-rates(electrons)
    ipt=PELEC;
    for(iinital=0;iinital<nband[ipt];iinital++)
    {                            //DO iinital =1,nband[ipt]
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
        {                        //DO ifinal  =1,nband[ipt]

            for(itab=0;itab<=MTAB;itab++)
            {                    //DO itab=0,MTAB
                ee=energy[itab];
                if(iifl)
                    sii=EIIRATE(ee);
                else
                    sii=0.0;
                //     ..calculate final energy after scattering
                //     .. formula by Taniguchi
                eeafter=(ee-sieg)/3.0;
                if(eeafter>0.0)
                {
                    dosd1=CALDOS(eeafter,ifinal+bandof[ipt]);
                    dosd2=CALDOSSUM(eeafter,ipt);
                    if(dosd2==0.0)
                    {
                        dosd1=0.0;
                        dosd2=1.0;
                    }            //endif
                }
                else
                {
                    dosd1=0.0;
                    dosd2=1.0;
                }                //endif
                scattiie[ifinal][iinital][itab]=dosd1/dosd2*sii;
            }                    //ENDDO

        }                        //ENDDO
    }                            //ENDDO

    //____II-Rate
    //____calculate scattering-rates(holes)
    ipt=PHOLE;

    for(iinital=0;iinital<nband[ipt];iinital++)
    {                            //DO iinital =1,nband[ipt]
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
        {                        //DO ifinal  =1,nband[ipt]
            for(itab=0;itab<=MTAB;itab++)
            {                    //DO itab=0,MTAB
                ee=energy[itab];
                if(iifl)
                    sii=HIIRATE(ee);
                else
                    sii=0.0;
                //     .. formula by Taniguchi
                eeafter=(ee-sieg)/3.0;
                if(eeafter>0.0)
                {
                    dosd1=CALDOS(eeafter,ifinal+bandof[ipt]);
                    dosd2=CALDOSSUM(eeafter,ipt);
                    if(dosd2==0e0)
                    {
                        dosd1=0.0;
                        dosd2=1.0;
                    }            //endif
                }
                else
                {
                    dosd1=0.0;
                    dosd2=1.0;
                }                //endif
                scattiih[ifinal][iinital][itab]=dosd1/dosd2*sii;
            }                    //ENDDO

        }                        //ENDDO
    }                            //ENDDO

    //____oxide electron constants
    scproe=3;
    besoe =1.0/(exp(tempoeop)-1.0);
    kintoe=2.0*PI*(dfoeelast*dfoeelast )/sirho/2e0/(siul*siul);
    kemoe=(dfoeop*dfoeop)*1.5e0*PI/tempoeop/sirho;
    kaboe=besoe*kemoe;
    kemoe=(besoe+1.0)*kemoe;

    //____calculate scattering-rates(oxide electrons)
    ipt=POXEL;
    for(ifinal=0;ifinal<nband[ipt];ifinal++)
    {                            //DO ifinal  =1,nband[ipt]
        for(itab=0;itab<=MTAB;itab++)
        {                        //DO itab=0,MTAB
            ee=energy[itab];

            //_______intra-valley elastic Phonon
            dosoe[0][ifinal][itab]=CALDOS(ee,ifinal+bandof[ipt]);

            //_______(Absorption)
            dosoe[1][ifinal][itab]=CALDOS(ee+tempoeop,ifinal+bandof[ipt]);

            //_______(Emission)
            dosoe[2][ifinal][itab]=CALDOS(ee-tempoeop,ifinal+bandof[ipt]);

        }                        //ENDDO
    }                            //ENDDO

    for(iinital=0;iinital<nband[ipt];iinital++)
    {                            //DO iinital =1,nband[ipt]
        for(ifinal=0;ifinal<nband[ipt];ifinal++)
        {                        //DO ifinal  =1,nband[ipt]

            //_______intra-valley elastic Phonon
            scattoe[0][ifinal][iinital]=kintoe;

            //_______(Absorption)
            scattoe[1][ifinal][iinital]=kaboe;

            //_______(Emission)
            scattoe[2][ifinal][iinital]=kemoe;

        }                        //ENDDO
    }                            //ENDDO

    //____Sum of scattering rates
    for(iinital=0;iinital<MNB;iinital++)
    {                            //DO iinital =1,MNB
        for(itab=0;itab<=MTAB;itab++)
        {                        //DO itab=0,MTAB
            sumscatt[itab][iinital]=0.0;
        }                        //ENDDO
    }                            //ENDDO

    //  ..electrons
    ipt=PELEC;
    for(itab=0;itab<=MTAB;itab++)
    {                            //DO itab=0,MTAB
        for(iinital=0;iinital<nband[ipt];iinital++)
        {                        //DO iinital =1,nband[ipt]
            for(ifinal=0;ifinal<nband[ipt];ifinal++)
            {                    //DO ifinal  =1,nband[ipt]
                for( iscat=0;iscat<scpre-1;iscat++)
                {                //DO iscat=1,scpre-1
                    sumscatt[itab][iinital+bandof[ipt]]=
                        sumscatt[itab][iinital+bandof[ipt]]
                        +scatte[iscat][ifinal][iinital]
                        *dose[iscat][ifinal][itab];
                }                //ENDDO
                sumscatt[itab][iinital+bandof[ipt]]=
                    sumscatt[itab][iinital+bandof[ipt]]
                    +scattiie[ifinal][iinital][itab];
            }                    //ENDDO
        }                        //ENDDO
    }                            //ENDDO

    //  ..holes
    ipt=PHOLE;
    for(itab=0;itab<=MTAB;itab++)
    {                            //DO itab=0,MTAB
        for(iinital=0;iinital<nband[ipt];iinital++)
        {                        //DO iinital =1,nband[ipt]
            for(ifinal=0;ifinal<nband[ipt];ifinal++)
            {                    //DO ifinal  =1,nband[ipt]
                for(iscat=0;iscat<scprh-1;iscat++)
                {                //DO iscat=1,scprh-1
                    sumscatt[itab][iinital+bandof[ipt]]=
                        sumscatt[itab][iinital+bandof[ipt]]
                        +scatth[iscat][ifinal][iinital]
                        *dosh[iscat][ifinal][itab];
                }                //ENDDO
                sumscatt[itab][iinital+bandof[ipt]]=
                    sumscatt[itab][iinital+bandof[ipt]]
                    +scattiih[ifinal][iinital][itab];
            }                    //ENDDO
        }                        //ENDDO
    }                            //ENDDO

    //  ..oxide electrons
    ipt=POXEL;
    for(itab=0;itab<=MTAB;itab++)
    {                            //DO itab=0,MTAB
        for(iinital=0;iinital<nband[ipt];iinital++)
        {                        //DO iinital =1,nband[ipt]
            for(ifinal=0;ifinal<nband[ipt];ifinal++)
            {                    //DO ifinal  =1,nband[ipt]
                for(iscat=0;iscat<scproe-1;iscat++)
                {                //DO iscat=1,scproe
                    sumscatt[itab][iinital+bandof[ipt]]=
                        sumscatt[itab][iinital+bandof[ipt]]
                        +scattoe[iscat][ifinal][iinital]
                        *dosoe[iscat][ifinal][itab];
                }                //ENDDO
            }                    //ENDDO
        }                        //ENDDO
    }                            //ENDDO

    //____cal. maximum scattering rate in each tetrahedron
    for(it=0;it<nt;it++)
    {                            //DO it=1,nt
        ipt=partyp[ibt[it]];
        itabl=int(Max(0,(int)(eek[tet[0][it]]/dtable)));
        itabu=int(Min(MTAB,(int)(eek[tet[3][it]]/dtable)+1));
        gamtet[it]=0.0;
        for(itab=itabl;itab<=itabu;itab++)
        {                        //DO itab=itabl,itabu
            gamtet[it]=Max(gamtet[it],sumscatt[itab][ibt[it]]);
        }                        //ENDDO
        if(gamtet[it]==0.0) gamtet[it]=sumscatt[MTAB][ibt[it]];
    }                            //ENDDO

    //____End of BUILDPHSCATT
    return;
}


////////havearest///////////////
//====

double Band::EIIRATE(double ee)
{
    //____Purpose : cal. electron II-rate
    double sii,eel,rtnv;

    //____account for temperature dependent band gap
    eel=ee-sieg+1.1241e0/eV0;
    sii=0.0;

    if(kaneiifl)
    {
        //__Kane
        if(eel*eV0>1.1241e0)
        {
            sii=7.216e11/scrt0*pow((eel*eV0-1.12410),3.540);
        }                        //endif
    }                            //endif

    if(thomiifl)
        //__Thoma
        if((eel*eV0>1.128e0)&&(eel*eV0<1.75e0))
    {
        sii=1.25e12/scrt0*pow((eel*eV0-1.128e0),3);
    }
    else if(eel*eV0>1.75e0)
    {
        sii=9.49e12/scrt0*(eel*eV0-1.572e0)*(eel*eV0-1.572e0);
    }                            //endif

    if(sanoiifl)
        //__sano
        if(eel*eV0>1.242e0)
    {
        sii=4.25e11/scrt0*pow((eel*eV0-1.242e0),4.188);
    }                            //endif

    if(fisciifl)
    {
        //___Fischetti
        if(eel*eV0>1.2e0)
        {
            sii=6.25e10/scrt0*pow(((eel*eV0-1.2e0)/1.20),2);
        }                        //endif
        if(eel*eV0>1.8e0)
        {
            sii=3e12/scrt0*pow(((eel*eV0-1.80)/1.80),2)+sii;
        }                        //endif
        if(eel*eV0>3.45e0)
        {
            sii=6.8e14/scrt0*pow(((eel*eV0-3.45e0)/3.45e0),2)+sii;
        }                        //endif
    }                            //endif

    rtnv=sii*iifacelec;

    //____End of EIIRATE
    return rtnv;
}


//====

double Band::HIIRATE(double ee)
{

    //____Purpose : cal. hole II-rate
    //____local variable
    double eel,rtnv;

    //____account for temperature dependent band gap
    eel=ee-sieg+1.1241e0/eV0;

    if(eel>hiithresh)
    {
        rtnv=iifachole*pow((eel-hiithresh),hiiexp);
    }
    else
    {
        rtnv=0.0;
    }                            //endif

    //____End of HIIRATE
    return rtnv;
}


double Band::CALSCATTSUM(double eed,int iinital)
{

    //____Purpose : interpolate the dos-rate between
    //             two tab-values
    //   Parameter : eed=energy of actual electron
    //             : ip=particle type

    int itab;
    double intp,rtnv;

    //____smallest allowed energy
    if(eed<=0.0)
    {
        rtnv=sumscatt[0][iinital];
        //____largest allowed energy
    }
    else if(eed>=emax)
    {
        rtnv=sumscatt[MTAB][iinital];
    }
    else
    {
        itab=(int)(eed/dtable);
        intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);
        rtnv=sumscatt[itab][iinital]+intp
            *(sumscatt[itab+1][iinital]
            -sumscatt[itab][iinital]);
    }                            //endif

    //____End of CALSCATTSUM
    return rtnv;
}


//====

void Band::CALSCATTOE(double eed,double*sca,int ib)
{

    //____Purpose : interpolate the scattering rates between
    //             two tab-values for oxide electrons
    //   Parameter : eed=energy of actual hole
    //             : sca=array of scattering rates for eed

    int itab,is,iinital,ifinal;
    double intp;                 //,sca(SCPROE*NBOE)

    //
    iinital=ib-bandof[POXEL];

    //____smallest allowed energy
    if(eed<=0e0)
    {
        for(ifinal=0;ifinal<nband[POXEL];ifinal++)
        {                        ////DO ifinal=1,nband[POXEL]
            for( is=0;is<scproe;is++)
            {                    //DO is=1,scproe
                sca[is+scproe*(ifinal)]=scattoe[is][ifinal][iinital]
                    *dosoe[is][ifinal][0];
            }                    //ENDDO
        }                        //ENDDO
        //____largest allowed energy
    }
    else if(eed>=emax)
    {
        for(ifinal=0;ifinal<nband[POXEL];ifinal++)
        {                        ////DO ifinal=1,nband[POXEL]
            for( is=0;is<scproe;is++)
            {                    //DO is=1,scproe
                sca[is+scproe*(ifinal)]=scattoe[is][ifinal][iinital]
                    *dosoe[is][ifinal][MTAB];
            }                    //ENDDO
        }                        //ENDDO
    }
    else
    {
        itab=(int)(eed/dtable);
        intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);

        for(ifinal=0;ifinal<nband[POXEL];ifinal++)
        {                        ////DO ifinal=1,nband[POXEL]
            for( is=0;is<scproe;is++)
            {                    //DO is=1,scproe
                sca[is+scproe*(ifinal-1)]=scattoe[is][ifinal][iinital]
                    *(dosoe[is][ifinal][itab]+intp
                    *(dosoe[is][ifinal][itab+1]
                    -dosoe[is][ifinal][itab]));
            }                    //ENDDO
        }                        //ENDDO
    }                            //endif

    //____End of CALSCATTOE
    return;
}


//====

void Band::CALSCATTH(double eed,double*sca,int ib)
{
    //____Purpose : interpolate the scattering rates between
    //             two tab-values for holes
    //   Parameter : eed=energy of actual hole
    //             : sca=array of scattering rates for eed

    int itab,is,iinital,ifinal;
    double intp;                 //,sca(SCPRH*NBH)

    //
    iinital=ib-bandof[PHOLE];

    //____smallest allowed energy
    if(eed<=0.0)
    {
        for(ifinal=0;ifinal<nband[PHOLE];ifinal++)
        {                        //         DO ifinal=1,nband[PHOLE]
            for(is=0;is<scprh-1;is++)
            {
                sca[is+scprh*(ifinal)]=scatth[is][ifinal][iinital]
                    *dosh[is][ifinal][0];
            }                    //ENDDO
            sca[scprh-1+scprh*(ifinal)]=scattiih[ifinal][iinital][0];
        }                        //ENDDO
        //____largest allowed energy
    }
    else if(eed>=emax)
    {
        for(ifinal=0;ifinal<nband[PHOLE];ifinal++)
        {                        //         DO ifinal=1,nband[PHOLE]
            for(is=0;is<scprh-1;is++)
            {                    //DO is=1,scprh-1
                sca[is+scprh*(ifinal)]=scatth[is][ifinal][iinital]
                    *dosh[is][ifinal][MTAB];
            }                    //ENDDO
            sca[scprh-1+scprh*(ifinal)]=scattiih[ifinal][iinital][MTAB];
        }                        //ENDDO
    }
    else
    {
        itab=(int)(eed/dtable);
        intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);

        for(ifinal=0;ifinal<nband[PHOLE];ifinal++)
        {                        //         DO ifinal=1,nband[PHOLE]
            for(is=0;is<scprh-1;is++)
            {                    //DO is=1,scprh-1
                sca[is+scprh*(ifinal)]=scatth[is][ifinal][iinital]
                    *(dosh[is][ifinal][itab]+intp
                    *(dosh[is][ifinal][itab+1]
                    -dosh[is][ifinal][itab]));
            }                    //ENDDO
            sca[scprh-1+scprh*(ifinal)]=
                scattiih[ifinal][iinital][itab]+intp
                *(scattiih[ifinal][iinital][itab+1]
                -scattiih[ifinal][iinital][itab]);
        }                        //ENDDO
    }                            //endif

    //____End of CALSCATTH
    return;
}


//====

void Band::CALSCATTE(double eed,double*sca,int ib)
{

    //____Purpose : interpolate the scattering rates between
    //             two tab-values for electrons
    //   Parameter : eed=energy of actual electron
    //             : sca=array of scattering rates for eed
    int itab,is,iinital,ifinal;
    double intp;                 //,sca(MNScaEle*NBE)

    //
    iinital=ib-bandof[PELEC];

    //____smallest allowed energy
    if(eed<=0e0)
    {
        for(ifinal=0;ifinal<nband[PELEC];ifinal++)
        {                        //DO ifinal=1,nband[PELEC]
            for(is=0;is<scpre-1;is++)
            {                    //DO is=1,scpre-1
                sca[is+scpre*(ifinal)]=scatte[is][ifinal][iinital]
                    *dose[is][ifinal][0];
            }                    //ENDDO
            sca[scpre-1+scpre*(ifinal)]=scattiie[ifinal][iinital][0];
        }                        //ENDDO
        //____largest allowed energy
    }
    else if(eed>=emax)
    {
        for(ifinal=0;ifinal<nband[PELEC];ifinal++)
        {                        //DO ifinal=1,nband[PELEC]
            for(is=0;is<scpre-1;is++)
            {                    //DO is=1,scpre-1
                sca[is+scpre*(ifinal)]=scatte[is][ifinal][iinital]
                    *dose[is][ifinal][MTAB];
            }                    //ENDDO
            sca[scpre-1+scpre*(ifinal)]=scattiie[ifinal][iinital][MTAB];
        }                        //ENDDO
    }
    else
    {
        itab=(int)(eed/dtable);
        intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);

        for(ifinal=0;ifinal<nband[PELEC];ifinal++)
        {                        //DO ifinal=1,nband[PELEC]
            for(is=0;is<scpre-1;is++)
            {                    //DO is=1,scpre-1
                sca[is+scpre*(ifinal)]=scatte[is][ifinal][iinital]
                    *(dose[is][ifinal][itab]+intp
                    *(dose[is][ifinal][itab+1]
                    -dose[is][ifinal][itab]));
            }                    //ENDDO
            sca[scpre-1+scpre*(ifinal)]=
                scattiie[ifinal][iinital][itab]+intp
                *(scattiie[ifinal][iinital][itab+1]
                -scattiie[ifinal][iinital][itab]);
        }                        //ENDDO
    }                            //endif

    //____End of CALSCATTE
    return;
}


//====

double Band::CALDOS(double eed,int ib)
{
    //____Purpose : interpolate the dos between
    //             two tab-values
    //   Parameter : eed=energy of actual electron
    //             : ib=bandindex
    int itab;
    double intp,rtnv;

    //____smallest allowed energy
    if(eed<=0e0)
    {
        rtnv=dos[ib][0];
        //____largest allowed energy
    }
    else if(eed>=emax)
    {
        rtnv=dos[ib][MTAB];
    }
    else
    {
        itab=(int)(eed/dtable);
        intp=(eed-energy[itab])/(energy[itab+1]-energy[itab]);
        rtnv=dos[ib][itab]+intp*(dos[ib][itab+1]-dos[ib][itab]);
    }                            //endif

    //____End of CALDOS
    return rtnv;
}


//====

double Band::CALDOSSUM(double eed,int ip)
{

    //____Purpose : interpolate the total dos between
    //             two tab-values
    //   Parameter : eed=energy of actual electron
    //             : ip=particle type

    int itab;
    double rtnv,intp;
    //____smallest allowed energy
    if(eed<=0e0)
        rtnv=sumdos[0][ip];
    //____largest allowed energy
    else if(eed>=emax)
        rtnv=sumdos[MTAB][ip];
    else
    {
        itab=(int)(eed/dtable);
        intp=(eed-energy[itab]) / (energy[itab+1]-energy[itab]);
        rtnv=sumdos[itab][ip]+intp
            *(sumdos[itab+1][ip]-sumdos[itab][ip]);
    }                            //endif

    //____End of CALDOSSUM
    return rtnv;
}


//====

double Band::CALDOSTETMAX(double eed,int ib)
{
    //____Purpose : the maximum area between two tab-values
    //   Parameter : eed=energy of actual electron
    //             : ib=bandindex
    int itab;
    double rtnv;
    //____smallest allowed energy
    if(eed<=0.0)
        rtnv=dostetmax[ib][0];
    //____largest allowed energy
    else if(eed>=emax)
        rtnv=dostetmax[ib][MTAB];
    else
    {
        itab=(int)(eed/dtable);
        rtnv=Max(dostetmax[ib][itab],dostetmax[ib][itab+1]);
    }                            //endif

    //____End of CALDOSTETMAX
    return rtnv;
}


//====

void Band::OVERLAP(int it,int&nct,int*ict)
{
    // Purpose:   Find cubes which contain the given tetrahedron.
    //--------  For each cube a list is built,which gives the tetrahedrons
    //            contained in the cube. The cube grid extends from
    //            0.7<=x<1.0,0<=y<0.075 and 0<=z<0.075 and is made only
    //            for the first conduction band.
    //____local variables
    bool Flag_Find;
    //int it,nct,ict(MCLE),
    int itc,ibase;
    int ic,icx,icy,icz;
    int icxmin,icymin,iczmin;
    int icxmax,icymax,iczmax;
    double xkl,xkh,ykl,ykh,zkl,zkh;
    double d1,d2,d3,d4;
    double ex,ey,ez,el,ll;
    double xs,ys,zs;
    double xk1,xk2,xk3,xk4;
    double yk1,yk2,yk3,yk4;
    double zk1,zk2,zk3,zk4;
    double xkmin,ykmin,zkmin;
    double xkmax,ykmax,zkmax;
    double dnt1,dnt2,dnt3,dnt4;
    double dnt5,dnt6,dnt7,dnt8;
    double dnt9,dnt10,dnt11,dnt12;
    double dnt13,dnt14,dnt15,dnt16;

    //____initializes numbers of cubes within the tetrahedron
    nct=0;

    //____check if tetrahedron is outside of all cubes
    Flag_Find=false;
    for(itc=0;itc<4;itc++)
    {                            //DO itc=1,4
        if(ykk[tet[itc][it]]<0.075*a0pi)Flag_Find=true;
    }                            //ENDDO
    if(!Flag_Find)               //!Flag_Find 1
    {                            //RETURN
        Flag_Find=false;
        for(itc=0;itc<4;itc++)   //DO itc=1,4
        {
            if(xkk[tet[itc][it]]>0.7*a0pi)Flag_Find=true;
        }                        //ENDDO
    }                            //!Flag_Find 1
    if(!Flag_Find)               ///!Flag_Find 2
    {                            //RETURN

        //____
        xk1=xkk[tet[0][it]];
        xk2=xkk[tet[1][it]];
        xk3=xkk[tet[2][it]];
        xk4=xkk[tet[3][it]];

        yk1=ykk[tet[0][it]];
        yk2=ykk[tet[1][it]];
        yk3=ykk[tet[2][it]];
        yk4=ykk[tet[3][it]];

        zk1=zkk[tet[0][it]];
        zk2=zkk[tet[1][it]];
        zk3=zkk[tet[2][it]];
        zk4=zkk[tet[3][it]];

        xkmin=Min(xk1,xk2,xk3,xk4);
        xkmax=Max(xk1,xk2,xk3,xk4);

        ykmin=Min(yk1,yk2,yk3,yk4);
        ykmax=Max(yk1,yk2,yk3,yk4);

        zkmin=Min(zk1,zk2,zk3,zk4);
        zkmax=Max(zk1,zk2,zk3,zk4);

        icxmin=int(Max(0,(int)((xkmin/a0pi-0.7e0)/0.3e0*(MNCubeX))));
        icymin=int(Max(0,(int)(ykmin/(a0pi*0.075e0)*(MNCubeY))));
        iczmin=int(Max(0,(int)(zkmin/(a0pi*0.075e0)*(MNCubeZ))));

        icxmax=int(Min(MNCubeX-1,(int)((xkmax/a0pi-0.7e0)/0.3e0*(MNCubeX))));
        icymax=int(Min(MNCubeY-1,(int)(ykmax/(a0pi*0.075e0)*(MNCubeY))));
        iczmax=int(Min(MNCubeZ-1,(int)(zkmax/(a0pi*0.075e0)*(MNCubeZ))));

        //ibase=16*(it-1)
        ibase=it;
        dnt1=datantlin[0][0][it];//(1+ibase)
        dnt2=datantlin[1][0][it];//(2+ibase)
        dnt3=datantlin[2][0][it];//(3+ibase)
        dnt4=datantlin[3][0][it];//(4+ibase)
        dnt5=datantlin[0][1][it];//(5+ibase)
        dnt6=datantlin[1][1][it];//(6+ibase)
        dnt7=datantlin[2][1][it];//(7+ibase)
        dnt8=datantlin[3][1][it];//(8+ibase)
        dnt9=datantlin[0][2][it];//(9+ibase)
                                 //(10+ibase)
        dnt10=datantlin[1][2][it];
                                 //(11+ibase)
        dnt11=datantlin[2][2][it];
                                 //(12+ibase)
        dnt12=datantlin[3][2][it];
                                 //(13+ibase)
        dnt13=datantlin[0][3][it];
                                 //(14+ibase)
        dnt14=datantlin[1][3][it];
                                 //(15+ibase)
        dnt15=datantlin[2][3][it];
                                 //(16+ibase)
        dnt16=datantlin[3][3][it];

        //____check every cube
        for(icx=icxmin;icx<=icxmax;icx++)
        {                        //DO icx=icxmin,icxmax
            for(icy=icymin;icy<=icymax;icy++)
            {                    //DO icy=icymin,icymax
                for(icz=iczmin;icz<=iczmax;icz++)
                {                //DO icz=iczmin,iczmax
                                 //+1 minim=0!!!!!!!!!!!ict[]=ic;ict[] is used as index???????
                    ic=icz+MNCubeZ*(icy+MNCubeY*icx);
                    //_______check if the cube under investigation contains the tetrahedron
                    Flag_Find=false;
                    xkl=(0.3e0/(MNCubeX)*(icx)+0.7e0)*a0pi;
                    if(xk1>xkl)Flag_Find=true;
                    if(xk2>xkl)Flag_Find=true;
                    if(xk3>xkl)Flag_Find=true;
                    if(xk4>xkl)Flag_Find=true;
                    if(!Flag_Find)continue;

                    Flag_Find=false;
                    xkh=(0.3e0/(MNCubeX)*(icx+1)+0.7e0)*a0pi;
                    if(xk1<xkh)Flag_Find=true;
                    if(xk2<xkh)Flag_Find=true;
                    if(xk3<xkh)Flag_Find=true;
                    if(xk4<xkh)Flag_Find=true;
                    if(!Flag_Find)continue;

                    Flag_Find=false;
                    ykl=(0.075e0/(MNCubeY)*(icy))*a0pi;
                    if(yk1>ykl)Flag_Find=true;
                    if(yk2>ykl)Flag_Find=true;
                    if(yk3>ykl)Flag_Find=true;
                    if(yk4>ykl)Flag_Find=true;
                    if(!Flag_Find)continue;

                    Flag_Find=false;
                    ykh=(0.075e0/(MNCubeY)*(icy+1))*a0pi;
                    if(yk1<ykh)Flag_Find=true;
                    if(yk2<ykh)Flag_Find=true;
                    if(yk3<ykh)Flag_Find=true;
                    if(yk4<ykh)Flag_Find=true;
                    if(!Flag_Find)continue;

                    Flag_Find=false;
                    zkl=(0.075e0/(MNCubeZ)*(icz))*a0pi;
                    if(zk1>zkl)Flag_Find=true;
                    if(zk2>zkl)Flag_Find=true;
                    if(zk3>zkl)Flag_Find=true;
                    if(zk4>zkl)Flag_Find=true;
                    if(!Flag_Find)continue;

                    Flag_Find=false;
                    zkh=(0.075e0/(MNCubeZ)*(icz+1))*a0pi;
                    if(zk1<zkh)Flag_Find=true;
                    if(zk2<zkh)Flag_Find=true;
                    if(zk3<zkh)Flag_Find=true;
                    if(zk4<zkh)Flag_Find=true;
                                 //itetcube
                    if(!Flag_Find)continue;

                    //_______check if one of the nodes of the tetrahedron is within the cube
                    Flag_Find=false;
                    if((xk1>=xkl)&&
                        (xk1<=xkh)&&
                        (yk1>=ykl)&&
                        (yk1<=ykh)&&
                        (zk1>=zkl)&&
                        (zk1<=zkh))Flag_Find=true;
                    if((xk2>=xkl)&&
                        (xk2<=xkh)&&
                        (yk2>=ykl)&&
                        (yk2<=ykh)&&
                        (zk2>=zkl)&&
                        (zk2<=zkh))Flag_Find=true;
                    if((xk3>=xkl)&&
                        (xk3<=xkh)&&
                        (yk3>=ykl)&&
                        (yk3<=ykh)&&
                        (zk3>=zkl)&&
                        (zk3<=zkh))Flag_Find=true;
                    if((xk4>=xkl)&&
                        (xk4<=xkh)&&
                        (yk4>=ykl)&&
                        (yk4<=ykh)&&
                        (zk4>=zkl)&&
                        (zk4<=zkh))Flag_Find=true;
                    if(Flag_Find)
                    {
                        nct=nct+1;
                        if(nct>MCLE)
                        {
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }

                    //_______check if the cube is sorounded by a tetrahedron
                    //cubeitet

                    //cube(xkl,ykl,zkl)itet
                    //   ..xkl,ykl,zkl
                    d1=dnt1*xkl+dnt2*ykl
                        +dnt3*zkl-dnt4;
                    d2=dnt5*xkl+dnt6*ykl
                        +dnt7*zkl-dnt8;
                    d3=dnt9*xkl+dnt10*ykl
                        +dnt11*zkl-dnt12;
                    d4=dnt13*xkl+dnt14*ykl
                        +dnt15*zkl-dnt16;
                    if(d1>0.0&&d2>0.0&&
                        d3>0.0&&d4>0.0)
                    {
                        nct=nct+1;
                        if(nct>MCLE)
                        {
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }

                    //   ..xkl,ykl,zkh
                    //cube(xkl,ykl,zkh)itet
                    d1=dnt1*xkl+dnt2*ykl
                        +dnt3*zkh-dnt4;
                    d2=dnt5*xkl+dnt6*ykl
                        +dnt7*zkh-dnt8;
                    d3=dnt9*xkl+dnt10*ykl
                        +dnt11*zkh-dnt12;
                    d4=dnt13*xkl+dnt14*ykl
                        +dnt15*zkh-dnt16;
                    if(d1>0.0&&d2>0.0&&
                        d3>0.0&&d4>0.0)
                    {            //THEN
                        nct=nct+1;
                        if(nct>MCLE)
                        {        //THEN
                            /*
                            WRITE(LUTTYO,*)'OVERLAP: nct >MCLE'
                            WRITE(LUOUT ,*)''
                            STOP
                            */
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }        //ENDIF
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }            //ENDIF

                    //   ..xkl,ykh,zkl
                    d1=dnt1*xkl+dnt2*ykh
                        +dnt3*zkl-dnt4;
                    d2=dnt5*xkl+dnt6*ykh
                        +dnt7*zkl-dnt8;
                    d3=dnt9*xkl+dnt10*ykh
                        +dnt11*zkl-dnt12;
                    d4=dnt13*xkl+dnt14*ykh
                        +dnt15*zkl-dnt16;
                    if(d1>0.0&&d2>0.0&&
                        d3>0.0&&d4>0.0)
                    {            //THEN
                        nct=nct+1;
                        if(nct>MCLE)
                        {
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }        //ENDIF
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }            //ENDIF

                    //   ..xkl,ykh,zkh
                    d1=dnt1*xkl+dnt2*ykh
                        +dnt3*zkh-dnt4;
                    d2=dnt5*xkl+dnt6*ykh
                        +dnt7*zkh-dnt8;
                    d3=dnt9*xkl+dnt10*ykh
                        +dnt11*zkh-dnt12;
                    d4=dnt13*xkl+dnt14*ykh
                        +dnt15*zkh-dnt16;
                    if(d1>0.0&&d2>0.0&&
                        d3>0.0&&d4>0.0)
                    {            //THEN
                        nct=nct+1;
                        if(nct>MCLE)
                        {
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }        //ENDIF
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }            //ENDIF

                    //   ..xkh,ykl,zkl
                    d1=dnt1*xkh+dnt2*ykl
                        +dnt3*zkl-dnt4;
                    d2=dnt5*xkh+dnt6*ykl
                        +dnt7*zkl-dnt8;
                    d3=dnt9*xkh+dnt10*ykl
                        +dnt11*zkl-dnt12;
                    d4=dnt13*xkh+dnt14*ykl
                        +dnt15*zkl-dnt16;
                    if(d1>0.0&&d2>0.0&&
                        d3>0.0&&d4>0.0)
                    {            //THEN
                        nct=nct+1;
                        if(nct>MCLE)
                        {
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }        //ENDIF
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }            //ENDIF

                    //   ..xkh,ykl,zkh
                    d1=dnt1*xkh+dnt2*ykl
                        +dnt3*zkh-dnt4;
                    d2=dnt5*xkh+dnt6*ykl
                        +dnt7*zkh-dnt8;
                    d3=dnt9*xkh+dnt10*ykl
                        +dnt11*zkh-dnt12;
                    d4=dnt13*xkh+dnt14*ykl
                        +dnt15*zkh-dnt16;
                    if(d1>0.0&&d2>0.0&&
                        d3>0.0&&d4>0.0)
                    {            //THEN
                        nct=nct+1;
                        if(nct>MCLE)
                        {
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }        //ENDIF
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }            //ENDIF

                    //   ..xkh,ykh,zkl
                    d1=dnt1*xkh+dnt2*ykh
                        +dnt3*zkl-dnt4;
                    d2=dnt5*xkh+dnt6*ykh
                        +dnt7*zkl-dnt8;
                    d3=dnt9*xkh+dnt10*ykh
                        +dnt11*zkl-dnt12;
                    d4=dnt13*xkh+dnt14*ykh
                        +dnt15*zkl-dnt16;
                    if(d1>0.0&&d2>0.0&&
                        d3>0.0&&d4>0.0)
                    {            //THEN
                        nct=nct+1;
                        if(nct>MCLE)
                        {
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }        //ENDIF
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }            //ENDIF

                    //   ..xkh,ykh,zkh
                    d1=dnt1*xkh+dnt2*ykh
                        +dnt3*zkh-dnt4;
                    d2=dnt5*xkh+dnt6*ykh
                        +dnt7*zkh-dnt8;
                    d3=dnt9*xkh+dnt10*ykh
                        +dnt11*zkh-dnt12;
                    d4=dnt13*xkl+dnt14*ykh
                        +dnt15*zkh-dnt16;
                    if(d1>0.0&&d2>0.0&&
                        d3>0.0&&d4>0.0)
                    {            //THEN
                        nct=nct+1;
                        if(nct>MCLE)
                        {
                            cout<<"error OVERLAP: nct >MCLE";exit(0);
                        }        //ENDIF
                                 //nct is number not inedx so use nct-1;
                        ict[nct-1]=ic;
                        continue;
                    }            //ENDIF

                    //_______intersection between edge 1 of tetrahedron and surface of cube
                    ex=xk2-xk1;
                    ey=yk2-yk1;
                    ez=zk2-zk1;
                    el=sqrt(ex*ex+ey*ey+ez*ez);
                    ex=ex/el;
                    ey=ey/el;
                    ez=ez/el;

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkl-xk1)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk1;
                            zs=ll*ez+zk1;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkh-xk1)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk1;
                            zs=ll*ez+zk1;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykl-yk1)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            zs=ll*ez+zk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykh-yk1)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            zs=ll*ez+zk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF
                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkl-zk1)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            ys=ll*ey+yk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkh-zk1)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            ys=ll*ey+yk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    //_______intersection between edge 2 of tetrahedron and surface of cube
                    ex=xk3-xk1;
                    ey=yk3-yk1;
                    ez=zk3-zk1;
                    el=sqrt(ex*ex+ey*ey+ez*ez);
                    ex=ex/el;
                    ey=ey/el;
                    ez=ez/el;

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkl-xk1)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk1;
                            zs=ll*ez+zk1;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkh-xk1)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk1;
                            zs=ll*ez+zk1;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykl-yk1)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            zs=ll*ez+zk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykh-yk1)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            zs=ll*ez+zk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkl-zk1)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            ys=ll*ey+yk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkh-zk1)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            ys=ll*ey+yk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    //_______intersection between edge 3 of tetrahedron and surface of cube
                    ex=xk4-xk1;
                    ey=yk4-yk1;
                    ez=zk4-zk1;
                    el=sqrt(ex*ex+ey*ey+ez*ez);
                    ex=ex/el;
                    ey=ey/el;
                    ez=ez/el;

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkl-xk1)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk1;
                            zs=ll*ez+zk1;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkh-xk1)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk1;
                            zs=ll*ez+zk1;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykl-yk1)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            zs=ll*ez+zk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykh-yk1)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            zs=ll*ez+zk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkl-zk1)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            ys=ll*ey+yk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkh-zk1)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk1;
                            ys=ll*ey+yk1;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    //_______intersection between edge 4 of tetrahedron and surface of cube
                    ex=xk3-xk2;
                    ey=yk3-yk2;
                    ez=zk3-zk2;
                    el=sqrt(ex*ex+ey*ey+ez*ez);
                    ex=ex/el;
                    ey=ey/el;
                    ez=ez/el;

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkl-xk2)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk2;
                            zs=ll*ez+zk2;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkh-xk2)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk2;
                            zs=ll*ez+zk2;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykl-yk1)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk2;
                            zs=ll*ez+zk2;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykh-yk1)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk2;
                            zs=ll*ez+zk2;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkl-zk2)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk2;
                            ys=ll*ey+yk2;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkh-zk2)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk2;
                            ys=ll*ey+yk2;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    //_______intersection between edge 5 of tetrahedron and surface of cube
                    ex=xk4-xk2;
                    ey=yk4-yk2;
                    ez=zk4-zk2;
                    el=sqrt(ex*ex+ey*ey+ez*ez);
                    ex=ex/el;
                    ey=ey/el;
                    ez=ez/el;

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkl-xk2)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk2;
                            zs=ll*ez+zk2;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkh-xk2)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk2;
                            zs=ll*ez+zk2;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykl-yk2)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk2;
                            zs=ll*ez+zk2;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykh-yk2)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk2;
                            zs=ll*ez+zk2;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkl-zk2)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk2;
                            ys=ll*ey+yk2;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkh-zk2)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk2;
                            ys=ll*ey+yk2;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    //_______intersection between edge 6 of tetrahedron and surface of cube
                    ex=xk4-xk3;
                    ey=yk4-yk3;
                    ez=zk4-zk3;
                    el=sqrt(ex*ex+ey*ey+ez*ez);
                    ex=ex/el;
                    ey=ey/el;
                    ez=ez/el;

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkl-xk3)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk3;
                            zs=ll*ez+zk3;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ex)>0.0)
                    {            //THEN
                        ll=(xkh-xk3)/ex;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            ys=ll*ey+yk3;
                            zs=ll*ez+zk3;
                            if((ys>ykl)&&(ys<ykh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykl-yk3)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk3;
                            zs=ll*ez+zk3;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ey)>0.0)
                    {            //THEN
                        ll=(ykh-yk3)/ey;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk3;
                            zs=ll*ez+zk3;
                            if((xs>xkl)&&(xs<xkh)&&
                                (zs>zkl)&&(zs<zkh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkl-zk3)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk3;
                            ys=ll*ey+yk3;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    if(fabs(ez)>0.0)
                    {            //THEN
                        ll=(zkh-zk3)/ez;
                        if((ll>0.0)&&(ll<1e0))
                        {        //THEN
                            xs=ll*ex+xk3;
                            ys=ll*ey+yk3;
                            if((xs>xkl)&&(xs<xkh)&&
                                (ys>ykl)&&(ys<ykh))
                            {    //THEN
                                nct=nct+1;
                                if(nct>MCLE)
                                {
                                    cout<<"error OVERLAP: nct >MCLE";exit(0);
                                    //ENDIF
                                }
                                 //nct is number not inedx so use nct-1;
                                ict[nct-1]=ic;
                                continue;
                            }    //ENDIF
                        }        //ENDIF
                    }            //ENDIF

                    //10      CONTINUE
                }                //ENDDO
            }                    //ENDDO
        }                        //ENDDO
    }                            //!Flag_Find 2
    //____end of OVERLAP
    return;
}


//====

void Band::GETINDENS(void)
{
    //____Calculate intrinsic carrier density(DOS is only for one spin
    // direction)

    //____local variables
    int itab;
    double en0,ep0,eel,eelmo;    //,CALDOSSUM
    double dosl,doslmo;

    en0=0.0;
    ep0=0.0;

    //____calculate intrinsic carrier densities using Boltzmann statistics
    eelmo=energy[0];
    for(itab=1;itab<=MTAB;itab++)
    {                            //DO itab=1,MTAB
        eel=energy[itab];
        dosl=CALDOSSUM(eel,PELEC);
        doslmo=CALDOSSUM(eelmo,PELEC);
        en0=en0+((doslmo*eel-dosl*eelmo)
            *(exp(-eelmo)-exp(-eel))
            +(dosl-doslmo)*((eelmo+1e0)*exp(-eelmo)
            -(eel +1e0)*exp(-eel)))
            /(eel-eelmo);
        dosl  =CALDOSSUM(eel  ,PHOLE);
        doslmo=CALDOSSUM(eelmo,PHOLE);
        ep0=ep0+((doslmo*eel-dosl*eelmo)
            *(exp(-eelmo)-exp(-eel))
            +(dosl-doslmo)*((eelmo+1e0)*exp(-eelmo)
            -(eel +1e0)*exp(-eel  )))
            /(eel-eelmo);

        //    en0=en0+0.5e0
        //>   *(CALDOSSUM(eel,PELEC)+CALDOSSUM(eelmo,PELEC))
        //>   *(eel-eelmo)*exp(-0.5e0*(eel+eelmo))
        //    ep0=ep0+0.5e0
        //>   *(CALDOSSUM(eel,PHOLE)+CALDOSSUM(eelmo,PHOLE))
        //>   *(eel-eelmo)*exp(-0.5e0*(eel+eelmo))
        eelmo=eel;
    }                            //ENDDO
    en0=2e0*en0;
    ep0=2e0*ep0;
    Ni=sqrt(en0*ep0*exp(-sieg));

    //____end of GETINDENS
    return;
}


void Band::GETDOS( bool calcdos)
{
    //get the density of states

    //LOGICAL calcdos
    int  itab,iptype,ip,iband;
    double dum, sum, maxaov;
    ifstream ftp;

    //calculate the DOS table
    if(calcdos)
    {
        ZUDI();
        //ENDIF
    }

    //read the DOS table from disk
    string filename;
    if(sifl)
      {
        filename = pathname + "/zd.si.asc";
        ftp.open(filename.c_str());
        assert(ftp);
      }
    else if(gefl)
      {
        filename = pathname +  "/zd.ge.trs";
        ftp.open(filename.c_str());
        assert(ftp);
      }
    else if(gaasfl)
      {
        filename = pathname + "/zd.gaas.asc";
        ftp.open(filename.c_str());
        assert(ftp);
      }
    
    //READ (lutmp) itab
        ftp>>itab;               //data is MTAB+1
        if(itab!=MTAB)
        {

            cout<<"error ETABF: wrong number of table entries in file";exit(0);
        }
        //read dos
        for(itab=0;itab<=MTAB;itab++)
        {
            //DO itab=0,MTAB
            for(iband=0;iband<MNB;iband++)
            {
                //DO iband=1, MNB
                dos[iband][itab]=0.0;
                //ENDDO
            }
            //READ (lutmp) dum, (dos[iband][itab],iband=1, nbt)
            ftp>>dum;
            for(iband=0;iband<nbt;iband++)
			{
				ftp>>dos[iband][itab];
				//cout<<iband<<endl<<itab<<endl<<dos[iband][itab]<<endl ;
			}
            //calculate the sum of dos
            for(iptype=0;iptype<NumParType;iptype++)
            {
                //DO iptype=1, NumParType
                sumdos[itab][iptype]=0.0;
                for(iband=bandof[iptype];iband<nband[iptype]+bandof[iptype];iband++)
                {
                    //DO iband=1+bandof[iptype],nband[iptype]+bandof[iptype]
                    sumdos[itab][iptype]=sumdos[itab][iptype]+dos[iband][itab];
                    //ENDDO
                }
                //ENDDO
            }
            //ENDDO
        }
        //read maximum area
        //READ (lutmp) itab
        ftp>>itab;
        for(itab=0;itab<=MTAB;itab++)
        {
            //DO itab=0,MTAB
            for(iband=0;iband<MNB;iband++)
            {
                //DO iband=1, MNB
                dostetmax[iband][itab]=0.0;
                //ENDDO
            }
            //READ (lutmp) dum, (dostetmax[iband][itab], iband=1, nbt)
            ftp>>dum;
            for(iband=0;iband<nbt;iband++)ftp>>dostetmax[iband][itab];
            //ENDDO
        }
        //CLOSE(lutmp)
        ftp.close();

                               

    //normalize density of states and maximum areas of intersection
    for(iptype=0;iptype<NumParType;iptype++)
    {
        //DO iptype=1, NumParType
        DOSMAX[iptype]=0.0;
        //ENDDO
    }

    for(itab=0;itab<=MTAB;itab++)
    {
        //DO itab=0,MTAB
        for(int iband=0;iband<nbt;iband++)
        {
            //DO iband=1,nbt
            dos[iband][itab]=dos[iband][itab]*eV0*pow(spr0,3.0);
            dostetmax[iband][itab]=dostetmax[iband][itab]*velo0/(spk0*spk0);
            //ENDDO
        }
        for(iptype=0;iptype<NumParType;iptype++)
        {
            //DO iptype=1, NumParType
            sumdos[itab][iptype]=sumdos[itab][iptype]*eV0*pow(spr0,3.0);
            DOSMAX[iptype]=Max(DOSMAX[iptype],sumdos[itab][iptype]);
            //ENDDO
        }
        //ENDDO
    }

    //calculate maximum DOS in listfp for Fischetti phonon scattering
    for(iband=0;iband<nbt;iband++)
    {
        //DO iband=1, nbt
        for(itab=0;itab<MWLEFP;itab++)
        {
            //DO itab=1, MWLEFP
            sum=0.0;
            maxaov=0.0;
            for(ip=ptlistfpab[itab][iband];ip<=ptlistfpab[itab][iband]+ntlistfpab[itab][iband]-1;ip++)
            {
                //DO ip=ptlistfpab[itab][iband], ptlistfpab[itab][iband] + ntlistfpab[itab][iband] - 1
                sum   =sum+maxaovtet[tlistfpab[ip]];
                maxaov=Max(maxaov,maxaovtet[tlistfpab[ip]]);
                if(ibt[tlistfpab[ip]]!=iband)
                {

                    cout<<"error GETDOS: Wrong band";exit(0);
                    //ENDIF
                }
                //ENDDO
            }
            dosfpab[iband][itab]=48.0/pow((2.0*PI),3.0)*sum;
            maxaovfpab[iband][itab]=maxaov;
            sum=0.0;
            maxaov=0.0;
            for(ip=ptlistfpem[itab][iband];ip<=ptlistfpem[itab][iband]+ntlistfpem[itab][iband]-1;ip++)
            {
                //DO ip=ptlistfpem[itab][iband], ptlistfpem[itab][iband] + ntlistfpem[itab][iband] - 1
                sum   =sum +      maxaovtet[tlistfpem[ip]];
                maxaov=Max(maxaov,maxaovtet[tlistfpem[ip]]);
                if(ibt[tlistfpem[ip]]!=iband)
                {

                    cout<<"error GETDOS: Wrong band";exit(0);
                    //ENDIF
                }
                //ENDDO
            }
            dosfpem[iband][itab]=48.0/pow((2.0*PI),3.0)*sum;
            maxaovfpem[iband][itab]=maxaov;
            //ENDDO
        }
        //ENDDO
    }

    //End of GETDOS
    return;
    //END
}


//subroutine
//OVERLAP(it,nct,ict)
//corrected by dugang DEC 11;
void Band::BUILDLISTS(void)
{

    bool chsfl;
    int itab,ibeg,iend,ipoint,it;
    int ipoins,ipoinz,ipoinss,ipointfpab,ipointfpem;
    int nct,ict[MCLE],ipoinc;

    int iband,idir;
    //list of tetraheder (table spacing)
    dlist=eV0/0.040;
    dlists=eV0/0.0020;
    /*
          Build list of tetraheda within a certain energy range
          ntlist : Number of tetraheda within the energy range
          ptlist : pointer to the first tetrahedron in a list
           tlist : list of tetraheda
           slist : the same list containing only tetraheda from
          the surface of the BZ
           zlist : list containing only tetrahedra from
          the surface Kz = 0 (used for injection)
           clist : list containing tetrahedra which overlap
          with a given cube
    */
                                 //minim=0;maxim=MNB-1
    for(iband=0;iband<MNB;iband++)
    {
                                 //minm=0;maxim=MWLE-1
        for(itab=0;itab<MWLE;itab++)
        {
            ntlist[itab][iband] = 0;
            nslist[itab][iband] = 0;
        }
        nzlist[iband]=0;
        itab=1;
                                 //minm=0;maxim=MWLES-1
        for(itab=0;itab<MWLES;itab++)
        {
            ntlists[itab][iband] = 0;
        }
    }
    for(itab=0;itab<MCLE;itab++) //minim=0 maxim=MCLE-1
        nclist[itab]=0;
    for(it=0;it<nt;it++)         //minim=0;maxim=nt-1
    {
        //DO it= 1, nt
        //Count tetraheda per energy range
        //get minimum and maximum energy index within the tetrahedron
                                 // + 1;//minim=0;
        ibeg=int(eek[tet[0][it]]*dlist);
                                 // + 1;//minim=0;
        iend=int(eek[tet[3][it]]*dlist);
        if(iend>=MWLE)
        {
            /* 
            WRITE (LUTTYO,*) 'ETABF: iend >=MWLE', ' iend = ', iend
            WRITE (LUOUT ,*) '', ' iend = ', iend
            STOP
            */
            //ENDIF
            cout<<"error ETABF: iend >=MWLE";exit(0);
        }
        if(ibeg>=0)
        {
            for(itab=ibeg;itab<=iend;itab++)
            {
                //DO itab = ibeg,iend
                //add tetrahedron to each list within the energy interval
                ntlist[itab][ibt[it]] = ntlist[itab][ibt[it]] + 1;
                chsfl=false;
                //(-4) and (-5) are the surfaces of the wedge which belong to the
                //surface of the BZ.
                //maybe (-5) and (-6)
                for(idir=0;idir<4;idir++)
                {
                    if((NeibrTet[idir][it]==-5)||(NeibrTet[idir][it]==-6))
                    {
                        chsfl =true;
                    }
                }
                if(chsfl)
                {
                    nslist[itab][ibt[it]] = nslist[itab][ibt[it]] + 1;
                }
                //ENDDO
            }
            //ENDIF
        }
        chsfl =false;
        //(-1) is the surface with Kz=0
        //maybe -2
        for(idir=0;idir<4;idir++)
        {
            //DO idir = 1, 4
            if(NeibrTet[idir][it]==(-2))
            {
                chsfl =true;
            }
        }
        if(chsfl)
        {
            nzlist[ibt[it]] = nzlist[ibt[it]] + 1;
        }
        //list for low energies with finer spacing than tlist
                                 // ) + 1;
        ibeg =int(eek[tet[0][it]]*dlists);
                                 // ) + 1;
        iend =int(eek[tet[3][it]]*dlists);
        if(iend>=MWLES)
        {
            iend=MWLES-1;        //maxim=
        }
        if(ibeg>=0)
        {
            for(itab=ibeg;itab<=iend;itab++)
            {
                //DO itab = ibeg,iend
                ntlists[itab][ibt[it]] = ntlists[itab][ibt[it]] + 1;
            }
            //ENDIF
        }

        //list of tetrahedra in regular cubes (to find states in the X-Minimum)
                                 ////if(ibt[it]==bandof[PELEC]+1)
        if(ibt[it]==bandof[PELEC])
        {
            OVERLAP(it,nct,ict);
            for(itab=0;itab<nct;itab++)
            {
                //DO itab = 1, nct
                nclist[ict[itab]]+=1;
            }
        }
    }

    ipoint = 0;
    ipoins = 0;
    ipoinz = 0;
    ipoinss = 0;
    ipoinc = 0;

    //set pointer to point to the element behind the last element of the list

    for(iband=0;iband<nbt;iband++)
    {
        //DO iband=1,nbt
        for(itab=0;itab<MWLE;itab++)
        {
            //DO itab=1,MWLE
            ipoint=ipoint + ntlist[itab][iband];
                                 // + 1;
            ptlist[itab][iband]=ipoint;
            ipoins=ipoins+nslist[itab][iband];
                                 // + 1;
            pslist[itab][iband]=ipoins;
            if(ipoint>=MWLI)
            {

                cout<<"error ETABF: ipoint > MWLI";exit(0);
                //ENDIF
            }
            if(ipoins>=MSLI)
            {

                cout<<"error ETABF: ipoins > MSLI";exit(0);
            }
        }
        ipoinz = ipoinz + nzlist[iband];
        pzlist[iband] = ipoinz;  // + 1;
        if(ipoinz>=MZLI)
        {

            cout<<"error ETABF: ipoinz > MZLI";exit(0);
        }
        for(itab=0;itab<MWLES;itab++)
        {
            //DO itab=1,MWLES
            ipoinss = ipoinss + ntlists[itab][iband];
                                 // + 1;
            ptlists[itab][iband] = ipoinss;
            if(ipoinss>=MWLIS)
            {

                cout<<"error ETABF: ipoinss > MWLIS";exit(0);
            }
        }
    }
    for(itab=0;itab<MCLE;itab++)
    {
        //DO itab=1,MCLE
        ipoinc = ipoinc + nclist[itab];
        pclist[itab] = ipoinc;   // + 1;
        if(ipoinc>=MCLI)
        {

            cout<<"error ETABF: ipoinc > MCLI";exit(0);
        }
    }

    for(it=0;it<nt;it++)
    {
        //DO it= 1, nt
        //make list
                                 // + 1;minim=0;
        ibeg=int(eek[tet[0][it]]*dlist);
                                 // + 1;minim=0;
        iend=int(eek[tet[3][it]]*dlist);
        //reduce pointer by one and assign the elements
        //after this loop the pointer point to the first element of each list
        if(ibeg>=0)
        {
            for(itab=ibeg;itab<=iend;itab++)
            {
                //DO itab = ibeg,iend
                ipoint = ptlist[itab][ibt[it]] - 1;
                ptlist[itab][ibt[it]] = ipoint;
                tlist[ipoint] = it;
                chsfl =false;
                for(idir=0;idir<4;idir++)
                {
                    //DO idir = 1, 4
                    //if(NeibrTet[idir][it]==(-4)||NeibrTet[idir][it]==(-5))
                                 //be careful!!!!!!!!
                    if(NeibrTet[idir][it]==(-5)||NeibrTet[idir][it]==(-6))
                    {
                        chsfl =true;
                    }
                }
                if(chsfl)
                {
                    ipoins = pslist[itab][ibt[it]] - 1;
                    pslist[itab][ibt[it]] = ipoins;
                    slist[ipoins]=it;
                    //ENDIF
                }
            }
            //ENDIF
        }
        chsfl =false;
        for(idir=0;idir<4;idir++)
        {
            //DO idir = 1, 4
            //if(NeibrTet[idir][it]==-1)
                                 //NeibrTet[][] be careful!!!!!!!!!!!!!
            if(NeibrTet[idir][it]==-2)
            {
                chsfl =true;
            }
        }
        if(chsfl)
        {
            ipoinz = pzlist[ibt[it]] - 1;
            pzlist[ibt[it]] = ipoinz;
            zlist[ipoinz] = it;
        }
                                 // ) + 1;
        ibeg=int(eek[tet[0][it]]*dlists);
                                 // ) + 1;
        iend=int(eek[tet[3][it]]*dlists);
        if(iend>=MWLES)
        {
            iend=MWLES-1;        //maxim=MWLES-1;
        }
        if(ibeg>=0)
        {
            for(itab=ibeg;itab<=iend;itab++)
            {
                //DO itab = ibeg,iend
                ipoinss = ptlists[itab][ibt[it]] - 1;
                ptlists[itab][ibt[it]] = ipoinss;
                tlists[ipoinss] = it;
                //ENDDO
            }
            //ENDIF
        }
                                 /////////////be careful 3.29
        if(ibt[it]==bandof[PELEC])
        {
            OVERLAP(it,nct,ict);
            for(itab=0;itab<nct;itab++)
            {
                //DO itab = 1, nct
                ipoinc = pclist[ict[itab]] - 1;
                pclist[ict[itab]] = ipoinc;
                clist[ipoinc] = it;
            }
            //ENDIF
        }
    }
    //list of tetraheder
    dlistfp = 1.0 / efapeet;
    /*
    C_____Build list of tetraheda within a certain energy range
    C     listfpab : This list contains all tetrahedrons which are in an
    C     energy interval of the width two times the transversal phonon
    C     energy. The spacing of the list is the energy of the transversal phonon.
    C     In the case of absorbtion all tetrahedrons are included which contain
    C     the inital energy plus the phonon energy. In the case of emisson
    C     it is the inital energy minus the phonon energy.
    */
    ///////////////////////15:53
    for(iband=0;iband<MNB;iband++)
    {
        //DO iband=1,MNB
        for(itab=0;itab<MWLEFP;itab++)
        {
            //DO itab=1,MWLEFP
            ntlistfpab[itab][iband]=0;
            ntlistfpem[itab][iband]=0;
        }
    }
    //Count tetraheda per energy range
    for(it=0;it<nt;it++)
    {
        //DO it= 1, nt
        //get minimum and maximum energy index within the tetrahedron
                                 // ) + 1;
        ibeg =int(eek[tet[0][it]]*dlistfp);
                                 // ) + 1;
        iend =int(eek[tet[3][it]]*dlistfp);
        if(iend>=MWLEFP)
        {

            cout<<"error ETABF: iend >=MWLEFP";exit(0);
        }
        if(ibeg>=0)
        {
            for(itab=ibeg;itab<=iend;itab++)
            {
                //DO itab = ibeg,iend
                //add tetrahedron to each list within the energy interval
                ntlistfpab[itab][ibt[it]] = ntlistfpab[itab][ibt[it]] + 1;
                ntlistfpem[itab][ibt[it]] = ntlistfpem[itab][ibt[it]] + 1;
            }
            if(ibeg-1>=0)        //be careful!!!!!!!!!!!!!!!!!!
            {
                //add tetrahedron to each list within the energy interval
                ntlistfpab[ibeg-1][ibt[it]] = ntlistfpab[ibeg-1][ibt[it]]+1;
                //ENDIF
            }
            if(iend+1<MWLEFP)    //be careful!!!!!!!!!!!
            {
                //add tetrahedron to each list within the energy interval
                ntlistfpem[iend+1][ibt[it]] = ntlistfpem[iend+1][ibt[it]]+1;
                //ENDIF
            }
            //ENDIF
        }
    }
    //set pointers
    ipointfpab = 0;
    ipointfpem = 0;

    //set pointer to point to the element behind the last element of the list
    for(iband=0;iband<nbt;iband++)
    {
        //DO iband=1,nbt

        for(itab=0;itab<MWLEFP;itab++)
        {
            //DO itab=1,MWLEFP
            ipointfpab = ipointfpab + ntlistfpab[itab][iband];
                                 // + 1;
            ptlistfpab[itab][iband]=ipointfpab;
            if(ipointfpab>=MWLIFP)
            {

                cout<<"error ETABF: ipointfpab > MWLIFP";exit(0);
                //ENDIF
            }
            ipointfpem = ipointfpem + ntlistfpem[itab][iband];
                                 // + 1;
            ptlistfpem[itab][iband] = ipointfpem;
            if(ipointfpem>=MWLIFP)
            {

                cout<<"error ETABF: ipointfpem > MWLIFP";exit(0);
            }
        }
    }

    //make list
    for(it=0;it<nt;it++)
    {
        //DO it= 1, nt
                                 // ) + 1;
        ibeg=int(eek[tet[0][it]]*dlistfp);
                                 // ) + 1;
        iend=int(eek[tet[3][it]]*dlistfp);
        //reduce pointer by one and assign the elements
        //after this loop the pointer point to the first element of each list
        if(ibeg>=0)
        {
            for(itab=ibeg;itab<=iend;itab++)
            {
                //DO itab = ibeg,iend
                ipointfpab = ptlistfpab[itab][ibt[it]] - 1;
                ptlistfpab[itab][ibt[it]] = ipointfpab;
                tlistfpab[ipointfpab] = it;
                ipointfpem = ptlistfpem[itab][ibt[it]] - 1;
                ptlistfpem[itab][ibt[it]] = ipointfpem;
                tlistfpem[ipointfpem] = it;
            }
            if(ibeg-1>=0)        //be careful !!!!!!!!!!!
            {
                ipointfpab = ptlistfpab[ibeg-1][ibt[it]] - 1;
                ptlistfpab[ibeg-1][ibt[it]] = ipointfpab;
                tlistfpab[ipointfpab] = it;
                //ENDIF
            }
            if(iend+1<MWLEFP)    //be careful !!!!!!!!!!
            {
                ipointfpem = ptlistfpem[iend+1][ibt[it]] - 1;
                ptlistfpem[iend+1][ibt[it]] = ipointfpem;
                tlistfpem[ipointfpem] = it;
                //ENDIF
            }
            //ENDIF
        }
    }

    //End of BUILDLISTS
    return;
}


void Band::IELEC(string path)
{
    //     process elec.scatt command
    //____local variables
    bool calcdos,calcband;
    int i,itab,iptype;
    double sck;
    // char cscat[CPCVL];
    
    //____set parameters for bandstructure and particle type
    //Typename[PELEC]="electron";
    //Typename[PHOLE]="    hole";
    //Typename[POXEL]="oxidelec";

    pathname = path;
    
    typemat[PELEC]=SILICON;
    typemat[PHOLE]=SILICON;
    typemat[POXEL]=OXIDE;

    //____set values for energy discretization [0:MTAB] -> [emin,emax]
    dtable=0.0010/eV0;
    emin=0.0;
    emax=(double)(MTAB)*dtable;

    
    for(itab=0;itab<=MTAB;itab++)
    {                            //DO itab=0,MTAB
        energy[itab]=(double)(itab)*dtable+emin;
    }                            //ENDDO

    //____read card parameters
    opcalc  =false;              //lval(1,linum)
    calcband =false;             //lval(2,linum)
    calcdos  =false;             //lval(3,linum)
    kaneiifl =false;             //lval(4,linum)
    thomiifl =true;              //lval(5,linum)
    fisciifl =false;             //lval(6,linum)
    sanoiifl =false;             //lval(7,linum)
    iifl    =true;               //lval(8,linum)
    seciifl =false;              //true;//lval(9,linum)
    tunxfl  =false;              //lval(10,linum)

    bhfl    =false;              //true;//lval(12,linum)
    frickfl =false;              //true;//lval(13,linum)
    bhmrtfl =false;              //lval(14,linum)
    injoxfl =false;              //lval(15,linum)
    jacophfl=true;               //lval(16,linum)
    fiscphfl=false;              //lval(17,linum)

    //____normalize phonon temperature and deformation potentials
    temptag=140.0/T0;            //dval( 1,linum)/T0
    templag=215.0/T0;            //dval( 2,linum)/T0
    templog=720.0/T0;            //dval( 3,linum)/T0
    temptaf=220.0/T0;            //dval( 4,linum)/T0
    templaf=550.0/T0;            //dval( 5,linum)/T0
    temptof=685.0/T0;            //dval( 6,linum)/T0
    dftag  =4.65e9/dpc0;         //dval( 7,linum)/dpc0
    dflag  =7.44e9/dpc0;         //dval( 8,linum)/dpc0
    dflog  =1.023e11/dpc0;       //dval( 9,linum)/dpc0
    dftaf  =2.79e9/dpc0;         //dval(10,linum)/dpc0
    dflaf  =1.86e10/dpc0;        //dval(11,linum)/dpc0
    dftof  =1.86e10/dpc0;        //dval(12,linum)/dpc0
    dfelast=8.7/eV0;             //dval(13,linum)/eV0

    temphop=680.0/T0;            //dval(14,linum)/T0
    dfhop  =4e10/dpc0;           //dval(15,linum)/dpc0
    dfhelast=9.3/eV0;            // dval(16,linum)/eV0

    tempoeop=735.0/T0;           //dval(17,linum)/T0
    dfoeop  =2.2e10/dpc0;        //dval(18,linum)/dpc0
    dfoeelast=25.0/eV0;          // dval(19,linum)/eV0

    iifacelec   =0.16;           //dval(20,linum)

    difpr[PELEC]=0.16;           //dval(21,linum)
    difpr[PHOLE]=0.35;           //dval(22,linum)
    difpr[POXEL]=0.16;           //dval(23,linum)

    sioxbgo     =3.2/eV0;        //dval(24,linum)/eV0
                                 //dval(25,linum)/eV0*SQRT(field0)
    beta        =2.15e-5/eV0*sqrt(field0);

    iifachole   =1.14e12;        //dval(26,linum)
    hiithresh=1.49/eV0;          //dval(27,linum)/eV0
    hiiexp   =3.4;               //dval(28,linum)
    iifachole  =iifachole*pow(eV0,hiiexp)*time0;

    ephb=1.6;                    //dval(29,linum)

    dmox  =0.5;                  //dval(30,linum)

    mell=0.9116;                 //dval(31,linum)//
    melt=0.1946;                 //dval(32,linum)
    meld=pow((mell*melt*melt),(1.0/3.0));

    efoplow =2.47e10/dpc0;       //dval(33,linum)/dpc0
    efophigh=2.97e10/dpc0;       //dval(34,linum)/dpc0
    efopee  =6.2e-2/eV0;         //dval(35,linum)/eV0

    efaplow =1.7/eV0;            //dval(36,linum)/eV0
    efaphigh=2.4/eV0;            //dval(37,linum)/eV0
    efapeet =4.43e-2/eV0;        //dval(38,linum)/eV0
    efapeel =2.21e-2/eV0;        //dval(39,linum)/eV0

    dmc[PELEC]=0.289;            //dval(40,linum)
    dmc[PHOLE]=0.349;            //dval(41,linum)
    dmc[POXEL]=0.5;              //dval(42,linum)

    //cscat ="scatt.dat";//cval(1,linum)

    //____parameters for surface scattering

    //____only one impact ionization model can be specified
    i=0;
    if (kaneiifl) i=i+1;
    if (thomiifl) i=i+1;
    if (sanoiifl) i=i+1;
    if (fisciifl) i=i+1;
    if (i>1)
    {
    }                            //ENDIF

    //____calculate certain bandstructure values
    if (calcband)
    {                            //THEN
        /*
        NEIB
        WRITE (LUTTYO,*)
        >      'IELEC: calculation of bandstruc-values is  finished'
        WRITE (LUOUT ,*)
        >      'IELEC: calculation of bandstruc-values is  finished'
        */
    }                            //ENDIF

    //____read bandstruc-values from files

    READBS();

    //____build lists for finding states in k-space
    BUILDLISTS();
    //____get density of states
    GETDOS(calcdos);

    //____build lists for phonon scattering
    BUILDPHSCATT();

    //____intrinsic carrier density
    GETINDENS();

    //____find the maximum of scattering rate
    for(iptype=0;iptype<NumParType;iptype++)
    {                            //DO iptype=1,NumParType
        if (nband[iptype]!=0)
        {                        //THEN
            gamma[iptype]=0.0;
            for(itab=0;itab<=MTAB;itab++)
            {                    //DO itab=0,MTAB
                sck=CALSCATTSUM(energy[itab],bandof[iptype]);
                //sck=CALSCATTSUM(energy[itab],1+bandof[iptype])
                gamma[iptype]=Max(gamma[iptype],sck);
            }                    //ENDDO
            /*
            WRITE (LUTTYO,1000) Typename[iptype],gamma[iptype]/time0
            WRITE (LUOUT ,1000) Typename[iptype],gamma[iptype]/time0
            */
        }                        //ENDIF
    }                            //ENDDO

    init_inject_and_density_table() ;
    //____end of IELEC
    return;
}

bool Band::out_of_range(double Ef) {
  if ((Ef > Efmax) || (Ef < Efmin))
    return true;
  return false;
}

double Band::density_to_Ef(double dens) {
  int idx;
  idx = lower_bound(Ef_density.begin(), Ef_density.end(), dens) - Ef_density.begin();
//  cout << Ef_value[idx] << ' ' << Ef_density[idx] * conc0 << endl;
  return Ef_value[idx];
}

double Band::Ef_to_density(double Ef) {
  int idx;

  idx = (int) ((Ef - Efmin) / deltaEf);

  if (Ef > Efmax) idx = NEf - 1;
  if (Ef < Efmin) idx = 0;
  
  return Ef_density[idx];
}

double Band::Ef_to_cross_number(double Ef) {
  int idx;

  idx = (int) ((Ef - Efmin) / deltaEf);

  if (Ef > Efmax) idx = NEf - 1;
  if (Ef < Efmin) idx = 0;
  
  return Ef_cross_number[idx];
}

void Band::init_inject_and_density_table() {
  ifstream ftp;
  string filename;
  int i;

  filename = pathname + "/density_cross_number.txt";

  ftp.open(filename.c_str());
  ftp >> NEf;
  ftp >> Efmin >> Efmax;
  Efmin /= eV0;
  Efmax /= eV0;
  deltaEf = (Efmax - Efmin) / NEf;

  Ef_value.resize(NEf);
  Ef_density.resize(NEf);
  Ef_cross_number.resize(NEf);

  for (i = 0;i < NEf; i ++){
    ftp >> Ef_value[i] >> Ef_density[i] >> Ef_cross_number[i];
    Ef_value[i] /= eV0;
    Ef_density[i] = Ef_density[i] * 1e27 / conc0;
    Ef_cross_number[i] = Ef_cross_number[i] * 1e27 * spr0 * spr0 * time0; 
  }

}

void Band::output_tet(int itet) {
  int i;
  cout << "vertex of tetrahedron " << itet << endl;
  for (i = 0;i < 4; i ++)
    cout << "( " << xkk[tet[i][itet]] << ','
	 << ykk[tet[i][itet]] << ','
	 << zkk[tet[i][itet]] << ')' << endl;

}
