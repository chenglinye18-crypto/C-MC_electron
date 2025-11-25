#ifndef __CONST_GLOBAL_MATH_FUNC__
#define __CONST_GLOBAL_MATH_FUNC__
#include <math.h>
#include<assert.h>
#include <time.h>
#include<string.h>
#include <iostream>
#include<fstream>
#include <iomanip>
using namespace std;
//#include "petscksp.h"
#include<string.h>

//mathematical constants
const double THIRD =3.333333333333333333333333e-1;
const double TWOTHI=6.666666666666666666666667e-1;
const double PI    =3.141592653589793238462643e0;
const double TWOPI =2.0*PI;
const double SQRTPI=1.772453850905516027298167e0;
const double LNPI  =1.144729885849400174143427e0;
const double SQRT2 =1.4142135623730950488e0;
const double SQRT2I=1.0/SQRT2;
const double SQRT3 =1.7320508075688772935e0;


//silicon bulk parameters:
//	sia0=5.43d-10  ! m        ! lattice parameter
//	sirho=2.33d3   ! kg/m**3  ! crystal density
//	siul=9.0d3     ! m/s      ! longitudinal sound velocity
//	siut=5.3d3     ! m/s      ! transverse sound velocity
//	a0pi           ! 1/m      ! 2*PI/sia0
//	sieg           ! eV       ! silicon band gap
extern double sia0,sirho,siul,siut,a0pi,sieg;


//Common for normalization
//	T0     = temperature . . . . . . . . . [K]
//	eV0    = energy. . . . . . . . . . . . [eV]
//	em0    = electon rest mass . . . . . . [kg]
//	hq0    = Planck constant . . . . . . . [eVs]
//	ec0    = electron charge . . . . . . . [As]
//	rmom0  = momentum. . . . . . . . . . . [eVs/m]
//	spr0   = r-space . . . . . . . . . . . [m]
//	spk0   = k-space . . . . . . . . . . . [1/m]
//	time0  = time. . . . . . . . . . . . . [s]
//	velo0  = velocity. . . . . . . . . . . [m/s]
//	cvr    = speed of light / velo0
//	pot0   = el. potential . . . . . . . . [V]
//	field0 = el. field . . . . . . . . . . [V/m]
//	conc0  = concentration . . . . . . . . [1/m**3]
//  dens   = mass density. . . . . . . . . [kg/m**3]
//  dpc0   = optical deformationpotential. [eV/m]
//  scrt0  = scattering rate . . . . . . . [1/s]
//  curr0  = terminal current (2D) . . . . [A/m]

extern double T0,eV0,em0,hq0,ec0,rmom0,spr0,spk0,time0,velo0;
extern double cvr,pot0,field0,conc0,dens0,dpc0,scrt0,curr0, Nc, Ncc, psi_si, phi_top, N_cur;
extern double mx[3], my[3], mz[3], Nccc[3];
extern double ml, mt;

//physical constants:
//	ec     = electron charge         [ A*s ]
//  em     = electron rest mass      [ kg ]
//  planck = Planck s constant       [ eV*s ]
//  clight = speed of ligth          [ m/s ]
//  boltz  = Boltzmann sconstant     [ eV/K ]
//  eps0   = dielectric constant     [ (A*s)/(V*m) ]
//  fsc    = fine structure constant
const double EC     = 1.602176e-19;//updated
const double EM     = 9.109383e-31;//updated
const double PLANCK = 6.582119e-16;//	g*j`gv//updated
const double CLIGHT = 299792458;//updated
const double BOLTZ  = 8.617343e-5;//updated
const double EPS0   = 8.854187817e-12;//updated
const double FSC    = EC/(PLANCK*CLIGHT*4.0*PI*EPS0);




//------------------------Kzz!jWSpe------------------------
//,geg`bTW[9eN9ev *YA~tN
//maximum number of conduction bands
const int NBE=4;
//maximum number of valence bands
const int NBH=3;
//maximum number of conduction bands in the oxide
const int NBOE=1;
//maximum number energy bands
const int MNB=NBE+NBH+NBOE;


//maximum number of tetraeder in the onedimensional list
const int MWLI=470000;
//maximum number of energy-fields for tetrdaederlist
const int MWLE=325;

//maximum number of tetraeder in the onedimensional list
//of tetrahedrons at the surface of the wedge
const int MSLI=16000;

//number of cubes in x-direction
const int MNCubeX=100;
//number of cubes in y-direction
const int MNCubeY=100;
//number of cubes in z-direction
const int MNCubeZ=100;

//the same as above for Fischetti type phonon scattering
const int MWLIFP=480000;
//the same as above for Fischetti type phonon scattering
const int MWLEFP=325;

//the same as above for a finer list
const int MWLIS=101000;
//the same as above for a finer list
const int MWLES=60;

//the same as above for the cube oriented list
const int MCLI=50000;
//the same as above for the cube oriented list
const int MCLE=MNCubeX*MNCubeY*MNCubeZ;

//Maximum number of k-space points
const int MNK=13000;

//Maximum number of tetraheda
const int MNTet=56000;

//------------------------[zz!jWSpe------------------------
//Maximum Number of Cuts in X-direction
const int MNCutX=100;
//Maximum Number of Cuts in Y-direction
const int MNCutY=100;
//Maximum Number of Cuts in Z-direction
const int MNCutZ=100;

//Maximum Number of Cells in X-direction
const int MNCellX=MNCutX-1;
//Maximum Number of Cells in Y-direction
const int MNCellY=MNCutY-1;
//Maximum Number of Cells in Z-direction
const int MNCellZ=MNCutZ-1;
//Maximum Number of Cells
const int MNCell=MNCellX*MNCellY*MNCellZ;

//Maximum Number of Points in X-direction
const int MNPointX=MNCutX;
//Maximum Number of Points in Y-direction
const int MNPointY=MNCutY;
//Maximum Number of Points in Z-direction
const int MNPointZ=MNCutZ;
//Maximum Number of Point
const int MNPoint=MNPointX*MNPointY*MNPointZ;

//Maximum Number of Regions
const int MNRegion=16;
//Maximum Number of Material Type
const int MNMaterialType=6;
//Maximum Number of Surface//Surface refers to Si/SiO2 interface
const int MNSurface=30;
//Token for a material that represents nothing (Neumann boundary conditions
//on the surface for every thing
const int VOID=0;
//Token for vacuum
const int VACUUM=1;
//Token for oxide
const int OXIDE=2;
//Token for silicon
const int SILICON=3;
//Token for surface at Si02-Si
const int SILINE=4;
//Token for surface at Polygate-SiO2
const int GATELINE=5;
//Token for IGZO
const int IGZO=6;

//Maximum Number of Contacts
const int MNContact=16;
//Maximum Number of Contact Planes
const int MNContactPlane=10;
//Maximum Number of Contact Vapp Steps
const int MNContactVappStep=100;
//Token for no cont
const int NOCONT=0;
//Token for silicon contact
const int SICONT=1;
//Token for a gate contact
const int GATE=2;
//Token for a catch contact
const int CCONT=3;

//Token for direction
const int UP=0,DOWN=1,LEFT=2,RIGHT=3,FRONT=4,BACK=5;
//Token for corner
const int ULF=0,DLF=1,URF=2,DRF=3,ULB=4,DLB=5,URB=6,DRB=7;

//Maximum Number of interface charges
const int MNQsurf=16;

//Token for motion rule pass (particle passes from one quadrant to another)
const int PASS=1;
//Token for motion rule reflect (particle is reflected at the boundary)
const int REFLECT=2;
//Token for motion rule scattox (particle is reflected or scattered at the SiO2 boundary)
const int SCATTOX=3;
//Token for motion rule period (particle is moved to the opposite boundary of the quadrant)
const int PERIOD=4;
//Token for motion rule generate (new particle is generated at the other
//side of the boundary, original particle experiences period boundary condition)
const int GENERATE=5;
const int GENREF=8;
//Token for motion rule catch (particle is annihilated)
const int CATCH=6;
//Token for motion rule catchgate (particle is annihilated by a gate)
const int CATCHGATE=7;


//------------------------|P[Spe------------------------
//maximum number of particles
const int MNParticle=2000000;


//Parameters for particle type
const int NumParType=3,PELEC=0,PHOLE=1,POXEL=2;


//5uw&{S
extern double ChargeSign[NumParType];

//	gHe(
extern double MASSeff[NumParType];


//Maximum number of scattering processes (electrons)
const int MNScaEle=14;
//Maximum number of scattering processes (holes)
const int MNScaHole=4;
//Maximum number of scattering processes (oxide electrons)
const int MNScaOxEle=3;




//------------------------!jb0pef[{vsQ------------------------
//maximum number of band width points
const int MBW=MNPointX*MNPointZ;

//Maximum number of time slice
const int MNDt=1000000;



//------------------------Multiple Refresh------------------------
//Maximum Number of Energy Ranges
const int MNERange=70;
//Maximu Number of Phase Space Cells
const int MNPSC= MNCell * 10;
//Replace-N_v4Nepe~v'Y
const int MNRP=100000;

//maximum number of tetraeder in the onedimensional list
//of tetrahedrons at the kz=0 surface of the wedge
const int MZLI=4500;
//maxnumber of energy dependent dos-table and scattering-table
const int MTAB=7500;

//load GaAs band strucutre
extern bool gaasfl;
//load silicon band structure
extern bool sifl;
extern bool gefl;


extern double Max(double p1,double p2);

extern  double Max(double p1,double p2,double p3);

extern  double Max(double p1,double p2,double p3,double p4);

extern  double Max(double p1,double p2,double p3,double p4,double p5);

extern  double Min(double p1,double p2);

extern  double Min(double p1,double p2,double p3);

extern  double Min(double p1,double p2,double p3,double p4);

extern  double Min(double p1,double p2,double p3,double p4,double p5,double p6, int & flag);

#endif
