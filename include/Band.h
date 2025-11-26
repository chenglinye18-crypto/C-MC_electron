
#ifndef __CLASS_BAND__
#define __CLASS_BAND__

#include "Constant_GlobalVariable_MathFunction.h"
#include "utils.h"
#include <string>
#include <vector>

class Band
{
    //    contains Name for each particle TYPE
    //CHARACTER*8 Typename[NumParType]
public:
    // 解析(类抛物线)模型开关与建表开关，仅作用于电子
    bool use_analytic_band = false;
    bool build_analytic_tables = true;

    // 解析模型预留表：E-k、DOS、v-k（仅声明，后续填充）
    std::vector<std::string> analytic_ek_table;
    std::vector<std::string> analytic_dos_table;
    std::vector<std::string> analytic_vk_table;

    char Typename[NumParType][9];
    //     number of scattering processes (electrons)
    int scpre;

    //     number of scattering processes (holes)
    int scprh;

    //     number of scattering processes (oxide electrons)
    int scproe;

    //    field with 48 Symetrie-operations
    int matsym[6][48];

    //    inverse of matsym
    // int indmat[-1:1,-1:1,-1:1,3,3,3]
    int indmat[3][3][3][3][3][3];

    //    matsym with a scalar argument
    //int matsymlin[6*48]

    //    number tetraeder in in the energy-fields MWLE
    int ntlist[MWLE][MNB];

    //    pointer to tlist
    int ptlist[MWLE][MNB];
    //     particle type of the energy bands
    int partyp[MNB];

    //    sequentiell list of tedraeder
    int tlist[MWLI];

    //    number of surface tetraeder in in the energy-fields MWLE
    int nslist[MWLE][MNB];

    //    pointer to slist
    int pslist[MWLE][MNB];

    //    sequentiell list of surface tedraeder
    int slist[MSLI];

    //    number of surface tetraeder in the energy-fields MWLE
    int nzlist[MNB];

    //    pointer to slist
    int pzlist[MNB];

    //    sequentiell list of surface tedraeder
    int zlist[MZLI];

    //    number of tetraeder in the energy-fields MWLES
    //    [finer energy spacing than tlist]
    int ntlists[MWLES][MNB];

    //    pointer to tlists
    int ptlists[MWLES][MNB];

    //    sequentiell list of tedraedra
    int tlists[MWLIS];

    //    number of tetraeder in each cube
    int nclist[MCLE];

    //    pointer to clist
    int pclist[MCLE];

    //    sequentiell list of tetrahedra
    int clist[MCLI];

    //    number tetraeder in the energy-fields MWLEFP [Fischetti, absorbtion]
    int ntlistfpab[MWLEFP][MNB];

    //    pointer to tlistfpab [Fischetti, absorbtion]
    int ptlistfpab[MWLEFP][MNB];

    //    sequentiell list of tedraeder [Fischetti, absorbtion]
    int tlistfpab[MWLIFP];

    //    number tetraeder in the energy-fields MWLEFP [Fischetti, emisson]
    int ntlistfpem[MWLEFP][MNB];

    //    pointer to tlistfpab [Fischetti, emisson]
    int ptlistfpem[MWLEFP][MNB];

    //    sequentiell list of tedraeder [Fischetti, emisson]
    int tlistfpem[MWLIFP];

    //    contains total number of bands
    int nbt;

    //    contains number of bands for each particle type
    int nband[NumParType];

    //    OFfset of BAND index
    int bandof[NumParType];

    //    contains material type for each particle TYPE
    int typemat[NumParType];

    //    indices(����) of the for nodes of each tetrahedron
    int tet[4][MNTet];

    //    number of k-space grid points
    int nk;

    //    number of tetraheda
    int nt;

    //    energy band index of each tetrahdron
    int ibt[MNTet];

    //    list of neighbour tetraheda
    int NeibrTet[4][MNTet];

    //    in the case of a change of the wedge newisym gives the new symmetry trans.
    int newisym[5][48];
    //     flag for writing the scattering files to a rgraph file
    bool opcalc;

    //     flag for Brooks-Herring impurity scattering
    bool bhfl;
    //     flag for Sano formula for electron impact ionization rate
    bool sanoiifl;
    //     flag for use of the inverse microscopic relaxation time instead
    //     of scattering rate
    bool bhmrtfl;                ///taf.f scatt.f

    //     flag for Kane formula for electron impact ionization rate
    bool kaneiifl;

    //     flag for Thoma formula for electron impact ionization rate
    bool thomiifl;

    //     flag for Fischetti formula for electron impact ionization rate
    bool fisciifl;

    //     flag for impact ionization scattering
    bool iifl;



    //    self scattering event
    // bool selfscfl;
    //     flag for generation of secondary particles
    bool seciifl;

    //     flag for tunneling of electrons in the oxide (only negative x-direction)
    bool tunxfl;
    //     flag for Jacoboni phonon system
    bool jacophfl;

    //     falg for fischetti phonon system
    bool fiscphfl;

    //     empirical correction factor for BH impurity scattering
    bool frickfl;

    //     flag for particle transition from silicon into oxide
    bool injoxfl;
    //    energy of the base of the tetrahedron
    double ebzp;

    //    groupvelocity within the tetrahedron [irreducible wedge]
    double vgx,vgy,vgz;

    //    maximum of density of states for each particle type
    double DOSMAX[NumParType];

    //    table of density of states for each band
    double dos[MNB][MTAB+1];

    //    table of the sum of DOS over all bands for each particle type
    double sumdos[MTAB+1][NumParType];

    //    DOS of all tetrahedra in the list, each tetrahedron
    //    contributing with its maximum DOS [Fischetti, absorbtion]
    double dosfpab[MNB][MWLEFP];

    //    DOS of all tetrahedra in the list, each tetrahedron [Fischetti, emisson]
    //    contributing with its maximum DOS
    double dosfpem[MNB][MWLEFP];

    //    maximum AOV for list [Fischetti, aborbtion]
    double maxaovfpab[MNB][MWLEFP];

    //    maximum AOV for list [Fischetti, emisson]
    double maxaovfpem[MNB][MWLEFP];

    //    energy for tables [scatt and DOS]
    double energy[MTAB+1];

    //    energy spacing of tables
    double dtable;

    //    minimum energy in table
    double emin;

    //    maximum energy in table
    double emax;

    //    energy spacing for treaheder list
    double dlist;

    //    energy spacing for treaheder list
    double dlists;

    //    energy spacing for treaheder list [Fischetti]
    double dlistfp;

    //    array of grid points in k-space
    double xkk[MNK],ykk[MNK],zkk[MNK];

    //    array of energy of grid points in k-space
    double eek[MNK];

    //    array of comp. of velocity in tetrahedron
    double vt[4*MNTet];

    //    absolute value of the group velocity
    double vgt[MNTet];

    //    center coordinates of each tetrahedron
    double xkct[MNTet],ykct[MNTet],zkct[MNTet];

    //    maximum dos in tet. for each band and energy
    double dostetmax[MNB][MTAB+1];

    //    maximum area over velocity for each tet.
    double maxaovtet[MNTet];
    //    inverse mass [1:xx, 2:xy, 3:xz, 4:yy, 5:yz, 6:zz]
    //    centered in the middle of the tetrahedron
    double massinv[6][MNTet];

    //    array of normal vectors of tetrahedra surfaces
    //    first index = tet. surface
    //    datantlin[0,*,*] = x-comp
    //    datantlin[1,*,*] = y-comp
    //    datantlin[2,*,*] = z-comp
    //    datantlin[3,*,*] = distances of tet. surfaces from origin
    //double datantlin[4][4][MNTet];
    //    same as datant with a scalar argument
    // MNTet: Maximum numeber of Tetrahedron
    double datantlin[4][4][MNTet];
    //    factor for FASTSURF for energy between e1 and e2
    double faclow[MNTet];
    //    factor for FASTSURF for energy between e2 and e3
    double facmedium[MNTet];
    //    factor for FASTSURF for energy between e3 and e4
    double fachigh[MNTet];

    //     prefactor for hole II
    double iifachole;
    //     hole II threshold energy
    double hiithresh;
    //     hole II energy exponent
    double hiiexp;
    //     enhancement factor for e-P scattering except for intraband scattering
    //     in the first conduction band
    double ephb;
    //     oxide electron mass
    double dmox;
    //     expectation value of the trace of the inverse mass tensor for equilibrium
    double dmc[NumParType];
    //     longitudinal mass in the minimum of the first conduction band
    double mell;
    //    transversal mass in the minimum of the first conduction band
    double melt;
    //    DOS mass in the minimum of the first conduction band
    double meld;
    //    temperature of the transversal acoustic g-phonon (electrons)
    double temptag;
    //    temperature of the longitudinal acoustic g-phonon (electrons)
    double templag;
    //    temperature of the longitudinal optical g-phonon (electrons)
    double templog;
    //    temperature of the transversal acoustic f-phonon (electrons)
    double temptaf;
    //    temperature of the longitudinal acoustic f-phonon (electrons)
    double templaf;
    //    temperature of the transversal optical f-phonon (electrons)
    double temptof;
    //    temperature of the optical phonon (holes)
    double temphop;
    //    temperature of the optical phonon (oxide electrons)
    double tempoeop;
    //    deformation potential constant of the transversal acoustic g-phonon
    //    (electrons)
    double dftag;
    //    deformation potential constant of the longitudinal acoustic g-phonon
    //    (electrons)
    double dflag;
    //    deformation potential constant of the longitudinal optical g-phonon
    //    (electrons)
    double dflog;
    //    deformation potential constant of the transversal acoustic f-phonon
    //    (electrons)
    double dftaf;
    //    deformation potential constant of the longitudinal acoustic f-phonon
    //    (electrons)
    double dflaf;
    //    deformation potential constant of the transversal optical f-phonon
    //    (electrons)
    double dftof;
    //    deformation potential for elastic acoustic phonon scattering (electrons)
    double dfelast;
    //    deformation potential constant for optical phonon scattering (holes)
    double dfhop;
    //    deformation potential for elastic acoustic phonon scattering (holes)
    double dfhelast;
    //    deformation potential constant for optical phonon scattering
    //    (oxide electrons)
    double dfoeop;
    //    deformation potential for elastic acoustic phonon scattering
    //    (oxide electrons)
    double dfoeelast;
    //    prefactor for electron impact ionization
    double iifacelec;

    //    coupling constant for optical phonons low band
    double efoplow;

    //    coupling constant for optical phonons high band
    double efophigh;

    //    energy for optical phonons
    double efopee;

    //    coupling constant for optical phonons low band
    double efaplow;

    //    coupling constant for optical phonons high band
    double efaphigh;

    //    energy for transversal acoustic phonons
    double efapeet;

    //    energy for longitudinal acoustic phonons
    double efapeel;
    //     intrinsic carrier concentration
    double Ni;
    //    lattice vectors for change of BZ
    double xbz[5][48],ybz[5][48],zbz[5][48];
    //     maximaum scattering rate for particle scattering
    double gamma[NumParType];
    //     maximum scattering rate in a tetrahedron
    double gamtet[MNTet];

    //     fraction of diffusive scattering at Si/SiO2 interface
    double difpr[NumParType];


    //     
    /**
     * @brief phonon scattering rate prefactor table for electrons
     * 即散射率公式中除 态密度 之外的部分
     * 
     * @details
     * MNScaEle: M_N_Sca_Ele: Maximum Number of scattering processes for electrons
     * NBE: maximum number of conduction bands
     */
    double scatte[MNScaEle][NBE][NBE];

    //     scattering rate table for electron impact ionization
    double scattiie[NBE][NBE][MTAB+1];

    //     phonon scattering rate prefactor table for holes
    double scatth[MNScaHole][NBH][NBH];

    //     scattering rate table for hole impact ionization
    double scattiih[NBH][NBH][MTAB+1];

    //     phonon scattering rate prefactor table for oxide electrons
    double scattoe[MNScaOxEle][NBOE][NBOE];

    //     DOS for the different phonon scattering processes final energies
    //     (electrons)
    double dose[MNScaEle][NBE][MTAB+1];

    //     DOS for the different phonon scattering processes final energies (holes)
    double dosh[MNScaHole][NBH][MTAB+1];

    //     DOS for the different phonon scattering processes final energies
    //     (oxide electrons)
    double dosoe[MNScaOxEle][NBOE][MTAB+1];

    //     sum of scattering rates for each band (II + Phonon)
    double sumscatt[MTAB+1][MNB];
    
    vector<double> Ef_value, Ef_cross_number, Ef_density;

    int NEf;
    double deltaEf, Efmin, Efmax;

    string pathname;
    
    void BUILDPHSCATT(void);
    double CALDOSTETMAX(double ee,int iband);
    void CALSCATTE(double ee,double*sca,int iband);
    void CALSCATTH(double ee,double*sca, int iband);
    void CALSCATTOE(double ee,double*sca,int iband);
    double CALSCATTSUM(double ee,int iband);
    void CALCFAC(void);
    double CALDOS(double eed,int ib);
    double CALDOSSUM(double eed,int ip);
    void CONFBZ(double &xkl,double &ykl,double &zkl);
    double EIIRATE(double ee);
    double FASTSURF(double ee,int itet);
    double FII(double x1,double y1,double z1,
        double x2,double y2,double z2,
        double x3,double y3,double z3);
    void GETMAXAOVTET(void);
    void GETINDENS(void);
    double HIIRATE(double ee);
    void HESSTET(void);
    void SETNEIBSYM(void);
    void SETMATSYM(void);
    double SURF(double eps,int it);
    void VELOTET(void);
    void ZUDI(void);
    void OutWedge(double &xout,double &yout,double &zout,
        double xin,double yin,double zin,int isym);
    void OVERLAP(int it,int&nct,int*ict);
    void  READBS(void);
    void GETDOS(bool calcdos);
    void BUILDLISTS(void);
    void init_inject_and_density_table() ;
    double Ef_to_cross_number(double Ef) ;
    double Ef_to_density(double Ef) ;
    double density_to_Ef(double dens) ;
    bool out_of_range(double Ef) ;
    void output_tet(int);
    // 解析(类抛物线)模型建表接口（仅电子）：E-k、DOS、v-k
    void BuildAnalyticTables(const string& path);
    void BuildAnalyticEKTable(const string& path);
    void BuildAnalyticDOSTable(const string& path);
    void BuildAnalyticVKTable(const string& path);
public :
    //    silicon to silicon oxide conduction band discontinuity
  double sioxbgo;
  //     image potential factor for Schottky barrier lowering
  double beta;
  
  void IELEC(string);
  
};

#endif  
