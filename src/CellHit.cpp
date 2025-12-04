
#include "mcmodel.h"
#include "Band.h"

double MeshQuantities::CellTime()
{
    //dt,
    //
    int idirx,idiry,idirz;
    double CellTime,tx,ty,tz;

    //time to reach border in x-direction
    
    if(vx > 0)
    {
        idirx = DOWN;
        tx=(lx[icell + 1] - x)/ vx;
    }

    else if(vx < 0)
    {
        idirx = UP;
        tx=-(x-lx[icell])/ vx;
    }

    else
    {
        idirx = UP;
        tx = 1e99;
    }

    //time to reach border in y-direction
    if(vy>0)
    {
        idiry=RIGHT;
        ty=(ly[jcell+1]-y)/vy;
    }

    else if(vy<0)
    {
        idiry=LEFT;
        ty=-(y-ly[jcell])/vy;
    }

    else
    {
        idiry=RIGHT;
        ty=1e99;
    }

    //time to reach border in y-direction
    if(vz>0)
    {
        idirz= BACK;
        tz=(lz[kcell+1]-z)/vz;
    }

    else if(vz<0)
    {
        idirz=FRONT;
        tz=-(z-lz[kcell])/vz;
    }

    else
    {
        idirz=FRONT;
        tz=1e99;
    }

    //chose the smallest time
    idir = idirx;
    CellTime = tx;

    if(ty < CellTime)
    {
        idir = idiry;
        CellTime = ty;
    }

    if(tz < CellTime)
    {
        idir = idirz;
        CellTime = tz;
    }

    if(CellTime < - MY_ZERO)
      {
	cout<<"CellTime<0"<<"   "<<par_id << ' ' << par_type << ' '
	  << CellTime << ' ' << icell<<"  "<<jcell<<"  " << x<<"  "<< y << "  " <<endl;
	CellTime=0;
      }

    return CellTime;
}


int MeshQuantities::HitCell()
{
    //Model has arrived at a boundary. Perform necessary action
  
    int imot;

    //get rules of motion for boundary
    imot=(*c_motion_rule)[C_LINDEX_GHOST(icell, jcell, kcell, idir, 6)];

    //pass border
    if(imot == PASS || imot == GENERATE || imot == GENREF) {
      if (((idir == LEFT) && (jcell == c_jbegin_ghost)) ||
          ((idir == RIGHT) && (jcell == c_jend_ghost)))
        return 0;
    }
    
    if (imot == PASS) {
      Pass();
      band.flag_cellhit=1;
    }
    
    // ..reflecte at border
    else if(imot==REFLECT){
        Reflect();
        band.flag_cellhit=2;
    }
    // ..reflecte or scatter at border
    else if(imot==SCATTOX)
      {
        band.flag_cellhit=3;
    	if(Flag_SurfaceScatter)
          Reflect();

    	else
          {
            //GETINTER();
        if(Random()>band.difpr[par_type])
            Reflect();
        else
            Diffuse();
        }
    }
    //periodic boundary condition in quadrant
    else if(imot==PERIOD){
      Period();
    }
    //periodic boundary condition in quadrant plus generation of a new particle
    else if(imot==GENERATE)
    {	
      //		GetCurrent(patch);//
        Generate();
	Reflect();
/*
     int icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];

     icont --;

      if (icont < 0){
	err_message(WRONG_CONTACT, "Attached_Contact wrong , in GENERATE.");
      }

     contact[icont].NumParGen ++ ;
     contact[icont].CharGen += charge;
     contact[icont].EnergyGen += energy;

     gen_par ++;
     Pass(); 
     */
    }
    //periodic boundary condition in quadrant plus generation of a new particle
    else if(imot==GENREF)
    {
      //		GetCurrent(patch);//
        
        Generate();

        Reflect();
    }

    //catch particle at contact
    else if(imot==CATCH)
    {
      //      GetCurrent(patch);
      CatchAtContact();
      Flag_Catch=true;
    }

    //catch particle at a gate contact
    else if(imot==CATCHGATE)
    {
      //GetCurrent(patch);
        CatchAtGate();
	Flag_Catch=true;
    }
    else
    {
      cout<<"error QUADHIT: No such rule of motion "<< mpi_rank << ' ' << imot << ' ' << p_box 
	  << ' ' << icell << ' ' << jcell << ' ' << kcell << endl;
        exit(0);
    }

    Rho = (*c_par_charge)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
  
    DA = (*c_da)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
    
    if (par_type == 0) {
      Ex  = ChargeSign[par_type] * (*c_field_x)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
	
      Ey  = ChargeSign[par_type] * (*c_field_y)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];

      Ez  = ChargeSign[par_type] * (*c_field_z)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
    } else {
      Ex  = ChargeSign[par_type] * (*c_h_field_x)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
	
      Ey  = ChargeSign[par_type] * (*c_h_field_y)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];

      Ez  = ChargeSign[par_type] * (*c_h_field_z)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];
    }

    return 1;
}


void MeshQuantities::Reflect()
{
    // 解析能带：翻转 k 分量，更新索引并查表速度
if (band.use_analytic_band) {
        // 1. 处理 X 方向反射
        if (idir == UP) { // 撞击左边界 (Lx[i])
            kx = -kx;
            vx = -vx;     // 暂时手动翻转
            x = lx[icell]; // [Fix] 强制吸附到边界
        } 
        else if (idir == DOWN) { // 撞击右边界 (Lx[i+1])
            kx = -kx;
            vx = -vx;
            x = lx[icell + 1]; // [Fix] 强制吸附到边界
        } 
        
        // 2. 处理 Y 方向反射
        else if (idir == LEFT) { // 撞击左边界 (Ly[j])
            ky = -ky;
            vy = -vy;
            y = ly[jcell]; // [Fix] 强制吸附到边界
        } 
        else if (idir == RIGHT) { // 撞击右边界 (Ly[j+1])
            ky = -ky;
            vy = -vy;
            y = ly[jcell + 1]; // [Fix] 强制吸附到边界
        } 
        
        // 3. 处理 Z 方向反射
        else if (idir == FRONT) { // 撞击前边界 (Lz[k])
            kz = -kz;
            vz = -vz;
            z = lz[kcell]; // [Fix] 强制吸附到边界
        } 
        else if (idir == BACK) { // 撞击后边界 (Lz[k+1])
            kz = -kz;
            vz = -vz;
            z = lz[kcell + 1]; // [Fix] 强制吸附到边界
        }

        static const double a_lattice = 5.43e-10;
        double to_pi = 1.0 / ((PI / a_lattice) * spr0);
        par_iter->kx_idx = band.GetAxisIndex_O1(kx * to_pi);
        par_iter->ky_idx = band.GetAxisIndex_O1(ky * to_pi);
        par_iter->kz_idx = band.GetAxisIndex_O1(kz * to_pi);

       //band.GetAnalyticV_By_Index(&(*par_iter));
       //vx = band.analytic_vx;
       //vy = band.analytic_vy;
       //vz = band.analytic_vz;
        return;
    }

    //  Refelct the particle at the boundary
    if(idir==UP||idir==DOWN)
    {
        kx=-kx;
        isym=band.indmat[-band.matsym[3][isym]+1][band.matsym[4][isym]+1]
            [band.matsym[5][isym]+1][band.matsym[0][isym]-1]
            [band.matsym[1][isym]-1][band.matsym[2][isym]-1];

    }
    else if(idir==RIGHT||idir==LEFT)
    {
        ky=-ky;
        isym=band.indmat[band.matsym[3][isym]+1][-band.matsym[4][isym]+1]
            [band.matsym[5][isym]+1][band.matsym[0][isym]-1]
            [band.matsym[1][isym]-1][band.matsym[2][isym]-1];
    }
    else if(idir==FRONT||idir==BACK)
    {
        kz=-kz;
        isym=band.indmat[band.matsym[3][isym]+1][band.matsym[4][isym]+1]
            [-band.matsym[5][isym]+1][band.matsym[0][isym]-1]
            [band.matsym[1][isym]-1][band.matsym[2][isym]-1];
    }

}


void MeshQuantities::Diffuse()
{
    //解析能带分支：等能面随机选取状态并保证朝内
    if (band.use_analytic_band) {
        Particle tmp_p{};
        tmp_p.energy = energy;
        band.SelectAnalyticKState(&tmp_p, energy);
        band.GetAnalyticV_By_Index(&tmp_p);

        double vx_new = band.analytic_vx;
        double vy_new =  band.analytic_vy;
        double vz_new = band.analytic_vz;

        bool need_flip = false;
        if (idir == UP    && vx_new < 0) need_flip = true;
        if (idir == DOWN  && vx_new > 0) need_flip = true;
        if (idir == LEFT  && vy_new < 0) need_flip = true;
        if (idir == RIGHT && vy_new > 0) need_flip = true;
        if (idir == FRONT && vz_new < 0) need_flip = true;
        if (idir == BACK  && vz_new > 0) need_flip = true;

        if (need_flip) {
            tmp_p.kx = -tmp_p.kx; tmp_p.ky = -tmp_p.ky; tmp_p.kz = -tmp_p.kz;
            vx_new = -vx_new; vy_new = -vy_new; vz_new = -vz_new;

            // 翻转后重算索引并查表
            static const double a_lattice = 5.43e-10;
            double to_pi = 1.0 / ((PI / a_lattice) * spr0);
            tmp_p.kx_idx = band.GetAxisIndex_O1(tmp_p.kx * to_pi);
            tmp_p.ky_idx = band.GetAxisIndex_O1(tmp_p.ky * to_pi);
            tmp_p.kz_idx = band.GetAxisIndex_O1(tmp_p.kz * to_pi);
            band.GetAnalyticV_By_Index(&tmp_p);
            vx_new = band.analytic_vx;
            vy_new = band.analytic_vy;
            vz_new = band.analytic_vz;
        }

        // 同步状态
        kx = tmp_p.kx; ky = tmp_p.ky; kz = tmp_p.kz;
        vx = vx_new; vy = vy_new; vz = vz_new;
        par_iter->kx = kx; par_iter->ky = ky; par_iter->kz = kz;
        par_iter->kx_idx = tmp_p.kx_idx;
        par_iter->ky_idx = tmp_p.ky_idx;
        par_iter->kz_idx = tmp_p.kz_idx;
        par_iter->energy = tmp_p.energy;

        //// 轻推粒子离开边界，避免 CellTf=0 死锁
        //double eps_nudge = 1.0e-9;
        //if (idir == UP) x += eps_nudge;
        //else if (idir == DOWN) x -= eps_nudge;
        //else if (idir == LEFT) y += eps_nudge;
        //else if (idir == RIGHT) y -= eps_nudge;
        //else if (idir == FRONT) z += eps_nudge;
        //else if (idir == BACK) z -= eps_nudge;
        return;
    }

    //diffusive boundary condition

    int ibandold,itetold,isymold;
    double xkold,ykold,zkold;

    // ..save initial particle state
    ibandold=iband;
    itetold=itet;
    isymold=isym;
    xkold=kx;
    ykold=ky;
    zkold=kz;
    // ..chose random direction on equi energy surface
    TALIND(1);
    get_state();
    // ..if GetState fails restore inital particle state
    if(Flag_SelfScatter)
    {
        iband=ibandold;
        itet=itetold;
        isym=isymold;
        kx=xkold;
        ky=ykold;
        kz=zkold;
    }
    GetV();
    // ..chose particle velocity to point into actual device region
    if(idir==UP)
    {
        if(vx==0.0)
        {
            iband=ibandold;
            itet=itetold;
            isym=isymold;
            kx=xkold;
            ky=ykold;
            kz=zkold;
        }
        if(vx<0.0) //vx
        {
            kx=-kx;
            isym=band.indmat[-band.matsym[3][isym]+1][band.matsym[4][isym]+1]
                [band.matsym[5][isym]+1][band.matsym[0][isym]-1]
                [band.matsym[1][isym]-1][band.matsym[2][isym]-1];
        }
    }
    else if(idir==RIGHT)
    {
        if(vy==0.0)
        {
            iband=ibandold;
            itet=itetold;
            isym=isymold;
            kx=xkold;
            ky=ykold;
            kz=zkold;
        }
        if(vy>0.0)//vx
        {
            ky=-ky;
            isym=band.indmat[band.matsym[3][isym]+1][-band.matsym[4][isym]+1]
                [band.matsym[5][isym]+1][band.matsym[0][isym]-1]
                [band.matsym[1][isym]-1][band.matsym[2][isym]-1];
        }
    }
    else if(idir==DOWN)
    {
        if(vx==0.0)
        {
            iband=ibandold;
            itet=itetold;
            isym=isymold;
            kx=xkold;
            ky=ykold;
            kz=zkold;
        }
        if(vx>0.0)
        {
            kx=-kx;
            isym=band.indmat[-band.matsym[3][isym]+1][band.matsym[4][isym]+1]
                [band.matsym[5][isym]+1][band.matsym[0][isym]-1]
                [band.matsym[1][isym]-1][band.matsym[2][isym]-1];
        }
    }
    else if(idir==LEFT)
    {
        if(vy==0.0)
        {
            iband=ibandold;
            itet=itetold;
            isym=isymold;
            kx=xkold;
            ky=ykold;
            kz=zkold;
        }
        if(vy<0.0)
        {
            ky=-ky;
            isym=band.indmat[band.matsym[3][isym]+1][-band.matsym[4][isym]+1]
                [band.matsym[5][isym]+1][band.matsym[0][isym]-1]
                [band.matsym[1][isym]-1][band.matsym[2][isym]-1];
        }
    }
    else if(idir==FRONT)
    {
        if(vz==0.0)
        {
            iband=ibandold;
            itet=itetold;
            isym=isymold;
            kx=xkold;
            ky=ykold;
            kz=zkold;
        }
        if(vz<0.0)
        {
            kz=-kz;
            isym=band.indmat[band.matsym[3][isym]+1][band.matsym[4][isym]+1]
                [-band.matsym[5][isym]+1][band.matsym[0][isym]-1]
                [band.matsym[1][isym]-1][band.matsym[2][isym]-1];
        }
    }
    else if(idir==BACK)
    {
        if(vz==0.0)
        {
            iband=ibandold;
            itet=itetold;
            isym=isymold;
            kx=xkold;
            ky=ykold;
            kz=zkold;
       }
        if(vz>0.0)
        {
            kz=-kz;
            isym=band.indmat[band.matsym[3][isym]+1][band.matsym[4][isym]+1]
                [-band.matsym[5][isym]+1][band.matsym[0][isym]-1]
                [band.matsym[1][isym]-1][band.matsym[2][isym]-1];
        }
    }
}


void MeshQuantities::Pass()
{
    switch(idir)
    {
        case UP   :x=lx[icell];   icell--; break;
        case DOWN :x=lx[icell+1]; icell++; break;
        case RIGHT:{
          y = ly[jcell + 1];
          jcell++;
          rcurrent[jcell] += charge;
          break;
        }
        case LEFT :{
          // use current cell index for position and current accumulation,
          // then move to the left neighbor
          int jcur = jcell;
          y = ly[jcur];
          lcurrent[jcur] += charge;
          jcell--;
          break;
        }
        case FRONT:z=lz[kcell];   kcell--; break;
        case BACK :z=lz[kcell+1]; kcell++; break;
    }
}


void MeshQuantities::Period()
{
    //  periodic boundary condition
    if(idir==UP)
        x=lx[icell+1];
    else if(idir==DOWN)
        x=lx[icell];
    else if(idir==LEFT)
        y=ly[jcell+1];
    else if(idir==RIGHT)
        y=ly[jcell];
    else if(idir==FRONT)
        z=lz[kcell + 1];
    else if(idir==BACK)
        z=lz[kcell];

}


void MeshQuantities::Generate()
{
    //  generation boundary condition

    int icont;
    gen_par ++;
  
    icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell, kcell)];

    icont --;

    if (icont < 0){
      err_message(WRONG_CONTACT, "Attached_Contact wrong , in GENERATE.");
    }

    contact[icont].NumParGen++;
    contact[icont].CharGen += charge;
    contact[icont].EnergyGen += energy;

    if (idir == RIGHT) 
      rcurrent[jcell + 1] += charge;

    if (idir == LEFT) 
      lcurrent[jcell] += charge;

    Particle tmp_par;

    list<Particle> * gen_par_list;
    int ii,jj,kk;
    ii = icell;
    jj = jcell;
    kk = kcell;
    switch (idir) {
      case UP:
	gen_par_list = &par_list[C_LINDEX_GHOST_ONE(icell - 1  , jcell, kcell)];
	ii = icell - 1;
	break;
      case DOWN:
	gen_par_list = &par_list[C_LINDEX_GHOST_ONE(icell + 1  , jcell, kcell)];
	ii = icell + 1;
	break;
      case LEFT:
	gen_par_list = &par_list[C_LINDEX_GHOST_ONE(icell , jcell - 1, kcell)];
        jj = jcell - 1;
	break;
      case RIGHT:
	gen_par_list = &par_list[C_LINDEX_GHOST_ONE(icell , jcell + 1, kcell)];
        jj = jcell + 1;
	break;
      case FRONT: 
	gen_par_list = &par_list[C_LINDEX_GHOST_ONE(icell , jcell , kcell - 1)];
        kk = kcell - 1;
	break;
      case BACK: 
	gen_par_list = &par_list[C_LINDEX_GHOST_ONE(icell , jcell , kcell + 1)];
        kk = kcell + 1;
	break;
    }

    list<Particle>::iterator
      tmp_iter = gen_par_list->insert(gen_par_list->end(), tmp_par);
       
    OutPar(tmp_iter);
    tmp_iter->i = ii;
    tmp_iter->j = jj;
    tmp_iter->k = kk;
}


void MeshQuantities::CatchAtContact()
{
    //catch boundary condition
    
  int icont;

  /* assume we are simulating the MOSFET whose structure is simple */

  // 根据方向获取 相邻 cell 所属的 contact 索引
  if(idir==UP)  // Negative x direction
    icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell - 1, jcell,kcell)];
  else if(idir==DOWN) // Positive x direction
    icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell + 1, jcell,kcell)];
  else if(idir==LEFT) // Negative y direction
    icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell - 1,kcell)];
  else if(idir==RIGHT)  // Positive y direction
    icont= (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell + 1,kcell)];
  else if(idir==FRONT)  // Negative z direction
    icont= (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell,kcell - 1)];
  else if(idir==BACK)   // Positive z direction
    icont= (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell,kcell + 1)];

//  icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell,kcell)];
  // 调整 contact 的索引
  icont --;
  
  if (icont < 0){
    err_message(WRONG_CONTACT, "Attached_Contact wrong , in catchAtContact.");
  }

  // 粒子在y方向移动时，统计电流
  if (idir == RIGHT) 
    rcurrent[jcell + 1] += charge;  // particle 向右移动时，右边cell多一个电荷
  if (idir == LEFT) 
    lcurrent[jcell] += charge;      // particle 向左移动时，当前 cell 多一个电荷

  contact[icont].NumParCatch++;       // 当前 contact 上捕获的粒子数 + 1
  contact[icont].CharCatch += charge; // 更新当前 contact 上捕获的电荷量
  contact[icont].EnergyGen += energy; // 更新当前 contact 上因捕获粒子而产生的能量
}

void MeshQuantities::CatchAtGate()
{
	//catchgate boundary condition
	//
  int icont;

   if(idir==UP)
    icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell - 1, jcell,kcell)];
  else if(idir==DOWN)
    icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell + 1, jcell,kcell)];
  else if(idir==LEFT)
    icont = (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell - 1,kcell)];
  else if(idir==RIGHT)
    icont= (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell + 1,kcell)];
  else if(idir==FRONT)
    icont= (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell,kcell - 1)];
  else if(idir==BACK)
    icont= (*c_attached_contact)[C_LINDEX_GHOST_ONE(icell, jcell,kcell + 1)];

  icont --;
  
  if (icont < 0){
    err_message(WRONG_CONTACT, "Attached_Contact wrong , in catchAtContact.");
  }

  contact[icont].NumParCatch++;
  contact[icont].CharCatch+=charge;
  contact[icont].EnergyGen+=energy;
    
}
