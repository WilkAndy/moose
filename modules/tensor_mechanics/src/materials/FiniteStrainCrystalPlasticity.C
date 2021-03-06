#include "FiniteStrainCrystalPlasticity.h"
#include <cmath>

extern "C" void FORTRAN_CALL(dsyev) ( ... );

template<>
InputParameters validParams<FiniteStrainCrystalPlasticity>()
{
  InputParameters params = validParams<FiniteStrainMaterial>();

  params.addRequiredParam< int >("nss","Properties");
  params.addRequiredParam< std::vector<Real> >("gprops","Properties");
  params.addRequiredParam< std::vector<Real> >("hprops","Properties");
  params.addRequiredParam< std::vector<Real> >("flowprops","Properties");
  params.addRequiredParam<std::string>("slip_sys_file_name", "Name of the file containing the slip system");
  params.addParam<std::string>("euler_angle_file_name","", "Name of the file containing the euler angles");
  params.addParam<Real>("rtol",1e-8,"Constitutive stress residue tolerance");
  params.addParam<Real>("gtol",1e2,"Constitutive gss residue tolerance");
  params.addParam<Real>("slip_incr_tol",2e-2,"Constitutive gss residue tolerance");

  return params;
}

FiniteStrainCrystalPlasticity::FiniteStrainCrystalPlasticity(const std::string & name,
                                                             InputParameters parameters)
    : FiniteStrainMaterial(name, parameters),
      _nss(getParam<int>("nss")),
      _gprops(getParam<std::vector<Real> >("gprops")),
      _hprops(getParam<std::vector<Real> >("hprops")),
      _flowprops(getParam<std::vector<Real> >("flowprops")),
      _slip_sys_file_name(getParam<std::string>("slip_sys_file_name")),
      _euler_angle_file_name(getParam<std::string>("euler_angle_file_name")),
      _rtol(getParam<Real>("rtol")),
      _gtol(getParam<Real>("gtol")),
  _slip_incr_tol(getParam<Real>("slip_incr_tol")),
  _fp(declareProperty<RankTwoTensor>("fp")),
  _fp_old(declarePropertyOld<RankTwoTensor>("fp")),
  _pk2(declareProperty<RankTwoTensor>("pk2")),
  _pk2_old(declarePropertyOld<RankTwoTensor>("pk2")),
  _lag_e(declareProperty<RankTwoTensor>("lage")),
  _gss(declareProperty<std::vector<Real> >("gss")),
  _gss_old(declarePropertyOld<std::vector<Real> >("gss")),
  _acc_slip(declareProperty<Real>("acc_slip")),
  _acc_slip_old(declarePropertyOld<Real>("acc_slip")),
  _update_rot(declareProperty<RankTwoTensor>("update_rot")),
  _crysrot(declareProperty<RankTwoTensor>("crysrot")),
  _crysrot_old(declarePropertyOld<RankTwoTensor>("crysrot"))


{
}

void FiniteStrainCrystalPlasticity::initQpStatefulProperties()
{
  int ind,is,ie,flag;
  Real *data;
  RealTensorValue rot;

  _stress[_qp].zero();

  _fp[_qp].zero();
  _fp[_qp].addIa(1.0);

  _pk2[_qp].zero();

  _gss[_qp].resize(_nss);
  _gss_old[_qp].resize(_nss);


  _acc_slip[_qp]=0.0;

  _a0.resize(_nss);
  _xm.resize(_nss);

  data=_gprops.data();

  flag=0;
  ind=-1;
  while (flag==0)
  {
    ind=ind+1;
    is=data[ind];
    ind=ind+1;
    ie=data[ind];

    if (ie==_nss)
      flag=1;

    ind=ind+1;
    for (int i=is;i<=ie;i++)
    {
      _gss[_qp][i-1]=_gprops[ind];
    }
  }


  data=_hprops.data();

  _r=data[0];
  _h0=data[1];
  _tau_init=data[2];
  _tau_sat=data[3];

  data=_flowprops.data();


  flag=0;
  ind=-1;
  while (flag==0)
  {
    ind=ind+1;
    is=data[ind];
    ind=ind+1;
    ie=data[ind];

    if (ie==_nss)
      flag=1;

    for (int i=is;i<=ie;i++)
    {
      _a0[i-1]=_flowprops[ind+1];
      _xm[i-1]=_flowprops[ind+2];
    }

    ind=ind+2;
  }


  get_euler_ang();
  get_euler_rot();

  _mo.resize(_nss*3);
  _no.resize(_nss*3);

  get_slip_sys();

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      rot(i,j)=_crysrot[_qp](i,j);


  _elasticity_tensor[_qp].rotate(rot);

  _update_rot[_qp]=_crysrot[_qp];


}




void FiniteStrainCrystalPlasticity::computeQpStress()
{

  std::vector<Real> gss_prev(_nss), slip_incr(_nss), tau(_nss);
  RankTwoTensor pk2,fp_old_inv,fp_inv,resid,dpk2,sig,rot;
  Real gmax,gdiff,rnorm;
  RankFourTensor jac;
  Real slip_incr_max,fac;
  int iter,iterg,maxiter,maxiterg;


  maxiter=100;
  maxiterg=100;

  _crysrot[_qp]=_crysrot_old[_qp];

  fp_inv.zero();
  fp_old_inv=_fp_old[_qp].inverse();


  gmax=1.1*_gtol;

  for (int i=0;i<_nss;i++)
    _gss[_qp][i]=_gss_old[_qp][i];

  iterg=0;
  while (gmax > _gtol && iterg < maxiterg)
  {

    iter=0;

    pk2=_pk2_old[_qp];

    calc_resid_jacob(&pk2,&sig,&fp_old_inv,&fp_inv,&slip_incr[0],&tau[0],&resid,&jac);

    slip_incr_max=0.0;

    for (int i=0;i<_nss;i++)
    {
      if (fabs(slip_incr[i])>slip_incr_max)
      {
        slip_incr_max=fabs(slip_incr[i]);
      }

    }

    fac=1.0;
    if (slip_incr_max > _slip_incr_tol)
    {
      fac=0.1;
      printf("slip_incr=%d %d %f\n",_current_elem->id(),_qp,slip_incr_max);
      mooseError("Slip increment exceeds tolerance");
    }

    rnorm=resid.L2norm();


    //      printf("rnorm=%d %f\n",iter,rnorm);

    while (rnorm > _rtol && iter <  maxiter)
    {

      dpk2=jac.invSymm()*resid;

      pk2=pk2-dpk2*fac;

      calc_resid_jacob(&pk2,&sig,&fp_old_inv,&fp_inv,&slip_incr[0],&tau[0],&resid,&jac);

      slip_incr_max=0.0;

      for (int i=0;i<_nss;i++)
      {
        if (fabs(slip_incr[i])>slip_incr_max)
        {
          slip_incr_max=fabs(slip_incr[i]);
        }

      }

      fac=1.0;
      if (slip_incr_max > _slip_incr_tol)
      {
        fac=0.1;
        printf("slip_incr=%d %d %f\n",_current_elem->id(),_qp,slip_incr_max);
        mooseError("Slip increment exceeds tolerance");
      }

      rnorm=resid.L2norm();

      iter+=1;


      //    printf("rnorm=%d %f\n",iter,rnorm);

    }

    if (iter==maxiter)
      mooseError("Stress Integration error \n");


    for (int i=0;i<_nss;i++)
      gss_prev[i]=_gss[_qp][i];

    update_gss(&slip_incr[0]);

    gmax=0.0;
    for (int i=0;i<_nss;i++)
    {
      gdiff=fabs(gss_prev[i]-_gss[_qp][i]);

      if (gdiff > gmax)
        gmax=gdiff;

    }

    iterg+=1;

  }

  if (iterg==maxiterg)
  {
    printf("gmax=%f\n",gmax);
    mooseError("Hardness Integration error \n");
  }


  _fp[_qp]=fp_inv.inverse();
  _pk2[_qp]=pk2;
  _stress[_qp]=sig;


  rot=getmatrot(_dfgrd[_qp]);

  _update_rot[_qp]=rot*_crysrot[_qp];


}


void
FiniteStrainCrystalPlasticity::get_euler_rot()
{

  Real phi1, phi, phi2;
  Real cp, cp1, cp2, sp, sp1, sp2;
  RankTwoTensor RT;
  Real pi = 4.0*atan(1.0);

  //  printf("get euler=%d %d %f %f %f\n",_current_elem->id(),_qp,_euler_angle_1,_euler_angle_2,_euler_angle_3);

  phi1 = _euler_angle_1 * (pi/180.0);
  phi = _euler_angle_2 * (pi/180.0);
  phi2 = _euler_angle_3 * (pi/180.0);

  cp1 = cos(phi1);
  cp2 = cos(phi2);
  cp = cos(phi);

  sp1 = sin(phi1);
  sp2 = sin(phi2);
  sp = sin(phi);

  RT(0,0) = cp1 * cp2 - sp1 * sp2 * cp;
  RT(0,1) = sp1 * cp2 + cp1 * sp2 * cp;
  RT(0,2) = sp2 * sp;
  RT(1,0) = -cp1 * sp2 - sp1 * cp2 * cp;
  RT(1,1) = -sp1 * sp2 + cp1 * cp2 * cp;
  RT(1,2) = cp2 * sp;
  RT(2,0) = sp1 * sp;
  RT(2,1) = -cp1 * sp;
  RT(2,2) = cp;

  _crysrot[_qp]= RT.transpose();
  _crysrot_old[_qp]=_crysrot[_qp];

}


void
FiniteStrainCrystalPlasticity::get_slip_sys()
{

  std::vector<Real> sd(3*_nss), sn(3*_nss);
  Real vec[3];
  std::ifstream fileslipsys;

  fileslipsys.open(_slip_sys_file_name.c_str());

  if (!fileslipsys)
    mooseError("Can't open slip system input file ");

  for (int i=0;i<_nss;i++)
  {
    for (int j=0;j<3;j++)
      fileslipsys >> vec[j];

    Real mag;
    mag=pow(vec[0],2)+pow(vec[1],2)+pow(vec[2],2);
    mag=pow(mag,0.5);

    for (int j=0;j<3;j++)
      sn[i*3+j]=vec[j]/mag;


  }

  for (int i=0;i<_nss;i++)
  {
    for (int j=0;j<3;j++)
      fileslipsys >> vec[j];



    Real mag;
    mag=pow(vec[0],2)+pow(vec[1],2)+pow(vec[2],2);
    mag=pow(mag,0.5);

    for (int j=0;j<3;j++)
      sd[i*3+j]=vec[j]/mag;


  }

  fileslipsys.close();

  for (int i=0;i<_nss;i++)
  {

    for (int j=0;j<3;j++)
    {
      _mo[i*3+j]=0.0;
      for (int k=0;k<3;k++)
        _mo[i*3+j]=_mo[i*3+j]+_crysrot[_qp](j,k)*sd[i*3+k];

    }

    for (int j=0;j<3;j++)
    {
      _no[i*3+j]=0.0;
      for (int k=0;k<3;k++)
        _no[i*3+j]=_no[i*3+j]+_crysrot[_qp](j,k)*sn[i*3+k];

    }

  }


}

void
FiniteStrainCrystalPlasticity::update_gss(Real *slip_incr)
{
  std::vector<Real> hb(_nss);
  Real qab,val;
  Real *data;//Kalidindi

  data=_hprops.data();//Kalidindi
  Real a=data[4];//Kalidindi

  _acc_slip[_qp]=_acc_slip_old[_qp];
  for (int i=0;i < _nss; i++)
    _acc_slip[_qp]=_acc_slip[_qp]+fabs(slip_incr[i]);


  val=cosh(_h0*_acc_slip[_qp]/(_tau_sat-_tau_init));//Karthik
  val=_h0*pow(1.0/val,2.0);//Kalidindi


  for (int i=0;i < _nss; i++)
    //    hb[i]=val;
    hb[i]=_h0*pow(1.0-_gss[_qp][i]/_tau_sat,a);



  for (int i=0;i < _nss; i++)
  {

    _gss[_qp][i]=_gss_old[_qp][i];

    for (int j=0;j < _nss; j++)
    {
      int iplane,jplane;
      iplane=i/3;
      jplane=j/3;

      if (iplane==jplane)//Kalidindi
        qab=1.0;
      else
        qab=_r;


      _gss[_qp][i]=_gss[_qp][i]+qab*hb[j]*fabs(slip_incr[j]);


    }
  }

}


void
FiniteStrainCrystalPlasticity::calc_resid_jacob(RankTwoTensor* pk2,RankTwoTensor *sig, RankTwoTensor* fp_old_inv,RankTwoTensor* fp_inv,
                                                Real* slip_incr,Real *tau, RankTwoTensor* resid, RankFourTensor *jac)
{

  RankTwoTensor fe,ce,ee,iden,ce_pk2;
  RankTwoTensor eqv_slip_incr,pk2_new;
  std::vector<RankTwoTensor> s0(_nss),dtaudpk2(_nss),dfpinvdslip(_nss);
  RankTwoTensor temp2;
  std::vector<Real> dslipdtau(_nss);
  RankFourTensor dfedfpinv,deedfe,dfpinvdpk2,idenFour;
  RankFourTensor temp4;


  /*Calculate Residual*/
  iden.zero();
  iden.addIa(1.0);

  fe=_dfgrd[_qp]*(*fp_old_inv);

  ce=fe.transpose()*fe;
  ee=ce-iden;
  ee*=0.5;

  ce_pk2=ce*(*pk2);
  ce_pk2=ce_pk2/fe.det();
  //  ce_pk2=(*pk2);//Approximation


  for (int i=0;i<_nss;i++)
  {
    for (int j=0;j<3;j++)
      for (int k=0;k<3;k++)
        s0[i](j,k)=_mo[i*3+j]*_no[i*3+k];

    tau[i]=ce_pk2.doubleContraction(s0[i]);

  }


  get_slip_incr(tau,slip_incr,&dslipdtau[0]);//Calculate dslip,dslipdtau

  eqv_slip_incr.zero();
  for (int i=0;i<_nss;i++)
    eqv_slip_incr+=s0[i]*slip_incr[i];

  eqv_slip_incr=iden-eqv_slip_incr;

  (*fp_inv)=(*fp_old_inv)*eqv_slip_incr;

  fe=_dfgrd[_qp]*(*fp_inv);

  ce=fe.transpose()*fe;
  ee=ce-iden;
  ee*=0.5;

  pk2_new=_elasticity_tensor[_qp]*ee;

  _lag_e[_qp]=_dfgrd[_qp].transpose()*_dfgrd[_qp]-iden;
  _lag_e[_qp]=_lag_e[_qp]*0.5;

  (*resid)=(*pk2)-pk2_new;
  /*End Calculate Residual*/
  /*Calculate Jacobian*/

  for (int i=0;i<_nss;i++)
  {
    dtaudpk2[i]=s0[i];
    dfpinvdslip[i]=-(*fp_old_inv)*s0[i];
  }

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<3;k++)
        dfedfpinv(i,j,k,j)=_dfgrd[_qp](i,k);

  deedfe.zero();

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<3;k++)
      {
        deedfe(i,j,k,i)=deedfe(i,j,k,i)+fe(k,j)*0.5;
        deedfe(i,j,k,j)=deedfe(i,j,k,j)+fe(k,i)*0.5;
      }

  dfpinvdpk2.zero();
  for (int i=0;i<_nss;i++)
  {
    temp2=dfpinvdslip[i]*dslipdtau[i];
    temp4=outerProduct(temp2,dtaudpk2[i]);
    dfpinvdpk2+=temp4;

  }



  (*jac)=_elasticity_tensor[_qp]*deedfe*dfedfpinv*dfpinvdpk2;


  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      idenFour(i,j,i,j)=1.0;

  (*jac)=idenFour-(*jac);


  /*End Calculate Jacobian*/

  (*sig)=fe*(*pk2)*fe.transpose();
  (*sig)=(*sig)/fe.det();



}

void
FiniteStrainCrystalPlasticity::get_slip_incr(Real* tau,Real *slip_incr,Real *dslipdtau)
{
  for (int i=0;i<_nss;i++)
  {
    slip_incr[i]=_a0[i]*pow(fabs(tau[i]/_gss[_qp][i]),1.0/_xm[i])*copysign(1.0,tau[i])*_dt;

  }

  for (int i=0;i<_nss;i++)
  {
    dslipdtau[i]=_a0[i]/_xm[i]*pow(fabs(tau[i]/_gss[_qp][i]),1.0/_xm[i]-1.0)/_gss[_qp][i]*_dt;
  }

}

RankFourTensor FiniteStrainCrystalPlasticity::outerProduct(RankTwoTensor &a, RankTwoTensor &b)
{
  RankFourTensor result;

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<3;k++)
        for (int l=0;l<3;l++)
          result(i,j,k,l)=a(i,j)*b(k,l);

  return result;

}


RankTwoTensor FiniteStrainCrystalPlasticity::getmatrot(RankTwoTensor &a)
{
  RankTwoTensor rot;
  RankTwoTensor c,diag,evec;
  double cmat[3][3],w[3],work[10];
  int info;
  int nd=3;
  int lwork=10;

  c=a.transpose()*a;

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      cmat[i][j]=c(i,j);

  FORTRAN_CALL(dsyev)("V","U",&nd,cmat,&nd,&w,&work,&lwork,&info);

  diag.zero();

  for (int i=0;i<3;i++)
    diag(i,i)=pow(w[i],0.5);

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      evec(i,j)=cmat[i][j];

  rot=a*((evec.transpose()*diag*evec).inverse());


  return rot;

}


void
FiniteStrainCrystalPlasticity::get_euler_ang()
{

  std::ifstream fileeuler;
  Real vec[3];
  int elemno;

  if (_euler_angle_file_name.length()!=0)
  {

    fileeuler.open(_euler_angle_file_name.c_str());

    if (!fileeuler)
      mooseError("Can't open euler angle input file ");
    else
    {

      while (!fileeuler.eof())
      {
        fileeuler >> elemno;
        for (int i=0;i<3;i++)
          fileeuler >> vec[i];

        if (static_cast<dof_id_type>(elemno-1) == _current_elem->id())
        {
          _euler_angle_1=vec[0];
          _euler_angle_2=vec[1];
          _euler_angle_3=vec[2];

          fileeuler.close();
          return;

        }


      }

      fileeuler.close();

    }

  }


}

void FiniteStrainCrystalPlasticity::computeQpElasticityTensor()
{
  // Fill in the matrix stiffness material property
  _elasticity_tensor[_qp] = _Cijkl;

  ElasticityTensorR4 tan_mod;
  RankTwoTensor fe,fp_inv,pk2fet,fepk2;
  RankFourTensor dfedf,deedfe,dsigdpk2dfe,temp;

  fp_inv=_fp[_qp].inverse();
  fe=_dfgrd[_qp]*fp_inv;

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int l=0;l<3;l++)
        dfedf(i,j,i,l)=fp_inv(l,j);

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<3;k++)
      {
        deedfe(i,j,k,i)=deedfe(i,j,k,i)+fe(k,j)*0.5;
        deedfe(i,j,k,j)=deedfe(i,j,k,j)+fe(k,i)*0.5;
      }


  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<3;k++)
        for (int l=0;l<3;l++)
          temp(i,j,k,l)=fe(i,k)*fe(j,l);

  dsigdpk2dfe=temp*_elasticity_tensor[_qp]*deedfe;

  pk2fet=_pk2[_qp]*fe.transpose();
  fepk2=fe*_pk2[_qp];

  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int l=0;l<3;l++)
      {
        tan_mod(i,j,i,l)=tan_mod(i,j,i,l)+pk2fet(l,j);
        tan_mod(i,j,j,l)=tan_mod(i,j,j,l)+fepk2(i,l);
      }

  tan_mod+=dsigdpk2dfe;

  Real je=fe.det();
  tan_mod/=je;

  _Jacobian_mult[_qp] = tan_mod;

}
