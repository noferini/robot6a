#include "classi.h"
#include "TMath.h"
#include "TFile.h"
#include "TView.h"
#include "TPad.h"

TGeoManager *Geo::mManager = new TGeoManager("Robot","This is my robot");;
TGeoMaterial *Geo::mVacuum = new TGeoMaterial("vacuum",0,0,0);  
TGeoMaterial *Geo::mFe = new TGeoMaterial("Fe",55.845,26,7.87); 
TGeoVolume *Geo::mTop = nullptr;
TGeoMedium *Geo::mAir = new TGeoMedium("Vacuum",0,mVacuum);
TGeoMedium *Geo::mIron = new TGeoMedium("Iron",1,mFe);

//_________________________________
Link::Link(){
  mMatrice = new TMatrixT<double>(4,4);
  mMatriceTrunc = new TMatrixT<double>(4,4);
  mMatriceInv = new TMatrixT<double>(4,4);
}
//_________________________________
Link::~Link(){
  if(mMatrice) delete mMatrice; 
  if(mMatriceTrunc) delete mMatriceTrunc; 
  if(mMatriceInv) delete mMatriceInv; 
}
//_________________________________
void Link::setParameters(double length, double twist, double distance, double theta){
  double ctheta = cos(theta*TMath::DegToRad());
  double stheta = sin(theta*TMath::DegToRad());
  double ctw = cos(twist*TMath::DegToRad());
  double stw = sin(twist*TMath::DegToRad());
  (*mMatrice)[0][0] =  ctheta;
  (*mMatrice)[0][1] = -stheta * ctw;
  (*mMatrice)[0][2] =  stheta * stw;
  (*mMatrice)[0][3] =  ctheta * length;
  (*mMatrice)[1][0] =  stheta;
  (*mMatrice)[1][1] =  ctheta * ctw;
  (*mMatrice)[1][2] = -ctheta * stw;
  (*mMatrice)[1][3] =  stheta * length;
  (*mMatrice)[2][0] = 0;
  (*mMatrice)[2][1] =  stw;
  (*mMatrice)[2][2] =  ctw;
  (*mMatrice)[2][3] = distance;
  (*mMatrice)[3][3] = 1;

  (*mMatriceTrunc)[0][0] =  ctheta;
  (*mMatriceTrunc)[0][1] = -stheta * 0;
  (*mMatriceTrunc)[0][2] =  stheta * 1;
  (*mMatriceTrunc)[0][3] =  ctheta * length;
  (*mMatriceTrunc)[1][0] =  stheta;
  (*mMatriceTrunc)[1][1] =  ctheta * 0;
  (*mMatriceTrunc)[1][2] = -ctheta * 1;
  (*mMatriceTrunc)[1][3] =  stheta * length;
  (*mMatriceTrunc)[2][0] = 0;
  (*mMatriceTrunc)[2][1] =  1*1;
  (*mMatriceTrunc)[2][2] =  1*0;
  (*mMatriceTrunc)[2][3] = distance;
  (*mMatriceTrunc)[3][3] = 1;

  (*mMatriceInv)[0][0] =  ctheta;
  (*mMatriceInv)[0][1] =  stheta;
  (*mMatriceInv)[0][2] =  0;
  (*mMatriceInv)[0][3] =  -length;
  (*mMatriceInv)[1][0] =  -stheta * ctw;
  (*mMatriceInv)[1][1] =  ctheta * ctw;
  (*mMatriceInv)[1][2] =  stw;
  (*mMatriceInv)[1][3] =  -stw * distance;
  (*mMatriceInv)[2][0] =   stheta * stw;
  (*mMatriceInv)[2][1] =  -ctheta * stw;
  (*mMatriceInv)[2][2] =  ctw;
  (*mMatriceInv)[2][3] = -distance * ctw;
  (*mMatriceInv)[3][3] = 1;
}
//_________________________________
void Link::checkInv(){
  TMatrixT<double> mtemp = (*mMatrice) * (*mMatriceInv);
  mtemp.Print();
}
//_________________________________
void Robot6A::setParameters(int ilink, double length, double twist, double distance, double theta){
  mLink[ilink].setParameters(length, twist, distance, theta);
  mLength[ilink]=length;
  mTwist[ilink]=twist;
  mDistance[ilink] = distance; 
  mTheta[ilink]=theta;
}
//_________________________________
Robot6A::Robot6A(int n): mNaxis(n){  
  mMatriceTras = new TMatrixT<double>(4,4);
  (*mMatriceTras)[0][0] = 1;
  (*mMatriceTras)[1][1] = 1;
  (*mMatriceTras)[2][2] = 1;
  (*mMatriceTras)[3][3] = 1;

  mLink = new Link[mNaxis];
  mLength = new float[mNaxis];
  mTwist  = new float[mNaxis];
  mDistance = new float[mNaxis];
  mTheta = new float[mNaxis];

  mMatriceDelta = new TMatrixT<double>*[mNaxis];

  for(int i=0;i < mNaxis; i++){
    mMatriceDelta[i] = new TMatrixT<double>(4,4);
  }
}
//_________________________________
TGeoVolume* Geo::top() {
  if(!mTop){
    mTop = mManager->MakeBox("top",Geo::air(),100,100,100);
    mManager->SetTopVolume(mTop);
    mTop->SetVisibility(0);
  }
 return mTop;
}
//_________________________________
Robot6A::~Robot6A(){
  if(mMatriceTras) delete mMatriceTras;
  for(int i=0;i < mNaxis; i++)
    if(mMatriceDelta[i]) delete mMatriceDelta[i]; 

  delete[] mMatriceDelta;
  delete[] mLink;
  delete[] mLength;
  delete[] mTwist;
  delete[] mDistance;
  delete[] mTheta;
}
//_________________________________
void Robot6A::rotate(int iaxis,double alpha){
  TMatrixT<double> myrot(4,4);
  myrot[0][0] = 1;
  myrot[1][1] = 1;
  myrot[2][2] = 1;
  myrot[3][3] = 1;
  alpha *= TMath::DegToRad();
  if(iaxis == 0){
    myrot[1][1] = cos(alpha);
    myrot[1][2] = -sin(alpha);
    myrot[2][1] = sin(alpha);
    myrot[2][2] = cos(alpha);
  }
  else if(iaxis == 1){
    myrot[2][2] = cos(alpha);
    myrot[2][0] = -sin(alpha);
    myrot[0][2] = sin(alpha);
    myrot[0][0] = cos(alpha);
  }
  else if(iaxis == 2){
    myrot[0][0] = cos(alpha);
    myrot[0][1] = -sin(alpha);
    myrot[1][0] = sin(alpha);
    myrot[1][1] = cos(alpha);
  }
  *mMatriceTras *= myrot;
}
//_________________________________
void Robot6A::setGeometry(){
  if(! Geo::top()->GetNodes()) mFirstNodeInTopVolume = 0;
  else mFirstNodeInTopVolume = Geo::top()->GetNodes()->GetEntries();
  TMatrixT<double> mymatrix(4,4);
  mymatrix = *mMatriceTras;

  // links
  int colors[6] = {2,4,6,8,12,14};
  TGeoVolume *Link[6];
  TGeoVolume *Link2[6];

  TMatrixT<double> tt(4,4);
  tt[0][0]=1;
  tt[1][1]=1;
  tt[2][2]=1;
  tt[3][3]=1;

  TMatrixT<double> tt90(4,4);
  tt90[0][0]=0;
  tt90[0][2]=1;
  tt90[1][1]=1;
  tt90[2][0]=-1;
  tt90[2][2]=0;
  tt90[3][3]=1;

  for(int i=0; i < mNaxis; i++){
    TMatrixT<double> *mat = mLink[i].getMatrix();
    TMatrixT<double> *matTrunc = mLink[i].getMatrixTrunc();

    TMatrixT<double> mm = mymatrix;
    double angleX, angleY,angleZ;
    getAngles(mm, angleX, angleY, angleZ);

    tt[0][3] = 0;//-mDistance[i]*0.5*sin(angleX*TMath::DegToRad());
    tt[1][3] = 0;//-mDistance[i]*0.5*sin(angleX*TMath::DegToRad())*0;
    tt[2][3] = mDistance[i]*0.5;
    // tt[0][3] = -mDistance[i]*0.5*sin(angleX*TMath::DegToRad());
    // tt[1][3] = -mDistance[i]*0.5*sin(angleX*TMath::DegToRad())*0;
    // tt[2][3] = abs(mDistance[i]*0.5*cos(angleX*TMath::DegToRad()));

    mm *= tt;
    
    Link[i] = Geo::manager()->MakeTube(Form("Link%d",i+1),Geo::iron(),0,0.2,mDistance[i]*0.5);
    Link[i]->SetLineColor(colors[i]);
    Link[i]->SetFillColor(colors[i]);

    printf("angles %f %f %f\n",angleX,angleY,angleZ);

    TGeoRotation rot(Form("rot%d",i+1),0,0,0);
    rot.RotateX(angleX);
    rot.RotateY(angleY);
    rot.RotateZ(angleZ);
    TGeoTranslation tra(mm[0][3],mm[1][3],mm[2][3]);
    TGeoCombiTrans *trasf = new TGeoCombiTrans(tra, rot);

    Geo::top()->AddNodeOverlap(Link[i],1,trasf);

    mm = mymatrix;
    mm *= *matTrunc;
    mm *= tt90;
    mymatrix *= *mat;

    printf("angleX = %f angleY = %f angleZ = %f\n",angleX,angleY,angleZ);

    tt[0][3] = 0;
    tt[1][3] = 0;
    tt[2][3] = -mLength[i]*0.5;

    mm *= tt;

    getAngles(mm, angleX, angleY, angleZ);

    TGeoRotation rot2(Form("rot%d",i+1),0,0,0);
    // angleZ = TMath::ATan2(mm[1][0],mm[0][0])*TMath::RadToDeg();
    // angleX = TMath::ATan2(mm[1][2],mm[2][2])*TMath::RadToDeg();
    rot2.RotateX(angleX);
    rot2.RotateY(angleY);
    rot2.RotateZ(angleZ);
    //rot2.RotateZ(angleZ);
    //rot2.RotateX(-angleX);
    TGeoTranslation tra2(mm[0][3],mm[1][3],mm[2][3]);

    mm.Print();

    Link2[i] = Geo::manager()->MakeTube(Form("LinkB%d",i+1),Geo::iron(),0,0.4,mLength[i]*0.5);
    Link2[i]->SetLineColor(colors[i]);
    Link2[i]->SetFillColor(colors[i]);
    trasf = new TGeoCombiTrans(tra2, rot2);    
    Geo::top()->AddNodeOverlap(Link2[i],1,trasf);
  }
  
  changeGeometry();
}
//_________________________________
void Robot6A::setDeltaMatrix(){
  // define Delta Matrix
    
  for(int i=0; i< mNaxis; i++)
    *(mMatriceDelta[i]) = *mMatriceTras;

  for(int i=0; i< mNaxis; i++){
    TMatrixT<double> mDelta(4,4);

    double stheta = sin(mTheta[i]*TMath::DegToRad());
    double ctheta = cos(mTheta[i]*TMath::DegToRad());
    double stw = sin(mTwist[i]*TMath::DegToRad());
    double ctw = cos(mTwist[i]*TMath::DegToRad());
    mDelta[0][0] = -stheta;
    mDelta[0][1] = -ctheta*ctw;
    mDelta[0][2] =  ctheta*stw;
    mDelta[0][3] = -stheta*mLength[i];
    mDelta[1][0] =  ctheta;
    mDelta[1][1] = -stheta*ctw;
    mDelta[1][2] =  stheta*stw;
    mDelta[1][3] =  ctheta*mLength[i];

    for(int j=0; j< mNaxis; j++){
      if(i != j)
	*(mMatriceDelta[j]) *= *(mLink[i].getMatrix());
      else
	*(mMatriceDelta[j]) *= mDelta;
    }
  }
}
//_________________________________
void Robot6A::changeGeometry(){
  // define Delta Matrix
  TMatrixT<double> mId(4,4);
  mId[0][0] = 1;
  mId[1][1] = 1;
  mId[2][2] = 1;
  mId[3][3] = 1;

  setDeltaMatrix();
  
  TMatrixT<double> mymatrix(4,4);
  mymatrix = *mMatriceTras;

  // links
  int colors[mNaxis];
  colors[0]=2;
  for(int i=6; i < mNaxis;i++)
    colors[i] = colors[i-1] + 2;

  TGeoVolume *Link[mNaxis];
  TGeoVolume *Link2[mNaxis];

  TMatrixT<double> tt(4,4);
  tt[0][0]=1;
  tt[1][1]=1;
  tt[2][2]=1;
  tt[3][3]=1;

  TMatrixT<double> tt90(4,4);
  tt90[0][0]=0;
  tt90[0][2]=1;
  tt90[1][1]=1;
  tt90[2][0]=-1;
  tt90[2][2]=0;
  tt90[3][3]=1;

  for(int i=0; i < mNaxis; i++){
    TMatrixT<double> *mat = mLink[i].getMatrix();
    TMatrixT<double> *matTrunc = mLink[i].getMatrixTrunc();

    TMatrixT<double> mm = mymatrix;

    tt[0][3] = 0;//-mDistance[i]*0.5*sin(angleX*TMath::DegToRad());
    tt[1][3] = 0;//-mDistance[i]*0.5*sin(angleX*TMath::DegToRad())*0;
    tt[2][3] = mDistance[i]*0.5;
    // tt[0][3] = -mDistance[i]*0.5*sin(angleX*TMath::DegToRad());
    // tt[1][3] = -mDistance[i]*0.5*sin(angleX*TMath::DegToRad())*0;
    // tt[2][3] = abs(mDistance[i]*0.5*cos(angleX*TMath::DegToRad()));

    mm *= tt;

    TGeoNode *nn = Geo::top()->GetNode(i*2 + mFirstNodeInTopVolume);
    TGeoMatrix *gm = nn->GetMatrix();
    gm->SetDx(mm[0][3]);
    gm->SetDy(mm[1][3]);
    gm->SetDz(mm[2][3]);
    double *mtochange = const_cast<double*> (gm->GetRotationMatrix());
    mtochange[0] = mm[0][0];
    mtochange[1] = mm[0][1];
    mtochange[2] = mm[0][2];
    mtochange[3] = mm[1][0];
    mtochange[4] = mm[1][1];
    mtochange[5] = mm[1][2];
    mtochange[6] = mm[2][0];
    mtochange[7] = mm[2][1];
    mtochange[8] = mm[2][2];
    
    mm = mymatrix;
    mm *= *matTrunc;
    mm *= tt90;
    mymatrix *= *mat;

    tt[0][3] = 0;
    tt[1][3] = 0;
    tt[2][3] = -mLength[i]*0.5;

    mm *= tt;

    nn = Geo::top()->GetNode(i*2+1 + mFirstNodeInTopVolume);
    gm = nn->GetMatrix();
    gm->SetDx(mm[0][3]);
    gm->SetDy(mm[1][3]);
    gm->SetDz(mm[2][3]);
    mtochange = const_cast<double*> (gm->GetRotationMatrix());
    mtochange[0] = mm[0][0];
    mtochange[1] = mm[0][1];
    mtochange[2] = mm[0][2];
    mtochange[3] = mm[1][0];
    mtochange[4] = mm[1][1];
    mtochange[5] = mm[1][2];
    mtochange[6] = mm[2][0];
    mtochange[7] = mm[2][1];
    mtochange[8] = mm[2][2];
  }

  mCurrentPos[0] = mymatrix[0][3];
  mCurrentPos[1] = mymatrix[1][3];
  mCurrentPos[2] = mymatrix[2][3];

  printf("Current position = %f %f %f\n",mCurrentPos[0],mCurrentPos[1],mCurrentPos[2]);
}
//_________________________________
void Robot6A::draw(){
  Geo::top()->Draw("ogl");
}
//_________________________________
void Robot6A::write(){
  Geo::manager()->Export("geometry.root");
}
//_________________________________
void Robot6A::getAngles(TMatrixT<double> matrix,double &angleX,double &angleY,double &angleZ){
  double sy = sqrt(matrix[0][0]*matrix[0][0] + matrix[1][0]*matrix[1][0]);

  angleY = TMath::ATan2(-matrix[2][0],sy);
  if(sy > 1E-6){
    angleX = TMath::ATan2(matrix[2][1],matrix[2][2]);
    angleZ = TMath::ATan2(matrix[1][0],matrix[0][0]);
  }
  else{
    angleX = TMath::ATan2(-matrix[1][2],matrix[1][1]);
    angleZ = 0;
  }
  angleX *= TMath::RadToDeg();
  angleY *= TMath::RadToDeg();
  angleZ *= TMath::RadToDeg();
}
//_________________________________
void Robot6A::print(){
  for(int i=0; i < mNaxis;i++)
    printf("Link %i: a = %f, alpha = %f, d = %f, theta = %f\n",i,mLength[i],mTwist[i],mDistance[i],mTheta[i]);


  //  printf("X = %f, Y= %f, Z=%f\n",x,y,z);
  
  TMatrixT<double> global(4,4);
  global[0][0]=1;
  global[1][1]=1;
  global[2][2]=1;
  global[3][3]=1;

  for(int i=0; i < mNaxis; i++){
    TMatrixT<double> *mat = mLink[i].getMatrix();
    global *= *mat;
    printf("Position (after link %d): X=%f, Y=%f, Z=%f\n",i,global[0][3],global[1][3],global[2][3]);
  }
}
//_________________________________
void Robot6A::testIteration(double x, double y, double z){
  double toPos[3] = {x,y,z};
  double fromPos[3] = {mCurrentPos[0],mCurrentPos[1],mCurrentPos[2]};

  static int previousLink = -1;

  x -= mCurrentPos[0];
  y -= mCurrentPos[1];
  z -= mCurrentPos[2];
  double stepping = sqrt(x*x + y*y + z*z);
  x /= stepping;
  y /= stepping;
  z /= stepping;

  double sp1[mNaxis],comp1[mNaxis],sp1abs[mNaxis], norm;
  for(int i=0; i< mNaxis; i++){
    //mMatriceDelta[i]->Print();
    norm = sqrt((*mMatriceDelta[i])[0][3]*(*mMatriceDelta[i])[0][3] + (*mMatriceDelta[i])[1][3]*(*mMatriceDelta[i])[1][3] + (*mMatriceDelta[i])[2][3]*(*mMatriceDelta[i])[2][3]);
    comp1[i] = (*mMatriceDelta[i])[0][3] * x + (*mMatriceDelta[i])[1][3] * y + (*mMatriceDelta[i])[2][3] * z;
    sp1[i] = comp1[i]/norm;
    sp1abs[i] = abs(sp1[i]);
    comp1[i] = abs(comp1[i]);
    printf("%d) mDelta= %f %f %f -- dir = %f %f %f -> comp = %f, sp = %f\n",i,(*mMatriceDelta[i])[0][3],(*mMatriceDelta[i])[1][3],(*mMatriceDelta[i])[2][3],x,y,z,comp1[i],sp1[i]);
    printf("ord=%d) sp = %f\n",i,sp1[i]);
  }
  
  int order[mNaxis];
  int order2[mNaxis];
  TMath::Sort(mNaxis,sp1abs,order);
  TMath::Sort(mNaxis,comp1,order2);

  // if(order[0] == previousLink) order[0] = order[1];
  previousLink = order[0];

  double dtheta = 1./sp1[order[0]]/norm*stepping*TMath::RadToDeg();
  
  printf("testIteration: link=%d - scalar product = %f - deltatheta = %f\n",order[0],sp1[order[0]],dtheta);

  if(dtheta > 0.1) dtheta = 0.1;

  setParameters(order[0], mLength[order[0]], mTwist[order[0]], mDistance[order[0]],  mTheta[order[0]] + dtheta); 
  changeGeometry();
  printf("To %f %f %f - From %f %f %f (step = %f)\n",toPos[0],toPos[1],toPos[2],fromPos[0],fromPos[1],fromPos[2],stepping);
  printf("Now at: %f %f %f\n",mCurrentPos[0],mCurrentPos[1],mCurrentPos[2]);
}
//_________________________________
bool Robot6A::moveTo(double x, double y, double z, double stepping){
  double toPos[3] = {x,y,z};
  double fromPos[3] = {mCurrentPos[0],mCurrentPos[1],mCurrentPos[2]};
  x -= mCurrentPos[0];
  y -= mCurrentPos[1];
  z -= mCurrentPos[2];
  double modinv = sqrt(x*x + y*y + z*z);

  //  if(modinv > stepping) modinv = stepping;

  printf("distance = %f (step  = %f)\n",modinv,stepping);

  if(modinv < stepping*2) return 1; // don't move  

  modinv = 1./modinv;
  x *= modinv;
  y *= modinv;
  z *= modinv;

  printf("versor1 %f %f %f\n",x,y,z);

  double v1[3];
  double v2[3];

  if(abs(x) > abs(y) && abs(x) > abs(z)){
    if(abs(y) > 0 || abs(z) > 0){
      v1[0] = 0;
      v1[1] = -z;
      v1[2] = y;
      v2[0] = -(y*y + z*z)/x;
      v2[1] = y;
      v2[2] = z;
    }
    else{
      v1[0] = 0;
      v1[1] = 1;
      v1[2] = 0;
      v2[0] = 0;
      v2[1] = 0;
      v2[2] = 1;
    }
  }
  else if(abs(y) > abs(x) && abs(y) > abs(z)){
    if(abs(x) > 0 || abs(z) > 0){
      v1[0] = -z;
      v1[1] = 0;
      v1[2] = x;
      v2[0] = x;
      v2[1] = -(x*x + z*z)/y;
      v2[2] = z;
    }
    else{
      v1[0] = 1;
      v1[1] = 0;
      v1[2] = 0;
      v2[0] = 0;
      v2[1] = 0;
      v2[2] = 1;
    }
  }
  else{
    if(abs(x) > 0 || abs(y) > 0){
      v1[0] = -y;
      v1[1] = x;
      v1[2] = 0;
      v2[0] = x;
      v2[1] = y;
      v2[2] = -(x*x + y*y)/z;
    }
    else{
      v1[0] = 1;
      v1[1] = 0;
      v1[2] = 0;
      v2[0] = 0;
      v2[1] = 1;
      v2[2] = 0;
    }
  }

  double norm = 1./sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  v1[0] *= norm; 
  v1[1] *= norm; 
  v1[2] *= norm; 
  norm = 1./sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
  v2[0] *= norm; 
  v2[1] *= norm; 
  v2[2] *= norm; 

  printf("versor2 %f %f %f\n",v1[0],v1[1],v1[2]);
  printf("versor3 %f %f %f\n",v2[0],v2[1],v2[2]);

  bool exclusion[mNaxis];
  for(int i=0; i< mNaxis;i++) exclusion[i] = false;

  double sp1[mNaxis],sp2[mNaxis],sp3[mNaxis],sp1abs[mNaxis];
  double weight[100];
  int order[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};

  for(int i=0; i< mNaxis; i++){
    //mMatriceDelta[i]->Print();
    sp1[i] = (*mMatriceDelta[i])[0][3] * x + (*mMatriceDelta[i])[1][3] * y + (*mMatriceDelta[i])[2][3] * z;
    sp2[i] = (*mMatriceDelta[i])[0][3] * v1[0] + (*mMatriceDelta[i])[1][3] * v1[1] + (*mMatriceDelta[i])[2][3] * v1[2];
    sp3[i] = (*mMatriceDelta[i])[0][3] * v2[0] + (*mMatriceDelta[i])[1][3] * v2[1] + (*mMatriceDelta[i])[2][3] * v2[2];
    
    sp1abs[i] = abs(sp1[i]/sqrt(sp1[i]*sp1[i] + sp2[i]*sp2[i] + sp3[i]*sp3[i]));

    printf("%d) %f %f %f\n",i,(*mMatriceDelta[i])[0][3],(*mMatriceDelta[i])[1][3],(*mMatriceDelta[i])[2][3]);
    printf("%d) proj %f %f %f -> %f\n",i,sp1[i],sp2[i],sp3[i],sp1abs[i]);
  }

  bool isgood = false;

  double dtheta[100],maxdtheta;
  double normincr;

  int mode = 1;

  if(mode == 0){
    while(! isgood){
      int n = getAngleToBeMoved(sp1abs,order,exclusion);
      
      getWeight(sp1,sp2,sp3,order,weight);
      
      normincr = stepping/(weight[0]*sp1[order[0]] + weight[1]*sp1[order[1]] + weight[2]*sp1[order[2]]);
      
      maxdtheta=0;
      for(int i=0; i < 3; i++){
	printf("%d) order = %d, delta-theta = %f \n",i,order[i],dtheta[i]);
	dtheta[i] = normincr * weight[i] * TMath::RadToDeg();
	if(abs(dtheta[i]) > maxdtheta) maxdtheta = abs(dtheta[i]);
      }

      if(maxdtheta < 1 || n < 4) isgood = 1;
      else{
	for(int i=2; i >= 0; i--){
	  if(abs(dtheta[i]) > 1){
	    exclusion[order[i]] = true;
	    i = -1;
	  }
	}
      }
    }
  }
  else if(mode == 1){
   getBestWeight(sp1,sp2,sp3,weight);
   
   maxdtheta = 0;
   for(int i=0; i < mNaxis; i++){
     dtheta[i] = stepping * weight[i] * TMath::RadToDeg();
     if(abs(dtheta[i]) > maxdtheta) maxdtheta = abs(dtheta[i]);
   }
  }

  if(maxdtheta > 1){
    printf("MAXDTHETA = %f\n",maxdtheta);
    for(int i=0; i < mNaxis; i++){
     dtheta[i] *= 1/maxdtheta;
    }
  }


  printf("weights: %f %f %f (norm = %f)\n",weight[0],weight[1],weight[2],normincr);

  if (mode == 0){
    for(int i=0; i < 3; i++){
      mTheta[order[i]] += dtheta[i];
      setParameters(order[i], mLength[order[i]], mTwist[order[i]], mDistance[order[i]],  mTheta[order[i]]); 
    }
  }
  else if(mode == 1){
    for(int i=0; i < mNaxis; i++){
      //mTheta[i] += dtheta[i];
      printf("theta final: %d -> %f (increment %f)\n",i,mTheta[i], dtheta[i]);
    
      //setParameters(i, mLength[i], mTwist[i], mDistance[i],  mTheta[i]); 
    }
  }

  double realmov[3] = {0., 0., 0.};
  if(mode == 1){
    for(int i=0; i < mNaxis;i++){
      realmov[0] += (*mMatriceDelta[i])[0][3]*dtheta[i]*TMath::DegToRad();
      realmov[1] += (*mMatriceDelta[i])[1][3]*dtheta[i]*TMath::DegToRad();
      realmov[2] += (*mMatriceDelta[i])[2][3]*dtheta[i]*TMath::DegToRad();
    }
  }

  //  changeGeometry();
  checkDerivative(dtheta);
  printf("To %f %f %f - From %f %f %f (step = %f)\n",toPos[0],toPos[2],toPos[2],fromPos[0],fromPos[1],fromPos[2],stepping);
  printf("Now at: %f %f %f\n",mCurrentPos[0],mCurrentPos[1],mCurrentPos[2]);
  printf("calcolated moved of: %e %e %e\n",realmov[0],realmov[1],realmov[2]);
  printf("moved of: %e %e %e\n",mCurrentPos[0]-fromPos[0],mCurrentPos[1]-fromPos[1],mCurrentPos[2]-fromPos[2]);
  printf("expected: %e %e %e\n",x*stepping,y*stepping,z*stepping);
  return 0;
}
//_________________________________
int Robot6A::getAngleToBeMoved(double *sp1absOr, int *order, bool *exclusion){
  double sp1abs[mNaxis];
  int order1[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  int order2[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  int order3[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  int select[3] = {0,1,2};

  int neff = 0;
  for(int i=0; i < mNaxis; i++){
    order[neff] = i;
    sp1abs[neff] = sp1absOr[i];
    if(! exclusion[i]) neff++;
  }

  // first link to be moved -> the one with the movement closer to direction we want to follow
  TMath::Sort(neff,sp1abs,order1);
  select[0] = order[order1[0]];

  double orth1[mNaxis];
  for(int j=0; j< neff; j++){
    int i = order[j];
    double norm = sqrt((*mMatriceDelta[i])[0][3]*(*mMatriceDelta[i])[0][3] + (*mMatriceDelta[i])[1][3]*(*mMatriceDelta[i])[1][3] + (*mMatriceDelta[i])[2][3]*(*mMatriceDelta[i])[2][3]);
    orth1[j] = abs((*mMatriceDelta[select[0]])[0][3]*(*mMatriceDelta[i])[0][3] + (*mMatriceDelta[select[0]])[1][3]*(*mMatriceDelta[i])[1][3] + (*mMatriceDelta[select[0]])[2][3]*(*mMatriceDelta[i])[2][3]);
    orth1[j] /= norm;
    printf("orth1: %d -> %f\n",i,orth1[j]);
  }

  // take as second link the one most orthogonal to the first one
  TMath::Sort(neff,orth1,order2);
  select[1] = order[order2[neff-1]];
  for(int j=0; j< neff; j++){
    printf("order2(%d) = %d\n",j,order2[j]);
  }


  // take the direction orthogonal to both (vectorial product)
  double dir3[3];
  dir3[0] = (*mMatriceDelta[select[0]])[1][3]*(*mMatriceDelta[select[1]])[2][3] - (*mMatriceDelta[select[0]])[2][3]*(*mMatriceDelta[select[1]])[1][3];
  dir3[1] = (*mMatriceDelta[select[0]])[2][3]*(*mMatriceDelta[select[1]])[0][3] - (*mMatriceDelta[select[0]])[0][3]*(*mMatriceDelta[select[1]])[2][3];
  dir3[2] = (*mMatriceDelta[select[0]])[0][3]*(*mMatriceDelta[select[1]])[1][3] - (*mMatriceDelta[select[0]])[1][3]*(*mMatriceDelta[select[1]])[0][3];

  double orth2[mNaxis];

  for(int j=0; j< neff; j++){
    int i = order[j];
    double norm = sqrt((*mMatriceDelta[i])[0][3]*(*mMatriceDelta[i])[0][3] + (*mMatriceDelta[i])[1][3]*(*mMatriceDelta[i])[1][3] + (*mMatriceDelta[i])[2][3]*(*mMatriceDelta[i])[2][3]);
    orth2[j] = abs(dir3[0]*(*mMatriceDelta[i])[0][3] + dir3[1]*(*mMatriceDelta[i])[1][3] + dir3[2]*(*mMatriceDelta[i])[2][3]);
    orth2[j] /= norm;
    printf("orth2: %d -> %f\n",i,orth2[j]);
  }

  // take as third link the one most orthogonal to both
  TMath::Sort(neff,orth2,order3);
  select[2] = order[order3[0]];

  printf("dir3: %f %f %f\n",dir3[0],dir3[1],dir3[2]);
  printf("order: %d %d %d\n",order[0],order[1],order[2]);
  printf("select: %d %d %d\n",select[0],select[1],select[2]);

  order[0] = select[0];
  order[1] = select[1];
  order[2] = select[2];

  return neff;
}
//_________________________________
void Robot6A::getWeight(double *sp1,double *sp2,double *sp3,int *order, double weight[3]){
  if(sp1[order[0]] > 0) weight[0] = 1;
  else weight[0] = -1;

  weight[1] = (sp3[order[0]] * sp2[order[2]] / sp3[order[2]] - sp2[order[0]]) / (sp2[order[1]] - sp3[order[1]] * sp2[order[2]] / sp3[order[2]]) * weight[0];
  weight[2] = (sp3[order[0]] * sp2[order[1]] / sp3[order[1]] - sp2[order[0]]) / (sp2[order[2]] - sp3[order[2]] * sp2[order[1]] / sp3[order[1]]) * weight[0];
}
//_________________________________
void Robot6A::getBestWeight(double *sp1,double *sp2,double *sp3, double *weight){
  int nv = mNaxis;

  double vvx[100],vvy[100],vvz[100],vvxabs[100],vvyabs[100],vvzabs[100];
  double ww1[100],ww2[100],stepEff=0;
  double *ww3 = weight;
  int order[100];
  float sign[100];
  double vvxOr[100],vvyOr[100],vvzOr[100];

  for(int i=0; i < nv; i++){
    vvxOr[i] = sp1[i];
    vvyOr[i] = sp2[i];
    vvzOr[i] = sp3[i];
    vvx[i] = vvxOr[i];
    vvy[i] = vvyOr[i];
    vvz[i] = vvzOr[i];
    vvxabs[i] = abs(vvxOr[i]);
    vvyabs[i] = abs(vvyOr[i]);
    vvzabs[i] = abs(vvzOr[i]);
  }

  // step 1
  TMath::Sort(nv,vvzabs,order);
  int izsub = order[0];

  for(int i=0; i < nv; i++){
    if(i!= izsub) ww1[i] = vvz[i]/vvz[izsub];
    else ww1[i] = 0;
    vvx[i] -= ww1[i] * vvx[izsub];
    vvy[i] -= ww1[i] * vvy[izsub];
    vvz[i] -= ww1[i] * vvz[izsub];
    vvyabs[i] = abs(vvy[i]);
  }
  vvx[izsub] = 0;
  vvy[izsub] = 0;
  vvyabs[izsub] = 0;
  vvz[izsub] = 0;

  // step 2
  TMath::Sort(nv,vvyabs,order);
  int iysub = order[0];

  for(int i=0; i < nv; i++){
    if(i!= iysub) ww2[i] = vvy[i]/vvy[iysub];
    else ww2[i] = 0;
    vvx[i] -= ww2[i] * vvx[iysub];
    vvy[i] -= ww2[i] * vvy[iysub];
    vvz[i] -= ww2[i] * vvz[iysub];
    ww3[i] = vvx[i] * vvx[i];

    if(vvx[i] > 0) sign[i] = 1;
    else sign[i] = -1;

    if(i!= iysub) stepEff += ww3[i] * vvx[i] * sign[i];
    printf("%f %f %f\n",vvx[i],vvy[i],vvz[i]);
  }
  vvx[iysub] = 0;
  vvxabs[izsub] = 0;
  vvy[iysub] = 0;
  vvz[iysub] = 0;
  ww3[iysub] = 0;
  
  // step 3
  for(int i=0; i < nv; i++){
    if(i != izsub && i != iysub){
      ww3[i] /= stepEff;
      ww3[i] *= sign[i];
      ww3[izsub] += ww3[i]*(ww1[iysub]*ww2[i] - ww1[i]);
      ww3[iysub] += -ww3[i] * ww2[i];
    }
  }

  float check[3] = {0,0,0};
  for(int i=0; i < nv; i++){
    check[0] += vvxOr[i]*ww3[i];
    check[1] += vvyOr[i]*ww3[i];
    check[2] += vvzOr[i]*ww3[i];
  }

  printf("check: %f %f %f\n",check[0],check[1],check[2]);
 
}
//_________________________________
void Robot6A::checkDerivative(double *dtheta){
  double x0 = mCurrentPos[0];
  double y0 = mCurrentPos[1];
  double z0 = mCurrentPos[2];

  double pos[100][3];
  double der[100][3];

  for(int i=0; i <  mNaxis; i++){
    der[i][0] = (*mMatriceDelta[i])[0][3]*dtheta[i]*TMath::DegToRad();
    der[i][1] = (*mMatriceDelta[i])[1][3]*dtheta[i]*TMath::DegToRad();
    der[i][2] = (*mMatriceDelta[i])[2][3]*dtheta[i]*TMath::DegToRad();
  }

  for(int i=0; i <  mNaxis; i++){
    x0 = mCurrentPos[0];
    y0 = mCurrentPos[1];
    z0 = mCurrentPos[2];
    setParameters(i, mLength[i], mTwist[i], mDistance[i],  mTheta[i] + dtheta[i]);
    changeGeometry();
    pos[i][0] = (mCurrentPos[0] - x0);
    pos[i][1] = (mCurrentPos[1] - y0);
    pos[i][2] = (mCurrentPos[2] - z0);
    printf("%d) %f %f %f  -  %f %f %f\n",i,pos[i][0],pos[i][1],pos[i][2],der[i][0],der[i][1],der[i][2]);
  }

}
