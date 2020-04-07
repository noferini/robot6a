#ifndef CLASSI_DEF
#define CLASSI_DEF

#include "TMatrix.h"
#include "TGeoManager.h"

class Link{
 public:
  Link();
  ~Link();
  void print() const {mMatrice->Print(); mMatriceInv->Print();}
  void setParameters(double length, double twist, double distance, double theta); // a_i, alpha_i, d_i, theta_i
  void checkInv();
  TMatrixT<double> *getMatrix() {return mMatrice;}
  TMatrixT<double> *getMatrixTrunc() {return mMatriceTrunc;}

 private:
  TMatrixT<double> *mMatrice;
  TMatrixT<double> *mMatriceTrunc;
  TMatrixT<double> *mMatriceInv;
};

class Robot6A{
 public:
  Robot6A(int n=6);
  ~Robot6A();
  void setGeometry();
  void changeGeometry();
  void draw();
  void write();
  void setParameters(int ilink, double length, double twist, double distance, double theta);
  void getAngles(TMatrixT<double> matrix,double &angleX,double &angleY,double &angleZ);
  bool moveTo(double x, double y, double z, double stepping);
  void setOrigin(double x, double y, double z){(*mMatriceTras)[0][3]=x;(*mMatriceTras)[1][3]=y;(*mMatriceTras)[2][3]=z;};
  void rotate(int iaxis,double alpha);
  void rotateX(double alpha) {rotate(0, alpha);}
  void rotateY(double alpha) {rotate(1, alpha);}
  void rotateZ(double alpha) {rotate(2, alpha);}
  void print();

  int getNaxis() const {return mNaxis;}
 
 private:
  const int mNaxis = 6;
  int mFirstNodeInTopVolume;
  TMatrixT<double> *mMatriceTras;
  Link *mLink;
  float *mLength;
  float *mTwist;
  float *mDistance;
  float *mTheta;

  double mCurrentPos[3];

  TMatrixT<double> **mMatriceDelta;
};

class Geo{
 public:
  static TGeoManager* manager() {return mManager;}
  static TGeoMedium* air() {return mAir;}
  static TGeoMedium* iron() {return mIron;}
  static void createTopVolume(float xsize, float ysize, float zsize) { mTop = mManager->MakeBox("top",Geo::air(),xsize,ysize,zsize); mManager->SetTopVolume(mTop); mTop->SetVisibility(0);}
  static TGeoVolume* top();
  static void closeGeo() {mManager->CloseGeometry();}
 private:
  static TGeoManager *mManager;
  static TGeoMaterial *mVacuum;  
  static TGeoMaterial *mFe; 
  static TGeoVolume *mTop;
  static TGeoMedium *mAir;
  static TGeoMedium *mIron;

};
#endif
