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
  Robot6A();
  ~Robot6A();
  void setGeometry();
  void changeGeometry();
  void draw();
  void write();
  void setParameters(int ilink, double length, double twist, double distance, double theta);
  void getAngles(TMatrixT<double> matrix,double &angleX,double &angleY,double &angleZ);
  bool moveTo(double x, double y, double z, double stepping);

  void print();

 private:
  Link mLink[6];
  float mLength[6];
  float mTwist[6];
  float mDistance[6];
  float mTheta[6];

  double mCurrentPos[3];

  TMatrixT<double> *mMatriceDelta[6];
 
  TGeoManager *mManager;
  TGeoVolume *mTop;
  TGeoMedium *mAir;
  TGeoMedium *mIron;
};

#endif
