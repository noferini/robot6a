#ifndef CLASSI_DEF
#define CLASSI_DEF

#include "TMatrix.h"
#include "TGeoManager.h"
#include "TMultiLayerPerceptron.h"
#include "TMath.h"

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
  void setDeltaMatrix();
  void testIteration(double x, double y, double z);
  void draw();
  void write();
  void setParameters(int ilink, double length, double twist, double distance, double theta);
  void setTheta(int ilink, double theta) {mTheta[ilink]=theta; setParameters(ilink,mLength[ilink],mTwist[ilink],mDistance[ilink],mTheta[ilink]);}
  void getAngles(TMatrixT<double> matrix,double &angleX,double &angleY,double &angleZ);
  bool moveTo(double x, double y, double z, double stepping);
  bool moveToNN(double x, double y, double z, double stepping);
  int getAngleToBeMoved(double *sp1absOr, int *order, bool *exclusion);
  void getWeight(double *sp1,double *sp2,double *sp3,int *order, double weight[3]);
  void getBestWeight(double *sp1,double *sp2,double *sp3, double *weight);
  void setOrigin(double x, double y, double z){(*mMatriceTras)[0][3]=x;(*mMatriceTras)[1][3]=y;(*mMatriceTras)[2][3]=z;};
  void rotate(int iaxis,double alpha);
  void rotateX(double alpha) {rotate(0, alpha);}
  void rotateY(double alpha) {rotate(1, alpha);}
  void rotateZ(double alpha) {rotate(2, alpha);}
  void print();

  double getDeltaMatrix(int ilink, int coord);
  
  double getTheta(int i) {return mTheta[i];}

  int getNaxis() const {return mNaxis;}

  void checkDerivative(double *dtheta);
  void loadMlp(const char* filename);

  double getX() const {return mCurrentPos[0];}
  double getY() const {return mCurrentPos[1];}
  double getZ() const {return mCurrentPos[2];}

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

  TMultiLayerPerceptron *mMlp;
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

class NN{
 public:
  NN();
  ~NN();
  static double sigmoid(double x) {return 1./(1 + TMath::Exp(-x));}
  void setNinput(int n);
  void setNoutput(int n);
  void setNlayer(int n);
  void setNhidden(int i, int n);
  void createNN();
  void print();
  double* eval(double *in);
  double* testEvalW(double *in,int layer, int link, double ww=1.);
  double* testEvalN(double *in,int layer, int link, double ww=1.);
  void updateNNw(int layer, int link,double ww, double w0=mDeltaTrain);
  void updateNNn(int layer, int node,double ww, double w0=mDeltaTrain);
  void loadWeights(const char *file);
  void dumpWeights(const char *fout="");
  void draw();
  void randomize();
  int getNlayer() const {return mNlayer;}
  int getNweight(int layer) const;
  int getNnodes(int i) const {return mNhidden[i];}
 private:
  int mNinput = 0;
  int mNoutput = 0;
  int mNlayer = 0;
  int mReady = 0;
  int *mNhidden;
  double **mWeights=0;
  double **mDeltaWeights=0;
  double **mEvalNode=0;
  double *mInputOffset=0;
  double *mInputScale=0;
  double *mInputWeight=0;
  double *mOutputOffset=0;
  double *mOutputScale=0;
  double *mOutputWeight=0;
  double **mNodeWeight=0;

  int mLayTrain;
  int mWeightTrain;
  static double mDeltaTrain;
};
#endif
