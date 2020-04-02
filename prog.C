#include <stdlib.h>
#include "TApplication.h"
#include "TRint.h"

#include "classi.h"
#include "TCanvas.h"
#include "TGLViewer.h"
#include "TGLOrthoCamera.h"
#include "TTimer.h"
#include "TRandom.h"
#include "TVirtualPad.h"
#include "TGLLightSet.h"
#include "TGLSAViewer.h"
#include "TGFrame.h"

using namespace std;

Robot6A rob;
TGLViewer::ECameraType camera;
TTimer timer(25);

int main(int argc, char **argv){
  TRint *theApp = new TRint("eeeroot", &argc, argv);

//  Robot6A rob;

  double a1=1,a2=4,a3=8,a4=1,a5=1,a6=1;
  double d1=1,d2=4, d3=1, d4=2, d5=1, d6=3;
  double th1=30,th2=30,th3=0,th4=20,th5=14,th6=40;

  rob.setParameters(0, a1, 90, d1, th1);
  rob.setParameters(1, a2, 0, d2, th2);
  rob.setParameters(2, a3, 90, d3, th3);
  rob.setParameters(3, a4, -90, d4, th4);
  rob.setParameters(4, a5, 0, d5, th5);
  rob.setParameters(5, a6, 0, d6, th6);

  rob.setGeometry();
  rob.changeGeometry();

  rob.setParameters(3, a4, -90, d4, th4+10);

  rob.draw();

  TGLViewer * v = (TGLViewer *)gPad->GetViewer3D();

  // Random draw style
  Int_t style = 0;
  switch (style)
  {
    case 0: v->SetStyle(TGLRnrCtx::kFill); break;
    case 1: v->SetStyle(TGLRnrCtx::kOutline); break;
    case 2: v->SetStyle(TGLRnrCtx::kWireFrame); break;
  }

  TGLLightSet* ls = v->GetLightSet();
  if (0)
      ls->SetLight(TGLLightSet::kLightLeft, kFALSE);
  if (0)
      ls->SetLight(TGLLightSet::kLightRight, kFALSE);
  if (0)
      ls->SetLight(TGLLightSet::kLightTop, kFALSE);
  if (0)
      ls->SetLight(TGLLightSet::kLightBottom, kFALSE);

  int id = 1;
  camera = (TGLViewer::ECameraType) 1; // 0-6
//  v->SetCurrentCamera(camera);
//  v->CurrentCamera().SetExternalCenter(kTRUE);
  if (id > 2) {
      //0, 1, and 2 - are different 'perspective' camers.
      TGLOrthoCamera& o = static_cast<TGLOrthoCamera &>(v->CurrentCamera());
      o.SetEnableRotate(kTRUE);
  }

  // Now animate the camera
  TGLSAViewer* sav = dynamic_cast<TGLSAViewer*>(v);
  if (sav)
     sav->GetFrame()->Connect("CloseWindow()", "TTimer", &timer, "TurnOff()");
//  timer.SetCommand("AnimateCamera()");
  timer.TurnOn();

  rob.changeGeometry();

  rob.print();

  rob.write();

  theApp->Run();  

  return 0;
}
