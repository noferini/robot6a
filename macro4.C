using namespace std;

Robot6A rob;
TGLViewer::ECameraType camera;
TTimer timer(5);

double a1=1,a2=4,a3=8,a4=1,a5=1,a6=1;
double d1=1,d2=4, d3=1, d4=2, d5=1, d6=3;
double th1=30+180,th2=30,th3=0,th4=20,th5=14,th6=40;

int increm = -100;

double pos[5][3] = {{-15,-1,4},{-14,-2,5},{-13,-1,5},{-14,0,4},{-15,-1,4}};

void AnimateCamera(){
  TGLViewer * v = (TGLViewer *)gPad->GetViewer3D();

  if(increm < 0){
    increm++;
  }
  else{
    for(int jj=0;jj < 20;jj++){
      if(increm < 5){
	if(rob.moveTo(pos[increm][0],pos[increm][1],pos[increm][2],0.001))
	  increm++;
      }
    }
  }

  v->UpdateScene();
  
}

void macro4(){

//  Robot6A rob;
  rob.setParameters(0, a1, 90, d1, th1);
  rob.setParameters(1, a2, 0, d2, th2);
  rob.setParameters(2, a3, 90, d3, th3);
  rob.setParameters(3, a4, -90, d4, th4);
  rob.setParameters(4, a5, 90, d5, th5);
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
  if (1)
      ls->SetLight(TGLLightSet::kLightLeft, kFALSE);
  if (0)
      ls->SetLight(TGLLightSet::kLightRight, kFALSE);
  if (0)
      ls->SetLight(TGLLightSet::kLightTop, kFALSE);
  if (0)
      ls->SetLight(TGLLightSet::kLightBottom, kFALSE);

  int id = 1;
  camera = (TGLViewer::ECameraType) 2; // 0-6
  v->SetCurrentCamera(camera);
  v->CurrentCamera().SetExternalCenter(kTRUE);
  double box3[8][3] = {{-30,-30,-30},{-30,-30,30},{-30,30,30},{-30,30,-30},{30,30,-30},{30,-30,-30},{30,-30,30},{30,30,30}};
  v->CurrentCamera().Setup(box3);
  double c[3] = {0,0,5};
  v->CurrentCamera().Configure(0.001,1,c,0,0);

  // Now animate the camera
  TGLSAViewer* sav = dynamic_cast<TGLSAViewer*>(v);
  if (sav)
     sav->GetFrame()->Connect("CloseWindow()", "TTimer", &timer, "TurnOff()");
  timer.SetCommand("AnimateCamera()");
  timer.TurnOn();

  rob.print();

  rob.write();

}
