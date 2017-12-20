////////////////////////////////////////////////////////////////////////
//
//   Harvard University
//   CS175 : Computer Graphics
//   Professor Steven Gortler
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <list>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#define GL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED 1
#ifdef __MAC__
#   include <OpenGL/gl3.h>
#   include <GLUT/glut.h>
#else
#   include <GL/glew.h>
#   include <GL/glut.h>
#endif

#include "cvec.h"
#include "matrix4.h"
#include "geometrymaker.h"
#include "ppm.h"
#include "glsupport.h"
#include "quat.h"
#include "rigtform.h"
#include "arcball.h"
#include "asstcommon.h"
#include "scenegraph.h"
#include "drawer.h"
#include "picker.h"
#include "sgutils.h"
#include "geometry.h"
#include "mesh.h"
#include "hermitedata.hpp"

using namespace std;

#pragma mark - Globals
// G L O B A L S ///////////////////////////////////////////////////
const bool g_Gl2Compatible = false;

static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static const float g_frustNear = -0.1;    // near plane
static const float g_frustFar = -50.0;    // far plane
static const float g_groundSize = 10.0;   // half the ground length

static int g_windowWidth = 800;
static int g_windowHeight = 600;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static bool g_spaceDown = false;         // space state, for middle mouse emulation
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event

// --------- Materials
// This should replace all the contents in the Shaders section, e.g., g_numShaders, g_shaderFiles, and so on
static shared_ptr<Material> g_redDiffuseMat,
g_blueDiffuseMat,
g_bumpFloorMat,
g_arcballMat,
g_pickingMat,
g_lightMat;

shared_ptr<Material> g_overridingMaterial;


// --------- Geometry
typedef SgGeometryShapeNode MyShapeNode;

// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_sphere;

// --------- Scene

static shared_ptr<SgRootNode> g_world;
static shared_ptr<SgRbtNode> g_skyNode, g_groundNode,
g_robot1Node, g_robot2Node, g_light1Node, g_light2Node;
static shared_ptr<SgRbtNode> g_currentPickedRbtNode; // used later when you do picking

static weak_ptr<SgRbtNode> g_eyeNode;
static int g_eyeMode = 0;
static bool g_egoMode = false;
static double g_arcballScreenRadius = 0.25 * min(g_windowWidth, g_windowHeight);
static double g_arcballScale = 1.0f;
static Cvec2 g_arcBallScreenCoords;
static int g_mouseDownX, g_mouseDownY;
static bool g_pickMode = false;

typedef vector<RigTForm> Frame;

// Should be updated every time new objects are added or subtracted to scene
list<Frame> g_frames;
list<Frame>::iterator g_currentFrameIt;
vector<shared_ptr<SgRbtNode> > g_frameNodes;

static int g_msBetweenKeyFrames = 2000; // 2 seconds between keyframes
static int g_animateFramesPerSecond = 30; // frames to render per second during animation playback
static bool g_isPlaying = false;

// Material

static vector<shared_ptr<Material> > g_blobShellMats; // for blob shells

// New Geometry
static double g_radius = 6;
static double g_thickness = 4;
static int g_model = 0;
static double g_blobX = 12.0;
static double g_blobY = 12.0;
static bool g_wireframe = false;
static bool g_sharp = false;
static int g_debug = 0;
static double g_tolerance = 0.1;

static shared_ptr<SimpleGeometryPNX> g_blobGeometry;
static Mesh g_blobMesh;

// New Scene node
static shared_ptr<SgRbtNode> g_blobNode;


///////////////// END OF G L O B A L S //////////////////////////////////////////////////

#pragma mark - Blobby
// Specifying shell geometries based on g_tipPos, g_radius, and g_numShells.
// You need to call this function whenver the shell needs to be updated
static void updateBlobGeometry() {
    
  vector<VertexPNX> newVerts;
  const double scale = 0.1;
  DC::HermiteData<24, 24, 24> voxels(ToggleImplicit(Cvec3(g_blobX, g_blobY, 12.0), g_radius, g_thickness, g_model));
  if (g_debug == 0) {
    voxels.triangulateToVector(newVerts, scale, g_tolerance, g_sharp);
  } else if (g_debug == 1) {
    voxels.triangulateToVectorDebug(newVerts, scale, false);
  } else {
    voxels.triangulateToVectorDebug(newVerts, scale, true);
  }
  
  //
  g_blobGeometry->upload(&newVerts[0], newVerts.size());
}

#pragma mark - Initiation

// New function that loads the blob mesh and initializes the blob shell meshes
static void initBlobbyMeshes() {
  
  // Now allocate array of SimpleGeometryPNX to for shells, one per layer
  g_blobGeometry.reset(new SimpleGeometryPNX());
}

static void initGround() {
  int ibLen, vbLen;
  getPlaneVbIbLen(vbLen, ibLen);
  
  // Temporary storage for cube Geometry
  vector<VertexPNTBX> vtx(vbLen);
  vector<unsigned short> idx(ibLen);
  
  makePlane(g_groundSize*2, vtx.begin(), idx.begin());
  g_ground.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);
  
  // Temporary storage for cube Geometry
  vector<VertexPNTBX> vtx(vbLen);
  vector<unsigned short> idx(ibLen);
  
  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initSphere() {
  int ibLen, vbLen;
  getSphereVbIbLen(20, 10, vbLen, ibLen);
  
  // Temporary storage for sphere Geometry
  vector<VertexPNTBX> vtx(vbLen);
  vector<unsigned short> idx(ibLen);
  makeSphere(1, 20, 10, vtx.begin(), idx.begin());
  g_sphere.reset(new SimpleIndexedGeometryPNTBX(&vtx[0], &idx[0], vtx.size(), idx.size()));
}

// takes a projection matrix and send to the the shaders
inline void sendProjectionMatrix(Uniforms& uniforms, const Matrix4& projMatrix) {
  uniforms.put("uProjMatrix", projMatrix);
}

// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
static void updateFrustFovY() {
  if (g_windowWidth >= g_windowHeight)
    g_frustFovY = g_frustMinFov;
  else {
    const double RAD_PER_DEG = 0.5 * CS175_PI/180;
    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
           g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}

static void drawStuff(bool picking) {

  // Declare an empty uniforms
  Uniforms uniforms;
  
  // Build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(uniforms, projmat);
  
  const shared_ptr<SgRbtNode> eyeNode = g_eyeNode.lock();
  const RigTForm eyeRbt = getPathAccumRbt(g_world, eyeNode);
  const RigTForm invEyeRbt = inv(eyeRbt);

  // get world space coordinates of the light
  Cvec3 light1 = getPathAccumRbt(g_world, g_light1Node).getTranslation();
  Cvec3 light2 = getPathAccumRbt(g_world, g_light2Node).getTranslation();
  
  uniforms.put("uLight", Cvec3(invEyeRbt * Cvec4(light1, 1)));
  uniforms.put("uLight2", Cvec3(invEyeRbt * Cvec4(light2, 1)));
  
  if (!picking) {
    
    // Prepare blob for rendering
    /*if (g_shellNeedsUpdate) {
      updateBlobGeometry();
      g_shellNeedsUpdate = false;
    }*/
    
    // Draw world
    Drawer drawer(invEyeRbt, uniforms);
    g_world->accept(drawer);
    
    // Draw arcball as part of asst3
    if (g_currentPickedRbtNode) {
      const RigTForm location = getPathAccumRbt(g_world, g_currentPickedRbtNode);
      const Cvec3 sphereEye = Cvec3(invEyeRbt * Cvec4(location.getTranslation(), 1));
      if (sphereEye[2] < 0) {
        g_arcBallScreenCoords = getScreenSpaceCoord(sphereEye, projmat, g_frustNear, g_frustFovY, g_windowWidth, g_windowHeight);
        if (!(g_mouseLClickButton && g_mouseRClickButton) && !g_mouseMClickButton) {
          g_arcballScale = getScreenToEyeScale(sphereEye[2], g_frustFovY, g_windowHeight);
        }
        const double radius = g_arcballScale * g_arcballScreenRadius;
        
        const Matrix4 MVM = rigTFormToMatrix(invEyeRbt * location) * Matrix4::makeScale({radius, radius, radius});
        const Matrix4 NMVM = normalMatrix(MVM);
        
        sendModelViewNormalMatrix(uniforms, MVM, NMVM);
        g_arcballMat->draw(*g_sphere, uniforms);
      }
    }
  }
  else {
    Picker picker(invEyeRbt, uniforms);
    
    // Set overiding material to our picking material
    g_overridingMaterial = g_pickingMat;
    
    g_world->accept(picker);
    
    // Unset the overriding material
    g_overridingMaterial.reset();
    
    glFlush();
    g_currentPickedRbtNode = picker.getRbtNodeAtXY(g_mouseClickX, g_mouseClickY);
    if (g_currentPickedRbtNode == g_groundNode) {
      g_currentPickedRbtNode = shared_ptr<SgRbtNode>();   // set to NULL
    }
  }
}

static void pick() {
  // We need to set the clear color to black, for pick rendering.
  // so let's save the clear color
  GLdouble clearColor[4];
  glGetDoublev(GL_COLOR_CLEAR_VALUE, clearColor);
  
  glClearColor(0, 0, 0, 0);
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  // No more glUseProgram
  drawStuff(true); // no more curSS
  
  // Uncomment below and comment out the glutPostRedisplay in mouse(...) call back
  // to see result of the pick rendering pass
  // glutSwapBuffers();
  
  //Now set back the clear color
  glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);
  
  checkGlErrors();
}

#pragma mark GLUT Callbacks

static void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawStuff(false);
  glutSwapBuffers();
  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  g_arcballScreenRadius = 0.25 * min(g_windowWidth, g_windowHeight);
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  updateFrustFovY();
  glutPostRedisplay();
}

static void motion(const int x, const int y) {
  const double dx = x - g_mouseClickX;
  const double dy = g_windowHeight - y - 1 - g_mouseClickY;
  
  double translationAmount = 0.01;
  if (g_currentPickedRbtNode != nullptr && g_eyeMode == 0) {
    translationAmount = g_arcballScale;
  }
  
  RigTForm m;
  bool isTranslating = false;
  if (g_mouseLClickButton && !g_mouseRClickButton && !g_spaceDown) { // left button down?
    m = RigTForm(Quat::makeXRotation(-dy)) * RigTForm(Quat::makeYRotation(dx));
  }
  else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
    m = RigTForm(Cvec3(dx, dy, 0) * translationAmount);
    isTranslating = true;
  }
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton) || (g_mouseLClickButton && !g_mouseRClickButton && g_spaceDown) ) {  // middle or (left and right, or left + space) button down?
    m = RigTForm(Cvec3(0, 0, -dy) * translationAmount);
    isTranslating = true;
  }

  if (g_mouseClickDown) {
    const shared_ptr<SgRbtNode> manipulationTarget = g_currentPickedRbtNode == nullptr ? g_skyNode : g_currentPickedRbtNode;
    const RigTForm targetRbt = manipulationTarget->getRbt();
    const RigTForm targetRbtGlobal = getPathAccumRbt(g_world, manipulationTarget);
    const RigTForm skyRbtGlobal = getPathAccumRbt(g_world, g_skyNode);
    RigTForm newRbt;
    
    // Modify cube from sky camera
    if (g_currentPickedRbtNode != nullptr && g_eyeMode == 0 && !isTranslating) {
      
      // This is "A"
      const RigTForm refFrameGlobal = transFact(targetRbtGlobal) * linFact(skyRbtGlobal);
      
      const RigTForm parentFrameGlobal = getPathAccumRbt(g_world, manipulationTarget, 1);
      const RigTForm invParentFrameGlobal = inv(parentFrameGlobal);
      
      // This is "A_s"
      const RigTForm refFrame = invParentFrameGlobal * refFrameGlobal;
      //const RigTForm refFrame;
      const RigTForm invRefFrame = inv(refFrame);
      
      // Get position of the mouse on the screen-sphere
      const double r = g_arcballScreenRadius;
      const double glX = x;
      const double glY = g_windowHeight - y - 1;
      double offX = glX - g_arcBallScreenCoords[0];
      double offY = glY - g_arcBallScreenCoords[1];
      double offZ2 = r * r - offX * offX - offY * offY;
      double offZ = offZ2 >= 0 ? sqrt(offZ2) : 0;
      const Cvec3 arcPosition = {offX, offY, offZ};
      const Cvec3 arcPositionNorm = normalize(arcPosition);
      
      offX = g_mouseClickX - g_arcBallScreenCoords[0];
      offY = g_mouseClickY - g_arcBallScreenCoords[1];
      offZ2 = r * r - offX * offX - offY * offY;
      offZ = offZ2 >= 0 ? sqrt(offZ2) : 0;
      /* CLICK AWAY FUNCTIONALITY if (offZ2 < 0) {
        g_currentPickedRbtNode = shared_ptr<SgRbtNode>();
        g_pickMode = false;
        return;
      }*/
      const Cvec3 arcPositionOrig = {offX, offY, offZ};
      const Cvec3 arcPositionOrigNorm = normalize(arcPositionOrig);
      
      // Calculate
      const double phi = acos(dot(arcPositionOrigNorm, arcPositionNorm));
      const Cvec3 cp = cross(arcPositionOrigNorm, arcPositionNorm);
      const double cpNorm = norm(cp);
      if (cpNorm > CS175_EPS) {
        const Cvec3 k = cp / cpNorm;
        const Quat rotation(cos(phi), k * sin(phi));
        const RigTForm trans = RigTForm(rotation);
        
        newRbt = refFrame * trans * invRefFrame * targetRbt;
        manipulationTarget->setRbt(newRbt);
      }
      
    // Standard translation of objects
    }
    /* else */ if (g_currentPickedRbtNode != nullptr && g_eyeMode == 0 && isTranslating) {
      
      // Reference frame is the cube-sky frame
      const RigTForm refFrame = transFact(targetRbtGlobal) * linFact(skyRbtGlobal);
      const RigTForm invRefFrame = inv(refFrame);
      
      newRbt = refFrame * m * invRefFrame * targetRbt;
      manipulationTarget->setRbt(newRbt);
    }
    
    if (g_currentPickedRbtNode == nullptr && g_eyeMode == 0) {
    
      // Orbit
      if (!g_egoMode && !isTranslating) {
        const RigTForm skyRbt = g_skyNode->getRbt();
        const RigTForm refFrameX = linFact(skyRbt);
        const RigTForm invRefFrameX = inv(refFrameX);
        newRbt = refFrameX * RigTForm(Quat::makeXRotation(dy)) * invRefFrameX * targetRbt;
        newRbt = RigTForm(Quat::makeYRotation(-dx)) * newRbt;
        manipulationTarget->setRbt(newRbt);
        
      // Ego
      } else if (!isTranslating) {
        const RigTForm refFrameY = transFact(skyRbtGlobal);
        const RigTForm invRefFrameY = inv(refFrameY);
        newRbt = targetRbt * RigTForm(Quat::makeXRotation(dy));
        newRbt = refFrameY * RigTForm(Quat::makeYRotation(-dx)) * invRefFrameY * newRbt;
        manipulationTarget->setRbt(newRbt);
      
      // Translation
      } else {
        newRbt = targetRbt * inv(m);
        manipulationTarget->setRbt(newRbt);
      }
    }
    
    glutPostRedisplay(); // we always redraw if we changed the scene
  }

  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}


static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system
  g_mouseDownX = x;
  g_mouseDownY = g_mouseClickY;

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;
  
  if (g_pickMode) {
    pick();
    g_pickMode = false;
    cout << "Picking mode is off" << endl;
  }
  
  glutPostRedisplay();
}

static void keyboardUp(const unsigned char key, const int x, const int y) {
  switch (key) {
  case ' ':
    g_spaceDown = false;
    break;
  }
  glutPostRedisplay();
}

#pragma mark Animation

static void sceneToCurrentKeyframe() {
  assert(!g_frames.empty());
  
  const int numNodes = static_cast<int>(g_frameNodes.size());
  const int numFrameNodes = static_cast<int>(g_currentFrameIt->size());
  assert(numFrameNodes == numNodes);
  
  // Apply each of the transforms in the current from
  for (int i = 0; i < numNodes; i++) {
    
    const shared_ptr<SgRbtNode> targetNode = g_frameNodes[i];
    (*g_currentFrameIt)[i] = targetNode->getRbt();
  }
}

static void currentKeyframeToScene() {
  
  assert(!g_frames.empty());
  
  const int numNodes = static_cast<int>(g_frameNodes.size());
  const Frame currentFrame = *g_currentFrameIt;
  const int numFrameNodes = static_cast<int>(currentFrame.size());
  assert(numFrameNodes == numNodes);
  
  // Apply each of the transforms in the current from
  for (int i = 0; i < numNodes; i++) {
    
    const shared_ptr<SgRbtNode> targetNode = g_frameNodes[i];
    const RigTForm frameRbt = (*g_currentFrameIt)[i];
    targetNode->setRbt(frameRbt);
  }
}

static void insertKeyframeAfter() {
  
  RigTForm ident;
  vector<RigTForm> newFrame(g_frameNodes.size(), ident);
  
  if (!g_frames.empty()) {
    auto next = g_currentFrameIt;
    next++;
    g_frames.insert(next, newFrame);
    g_currentFrameIt++;
  } else {
    g_frames.insert(g_frames.end(), newFrame);
    g_currentFrameIt = g_frames.begin();
  }
  sceneToCurrentKeyframe();
}

static void loadAnimFile(string fileName) {
  g_frames.clear();
  
  // First read num
  ifstream inFile(fileName.c_str());
  int tformsPerFrame;
  inFile >> tformsPerFrame;
  inFile.ignore();
  
  // Read frames
  string line;
  int readTForms = 0;
  vector<RigTForm> readingFrame;
  
  while(getline(inFile, line)) {
    array<double, 7> transQuat;
    std::stringstream ss(line);
    double val;
    int i = 0;
    
    while (ss >> val) {
      transQuat[i] = val;
      i++;
      
      if (ss.peek() == ',') {
        ss.ignore();
      }
      if (ss.peek() == '\n') {
        ss.ignore();
        assert(i == 7);
        break;
      }
    }
    
    RigTForm tform(transQuat);
    readingFrame.push_back(tform);
    readTForms++;
    
    if (readTForms == tformsPerFrame) {
      assert(readingFrame.size() == tformsPerFrame);
      g_frames.push_back(readingFrame);
      readTForms = 0;
      readingFrame.clear();
    }
  }
  
  g_currentFrameIt = g_frames.begin();
  currentKeyframeToScene();
}

static void writeAnimFile(string fileName) {
  ofstream outFile;
  outFile.open(fileName.c_str());
  
  // Write number of nodes in a frame
  const int frameNodes = static_cast<int>(g_frameNodes.size());
  outFile << frameNodes << "\n";
  
  // Write all frames
  for (const auto & frame : g_frames) {
    for (int i = 0; i < frameNodes; i++) {
      RigTForm tform = frame[i];
      array<double, 7> posQuat = tform.getPosQuat();
      ostream_iterator<double> outputIt(outFile, ",");
      copy(posQuat.begin(), posQuat.end(), outputIt);
      outFile << endl;
    }
  }
  
  outFile.close();

}

double spline(float unitT, double c0, double d0, double e0, double c1) {
  const float minusT = 1 - unitT;
  const float f = minusT * c0 + unitT * d0;
  const float g = minusT * d0 + unitT * e0;
  const float h = minusT * e0 + unitT * c1;
  const float m = minusT * f + unitT * g;
  const float n = minusT * g + unitT * h;
  return minusT * m + unitT * n;
}

Cvec3 spline(float unitT, Cvec3 c0, Cvec3 d0, Cvec3 e0, Cvec3 c1) {
  return Cvec3(spline(unitT, c0[0], d0[0], e0[0], c1[0]),
               spline(unitT, c0[1], d0[1], e0[1], c1[1]),
               spline(unitT, c0[2], d0[2], e0[2], c1[2]));
}

Quat spline(float unitT, const Quat& q0, const Quat& d0, const Quat& e0, const Quat& q1) {
  const Quat f = slerp(q0, d0, unitT);
  const Quat g = slerp(d0, e0, unitT);
  const Quat h = slerp(e0, q1, unitT);
  const Quat m = slerp(f, g, unitT);
  const Quat n = slerp(g, h, unitT);
  return slerp(m, n, unitT);
}

// Splines between c1 and c2 using C.R.
Cvec3 splineCR(float unitT, Cvec3 c0, Cvec3 c1, Cvec3 c2, Cvec3 c3) {
  const Cvec3 d1 = (c2 - c0) / 6.0 + c1;
  const Cvec3 e1 = (c3 - c1) / -6.0 + c2;
  return spline(unitT, c1, d1, e1, c2);
}

Quat splineCRQuat(float unitT, const Quat& c0, const Quat& c1, const Quat& c2, const Quat& c3) {
  
  const Quat d1 = pow(cn(c2 * inv(c0)), 1.0 / 6.0) * c1;
  const Quat e1 = pow(cn(c3 * inv(c1)), -1.0 / 6.0) * c2;
  return spline(unitT, c1, d1, e1, c2);
}


// Given t in the range [0, n], perform interpolation and draw the scene
// for the particular t. Returns true if we are at the end of the animation
// sequence, or false otherwise.
bool interpolateAndDisplay(float t) {
  
  const int framePrevI = static_cast<int>(floor(t)) + 1;
  const int frameNextI = static_cast<int>(floor(t)) + 2;
  
  if (frameNextI >= g_frames.size() - 1) {
    cout << "Stopped" << endl;
    return true;
  }
  
  double alpha = t - floor(t);
  
  list<Frame>::iterator itPrevPrev = g_frames.begin();
  list<Frame>::iterator itPrev = g_frames.begin();
  list<Frame>::iterator itNext = g_frames.begin();
  list<Frame>::iterator itNextNext = g_frames.begin();
  advance(itPrevPrev, framePrevI - 1);
  advance(itPrev, framePrevI);
  advance(itNext, frameNextI);
  advance(itNextNext, frameNextI + 1);
  
  const int numNodes = static_cast<int>(g_frameNodes.size());
  for (int i = 0; i < numNodes; i++) {
    
    const auto rbtPrevPrev = (*itPrevPrev)[i];
    const auto rbtPrev = (*itPrev)[i];
    const auto rbtNext = (*itNext)[i];
    const auto rbtNextNext = (*itNextNext)[i];
    
    const Cvec3 transPrevPrev = rbtPrevPrev.getTranslation();
    const Cvec3 transPrev = rbtPrev.getTranslation();
    const Cvec3 transNext = rbtNext.getTranslation();
    const Cvec3 transNextNext = rbtNextNext.getTranslation();
    const Quat rotPrevPrev = rbtPrevPrev.getRotation();
    const Quat rotPrev = rbtPrev.getRotation();
    const Quat rotNext = rbtNext.getRotation();
    const Quat rotNextNext = rbtNextNext.getRotation();
    
    
    // LERP:
    // const Cvec3 newTrans = transPrev * (1.0 - alpha) + transNext * alpha;
    // const Quat newRot = slerp(rotPrev, rotNext, alpha);
    
    // SPLINE:
    const Cvec3 newTrans = splineCR(alpha, transPrevPrev, transPrev, transNext, transNextNext);
    const Quat newRot = splineCRQuat(alpha, rotPrevPrev, rotPrev, rotNext, rotNextNext);
    //const Quat newRot = spline(alpha, rotPrev, rotPrev, rotNext, rotNext);
    
    const RigTForm newRbt(newTrans, newRot);
    
    const shared_ptr<SgRbtNode> targetNode = g_frameNodes[i];
    targetNode->setRbt(newRbt);
  }
  
  return false;
}
// Interpret "ms" as milliseconds into the animation
static void animateTimerCallback(int ms) {
  float t = (float)ms/(float)g_msBetweenKeyFrames;
  bool endReached = false; //interpolateAndDisplay(t) || !g_isPlaying;
  if (!endReached) {
    g_blobY = 12.0 + sin(t * 8.0) * 3.0;
    updateBlobGeometry();
    glutTimerFunc(1000/g_animateFramesPerSecond,
                  animateTimerCallback,
                  ms + 1000/g_animateFramesPerSecond);
  } else {
    //g_isPlaying = false;
    //auto endIt = g_frames.end();
    //endIt--;
    //endIt--;
    //g_currentFrameIt = endIt;
    //currentKeyframeToScene();
  }
  glutPostRedisplay();
}

#pragma mark Keyboard
static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC
  case 'h':
    cout << " ============== H E L P ==============\n"
    << "h\t\thelp menu\n"
    << "p\t\tPick mode (slightly modified)\n"
    << "left/right\tradius\n"
    << "up/down\tnumber solver iterations\n"
    << "w\t\ttoggle wireframe\n"
    << "s\t\ttoggle hard edge detection\n"
    << "1/2/3/4\t\tswitch shapes\n"
    << "-/=\t\tchange the thickness of the torus\n"
      << "5/6/7\t\t(surface/hermite/normals) display\n\n" << endl;
    break;
  case '1':
      g_model = 0;
      updateBlobGeometry();
      break;
  case '2':
    g_model = 1;
    updateBlobGeometry();
    break;
  case '3':
    g_model = 2;
    g_sharp = true;
    updateBlobGeometry();
    break;
  case '4':
    g_model = 3;
    g_sharp = true;
    updateBlobGeometry();
    break;
  case '5':
      g_debug = 0;
    updateBlobGeometry();
    break;
  case '6':
    g_debug = 1;
    updateBlobGeometry();
    break;
      
  case '7':
    g_debug = 2;
    updateBlobGeometry();
    break;
      
  case '=':
      g_thickness *= 1.05;
      cout << "Thickness (torus only) = " << g_thickness << std::endl;
      updateBlobGeometry();
      break;
    case '-':
      g_thickness /= 1.05;
      cout << "Thickness (torus only) = " << g_thickness << std::endl;
      updateBlobGeometry();
      break;
      
    
  case 'w':
      g_wireframe = !g_wireframe;
      if (g_wireframe) {
        g_redDiffuseMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_LINE);
      } else {
        g_redDiffuseMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_FILL);
      }
      break;
  case 's':
      g_sharp = !g_sharp;
      updateBlobGeometry();
    break;
  case 'v':
    g_eyeMode = (g_eyeMode + 1) % 3;
    switch (g_eyeMode) {
      case 0:
        g_eyeNode = g_skyNode;
        cout << "Active eye is Sky" << endl;
        break;
      case 1:
        g_eyeNode = g_robot1Node;
        cout << "Active eye is Robot 1" << endl;
        break;
      case 2:
        g_eyeNode = g_robot2Node;
        cout << "Active eye is Robot 2" << endl;
        break;
    }
    break;
    
  case 'p':
    g_pickMode = true;
    cout << "Picking mode is on" << endl;
    break;

  case 'm':
    g_egoMode = !g_egoMode;
    if (g_egoMode) {
      cout << "Editing sky eye w.r.t. sky-sky frame" << endl;
    } else {
      cout << "Editing sky eye w.r.t. world-sky frame" << endl;
    }
    break;
    
  // Copy current key frame to scene graph
  case 'c':
    if (!g_frames.empty()) {
      currentKeyframeToScene();
      cout << "Cleared changes to keyframe" << endl;
    }
      
    break;
      
  // Copy scene graph to current frame
  case 'u':
    
    if (g_frames.empty()) {
      insertKeyframeAfter();
      cout << "Created keyframe" << endl;
    } else {
      sceneToCurrentKeyframe();
      cout << "Updated current keyframe" << endl;
    }
    
    break;
  
  case '.': {
    
    list<Frame>::iterator lastFrame = g_frames.end();
    lastFrame--;
    if (!g_frames.empty() && g_currentFrameIt != lastFrame) {
      g_currentFrameIt++;
      currentKeyframeToScene();
    }
    break;
  }
  case ',':
      
    if (!g_frames.empty() && g_currentFrameIt != g_frames.begin()) {
      g_currentFrameIt--;
      currentKeyframeToScene();
    }
    break;
  
  case 'd':
    if (!g_frames.empty()) {
      cout << "Deleted keyframe" << endl;
      auto prev = g_currentFrameIt;
      auto next = g_currentFrameIt;
      prev--;
      next++;
      bool isFirst = g_currentFrameIt == g_frames.begin();
      g_frames.erase(g_currentFrameIt);
      if (g_frames.empty()) {
        g_currentFrameIt = g_frames.end();
      } else {
        if (!isFirst) {
          g_currentFrameIt = prev;
        } else {
          g_currentFrameIt = next;
        }
        currentKeyframeToScene();
      }
    }
    break;
  
  case 'n': {
    insertKeyframeAfter();
    cout << "Inserted keyframe" << endl;
    break;
  }
      
  case 'i':{
    loadAnimFile("gabe.anim");
    cout << "Loaded from file" << endl;
    break;
  }
      
  case 'z':{
    writeAnimFile("gabe.anim");
    cout << "Saved to file" << endl;
    break;
  }
      
  case '+': {
    if (!g_isPlaying) {
      g_msBetweenKeyFrames += 250;
    }
    break;
  }
      
  case '_': {
    if (!g_isPlaying && g_msBetweenKeyFrames >= 500) {
      g_msBetweenKeyFrames -= 250;
    }
    break;
  }
      
  case 'y': {
    if ((g_frames.size() >= 4 && !g_isPlaying) || true) {
      cout << "Playing.." << endl;
      g_isPlaying = true;
      animateTimerCallback(0);
    } else if (g_isPlaying) {
      g_isPlaying = false;
    }
    break;
  }
    
  case ' ':
    g_spaceDown = true;
    break;
  }
  glutPostRedisplay();
}

// new  special keyboard callback, for arrow keys
static void specialKeyboard(const int key, const int x, const int y) {
  switch (key) {
    case GLUT_KEY_RIGHT:
      g_blobX += 0.5;
      break;
    case GLUT_KEY_LEFT:
      g_blobX -= 0.5;
      break;
    case GLUT_KEY_UP:
      //g_radius *= 1.05;
      g_tolerance *= 1.05;
      cout << "g_tolerance = " << g_tolerance << endl;
      break;
    case GLUT_KEY_DOWN:
      //g_radius /= 1.05;
      g_tolerance /= 1.05;
      cout << "g_tolerance = " << g_tolerance << endl;
      break;
  }
  updateBlobGeometry();
  glutPostRedisplay();
}

static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
#ifdef __MAC__
  glutInitDisplayMode(GLUT_3_2_CORE_PROFILE|GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH); // core profile flag is required for GL 3.2 on Mac
#else
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
#endif
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("Assignment 8");                       // title the window

  glutIgnoreKeyRepeat(true);                              // avoids repeated keyboard calls when holding space to emulate middle mouse

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
  glutKeyboardUpFunc(keyboardUp);
  glutSpecialFunc(specialKeyboard);
}

static void initGLState() {
  glClearColor(0, 0, 0, 0);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
  if (!g_Gl2Compatible)
    glEnable(GL_FRAMEBUFFER_SRGB);
}

static void initMaterials() {
  // Create some prototype materials
  Material diffuse("./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader");
  Material specular("./shaders/basic-gl3.vshader", "./shaders/specular-gl3.fshader");
  Material solid("./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader");
  
  // copy diffuse prototype and set red color
  g_redDiffuseMat.reset(new Material(specular));
  g_redDiffuseMat->getUniforms().put("uColor", Cvec3f(1, 0.65, 0));
  g_redDiffuseMat->getRenderStates().disable(GL_CULL_FACE);
  
  // copy diffuse prototype and set blue color
  g_blueDiffuseMat.reset(new Material(diffuse));
  g_blueDiffuseMat->getUniforms().put("uColor", Cvec3f(0, 0, 1));
  
  // normal mapping material
  g_bumpFloorMat.reset(new Material("./shaders/normal-gl3.vshader", "./shaders/normal-gl3.fshader"));
  g_bumpFloorMat->getUniforms().put("uTexColor", shared_ptr<ImageTexture>(new ImageTexture("Fieldstone.ppm", true)));
  g_bumpFloorMat->getUniforms().put("uTexNormal", shared_ptr<ImageTexture>(new ImageTexture("FieldstoneNormal.ppm", false)));
  
  // copy solid prototype, and set to wireframed rendering
  g_arcballMat.reset(new Material(solid));
  g_arcballMat->getUniforms().put("uColor", Cvec3f(0.27f, 0.82f, 0.35f));
  g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_LINE);
  
  // copy solid prototype, and set to color white
  g_lightMat.reset(new Material(solid));
  g_lightMat->getUniforms().put("uColor", Cvec3f(1, 1, 1));
  
  // pick shader
  g_pickingMat.reset(new Material("./shaders/basic-gl3.vshader", "./shaders/pick-gl3.fshader"));
};

static void initGeometry() {
  initGround();
  initCubes();
  initSphere();
  initBlobbyMeshes();
}

static void initScene() {
  g_world.reset(new SgRootNode());
  
  g_skyNode.reset(new SgRbtNode(RigTForm(Cvec3(0.0, 0.25, 4.0))));
  
  g_eyeNode = g_skyNode;
  
  g_groundNode.reset(new SgRbtNode(RigTForm(Cvec3(0, -3, 0))));
  //g_groundNode->addChild(shared_ptr<MyShapeNode>(
  //                                               new MyShapeNode(g_ground, g_bumpFloorMat, Cvec3(0, g_groundY, 0))));
  
  g_light1Node.reset(new SgRbtNode(RigTForm(Cvec3(-2, 1.5, 2))));
  g_light1Node->addChild(shared_ptr<MyShapeNode>(
                                                 new MyShapeNode(g_sphere, g_lightMat,
                                                   Cvec3(0, 0, 0),
                                                   Cvec3(0, 0, 0),
                                                   Cvec3(0.25, 0.25, 0.25))));

  g_light2Node.reset(new SgRbtNode(RigTForm(Cvec3(3, 2, 1))));
  g_light2Node->addChild(shared_ptr<MyShapeNode>(
                                                 new MyShapeNode(g_sphere, g_lightMat,
                                                   Cvec3(0, 0, 0),
                                                   Cvec3(0, 0, 0),
                                                   Cvec3(0.25, 0.25, 0.25))));
  
  g_robot1Node.reset(new SgRbtNode(RigTForm(Cvec3(-2, 1, -8))));
  g_robot2Node.reset(new SgRbtNode(RigTForm(Cvec3(2, 1, -8))));

  // create a single transform node for both the blob and the blob shells
  g_blobNode.reset(new SgRbtNode(RigTForm(Cvec3(0, 0, 0))));
  
  // add each shell as shape node
  g_blobNode->addChild(shared_ptr<MyShapeNode>(new MyShapeNode(g_blobGeometry, g_redDiffuseMat)));
  updateBlobGeometry();
  
  g_world->addChild(g_skyNode);
  g_world->addChild(g_groundNode);
  g_world->addChild(g_robot1Node);
  g_world->addChild(g_robot2Node);
  g_world->addChild(g_light1Node);
  g_world->addChild(g_light2Node);
  g_world->addChild(g_blobNode);
}

static void initAnimation() {
  dumpSgRbtNodes(g_world, g_frameNodes);
  g_currentFrameIt = g_frames.end();
}

int main(int argc, char * argv[]) {
  try {
    initGlutState(argc,argv);

    // on Mac, we shouldn't use GLEW.

#ifndef __MAC__
    glewInit(); // load the OpenGL extensions
#endif

    cout << (g_Gl2Compatible ? "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.5") << endl;

#ifndef __MAC__
    if ((!g_Gl2Compatible) && !GLEW_VERSION_3_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
    else if (g_Gl2Compatible && !GLEW_VERSION_2_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.0");
#endif
    
    initGLState();
    initMaterials();
    initGeometry();
    initScene();
    initAnimation();

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
