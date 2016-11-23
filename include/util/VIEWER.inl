/*
This file is part of ESFC (Eulerian Solid-FLuid Coupling).

ESFC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESFC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ESFC.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <util/SIMPLE_PARSER.h>
#include <libpng16/png.h>
template<class T>
void VIEWER<T>::init()
{
	Real moveSpeed = 1.0;
  Real eyeDistanceScale = 1.0;

  glvu.SetMoveSpeed( moveSpeed * glvu.GetMoveSpeed() );

  int argc = 0;
  glutInit(&argc, NULL);

  glvuVec3f boxCenter(0, 0, 0);
  glvuWindow( boxCenter, eyeDistanceScale );
}
template<class T>
void VIEWER<T>::displayFunc()
{
	static GLfloat mat_shininess[] = { 120.0 };

	static GLfloat mat_diffuse[] = { 0.3, 0.3, 0.3, 1.0 };
	static GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };

  Camera* camera = glvu.GetCurrentCam();
  glvuVec3f Eye, Lookat, Up;
  camera->GetLookAtParams(&Eye, &Lookat, &Up);

  // repack camera settings into arrays
  float eye[] = {Eye.x, Eye.y, Eye.z};
  float look[3];
  look[0] = Lookat.x - eye[0];
  look[1] = Lookat.y - eye[1];
  look[2] = Lookat.z - eye[2];
  float magnitude = 1.0f / sqrt(look[0] * look[0] + look[1] * look[1] + look[2] * look[2]);
  look[0] *= magnitude;
  look[1] *= magnitude;
  look[2] *= magnitude;

  glvu.BeginFrame();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
    
  if(showGrid)
    drawGrid();

  simulator->display();

  glvu.EndFrame();
}

template<class T>
void VIEWER<T>::drawGrid()
{
  int xRes = simulator->xRes();
  int yRes = simulator->yRes();
  int zRes = simulator->zRes();
  Real dh = simulator->dh();
  const VEC3F& origin = simulator->origin();
  Real xlength = xRes * dh;
  Real ylength = yRes * dh;
  Real zLength = zRes * dh;

  glColor4f(1, 1, 1, 1);
  glPushMatrix();
  glTranslatef(origin[0] + xlength / 2, origin[1] + ylength / 2, origin[2] + zLength / 2);
  glScalef(xlength, ylength, zLength);
  glutWireCube(1.0);
  glPopMatrix();
}

template<class T>
void VIEWER<T>::keyboardFunc(unsigned char Key, int x, int y)
{
  simulator->keyboardFunc(Key);
  switch(Key)
  {
  	case 'a':{
  		animate = !animate;
  		break;
  	}
    case 's':{
      step = !step;
      break;
    }
    case 'g':{
      showGrid = !showGrid;
      break;
    }
    case 'v':
      {
        Camera* camera = glvu.GetCurrentCam();
        glvuVec3f eye;
        glvuVec3f lookat;
        glvuVec3f up;
        camera->GetLookAtParams(&eye, &lookat, &up);
        cout << " Eye(" << eye[0] << ", " << eye[1] << ", " << eye[2] << "), " ;
        cout << " LookAtCntr(" << lookat[0] << ", " << lookat[1] << ", " << lookat[2] << "), " ;
        cout << " Up(" << up[0] << ", " << up[1] << ", " << up[2] << ");\n";

        cout << "eye = " << eye[0] << "," << eye[1] << "," << eye[2]
        		 << "\nlookatcntr = " << lookat[0] << "," << lookat[1] << "," << lookat[2]
        		 << "\nup = " << up[0] << "," << up[1] << "," << up[2] << "\n";

      }
      break;
    case 'q':
    case 'Q':
      TIMING_BREAKDOWN::printTimingBreakdown();
      exit(0);
      break;
  };

  glutPostRedisplay();
  if (Key != '=')
    glvu.Keyboard(Key,x,y);
}

template<class T>
void VIEWER<T>::idleFunc()
{
  static bool getScreenshot = SIMPLE_PARSER::getBool("screen shot", false);
  static int saveEvery = SIMPLE_PARSER::getInt("save every", 1);
	if(animate){
    if(getScreenshot){
      if(simulator->currentFrame % saveEvery == 0)
        screenshot(simulator->renderPath, simulator->currentFrame / saveEvery);
    }

		simulator->step();
    if(step)
      animate = false;
    glutPostRedisplay(); 
  }
}

template<class T>
void VIEWER<T>::mouseFunc(int button, int state, int x, int y)
{
  static float clickZ;
  int Modifiers = glutGetModifiers();
  if (button == GLUT_LEFT_BUTTON && 
      state == GLUT_DOWN &&
      Modifiers & GLUT_ACTIVE_SHIFT)
  {
    // retrieve and store the depth of this click
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glReadPixels(x, viewport[3] - y, 
                 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &clickZ);

    // get the world space coordinate
    VEC3F point = unproject(x,y,clickZ);

    // hand the coordinate to the tet mesh
    // simulator->click(point);
    // mouseClicked = true;
    return;
  }
  if (button == GLUT_LEFT_BUTTON && 
      state == GLUT_UP)
  {
    // mouseClicked = false;
    // simulator->unclick();
    // integrator->unclick();
    return;
  } 
	// pass through to default handler
  glvu.Mouse(button,state,x,y);
}

template<class T>
void VIEWER<T>::motionFunc(int x, int y)
{
	glvu.Motion(x,y);
}

template<class T>
int VIEWER<T>::glvuWindow(glvuVec3f bboxCenter, Real eyeDistanceScale)
{
	glvu.Init((char *)"GLVU Window",
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            windowStartX, windowStartY, windowWidth, windowHeight);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = {10.0, 10.0, 10.0, 1.0};
  GLfloat lightZeroColor[] = {0.8, 0.8, 0.8, 1.0};
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  //glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glShadeModel(GL_SMOOTH);
  //glShadeModel(GL_FLAT);
  glClearColor(0, 0, 0, 0);

  glutDisplayFunc(VIEWER<T>::displayFunc);
  glutMouseFunc(VIEWER<T>::mouseFunc);
  glutMotionFunc(VIEWER<T>::motionFunc);
  glutKeyboardFunc(VIEWER<T>::keyboardFunc);
  glutIdleFunc(VIEWER<T>::idleFunc);

  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10);

  VEC3F myEye = SIMPLE_PARSER::getVEC3F("eye", VEC3F(0, 0.3, 3));
  glvuVec3f Eye(myEye[0], myEye[1], myEye[2]);

  VEC3F myLookAtCntr = SIMPLE_PARSER::getVEC3F("lookatcntr", VEC3F(0, 0.2, 2));
  glvuVec3f LookAtCntr(myLookAtCntr[0], myLookAtCntr[1], myLookAtCntr[2]);

  VEC3F myUp = SIMPLE_PARSER::getVEC3F("up", VEC3F(0, 0.9, -0.2));
  glvuVec3f Up(myUp[0], myUp[1], myUp[2]);
  
  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov,Aspect, Near, Far);

  // reset center
  glvuVec3f center(0.0f, 0.0f, 0.0f);
  glvu.SetWorldCenter(center);
  
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

bool save_png_libpng(const char *filename, uint8_t *pixels, int w, int h)
{
  png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
  if (!png)
    return false;

  png_infop info = png_create_info_struct(png);
  if (!info) {
    png_destroy_write_struct(&png, &info);
    return false;
  }

  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    png_destroy_write_struct(&png, &info);
    return false;
  }

  png_init_io(png, fp);
  png_set_IHDR(png, info, w, h, 8 /* depth */, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
    PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  png_colorp palette = (png_colorp)png_malloc(png, PNG_MAX_PALETTE_LENGTH * sizeof(png_color));
  if (!palette) {
    fclose(fp);
    png_destroy_write_struct(&png, &info);
    return false;
  }
  png_set_PLTE(png, info, palette, PNG_MAX_PALETTE_LENGTH);
  png_write_info(png, info);
  png_set_packing(png);

  png_bytepp rows = (png_bytepp)png_malloc(png, h * sizeof(png_bytep));
  for (int i = 0; i < h; ++i)
    rows[i] = (png_bytep)(pixels + (h - i - 1) * w * 3);

  png_write_image(png, rows);
  png_write_end(png, info);
  png_free(png, palette);
  png_destroy_write_struct(&png, &info);

  fclose(fp);
  delete[] rows;
  return true;
}
//////////////////////////////////////////////////////////////////////////////
// Dump a frame with this number to this directory
//////////////////////////////////////////////////////////////////////////////
template<class T>
void VIEWER<T>::screenshot(string renderPath, int frame)
{
  // FILE *fp;

  string number = IO::itoPaddedString(frame);
  char FileName[256];
  sprintf(FileName,"%srenderGL.%s.png", renderPath.c_str(), number.c_str());

  GLint OldReadBuffer;
  glGetIntegerv(GL_READ_BUFFER,&OldReadBuffer);
  glReadBuffer(GL_FRONT);

  GLint OldPackAlignment;
  glGetIntegerv(GL_PACK_ALIGNMENT,&OldPackAlignment); 
  glPixelStorei(GL_PACK_ALIGNMENT,1);

  int WW = glutGet((GLenum)GLUT_WINDOW_WIDTH);
  int WH = glutGet((GLenum)GLUT_WINDOW_HEIGHT);
  int NumPixels = WW*WH;
  GLubyte* Pixels = new GLubyte[NumPixels*3];
  if (Pixels==NULL) { printf("UNABLE TO ALLOC PIXEL READ ARRAY!\n"); return; }
  glReadPixels(0,0,WW,WH,GL_RGB,GL_UNSIGNED_BYTE,Pixels);

  save_png_libpng(FileName, Pixels, WW, WH);

  delete[] Pixels;

  glPixelStorei(GL_PACK_ALIGNMENT,OldPackAlignment);
  glReadBuffer((GLenum)OldReadBuffer);
}

//////////////////////////////////////////////////////////////////////////////
// Translate screen space to world space
//
// Adapted from the NeHe page:
// http://nehe.gamedev.net/data/articles/article.asp?article=13
//////////////////////////////////////////////////////////////////////////////
template<class T>
VEC3F VIEWER<T>::unproject(float x, float y, float z)
{
  GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];

	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

  double worldX, worldY, worldZ;
	gluUnProject(x, viewport[3] - y, z,
               modelview, projection, viewport, 
               &worldX, &worldY, &worldZ);

  return VEC3F(worldX, worldY, worldZ);
}


