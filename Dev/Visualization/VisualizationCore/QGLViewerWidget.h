/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */

/*===========================================================================*\
 *                                                                           *             
 *   $Revision: 1258 $                                                         *
 *   $Date: 2015-04-28 15:07:46 +0200 (Di, 28 Apr 2015) $                   *
 *                                                                           *
\*===========================================================================*/


#ifndef QGLVIEWERWIDGET_H
#define QGLVIEWERWIDGET_H


//== INCLUDES =================================================================

#include "visualizationcore_global.h"
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <QGLWidget>
#include <string>
#include <vector>
#include <map>
#include <armadillo>
#include "qglpointselection.h"


//== FORWARD DECLARATIONS =====================================================

class QMenu;
class QActionGroup;
class QAction;

//== CLASS DEFINITION =========================================================

  
class VISUALIZATIONCORESHARED_EXPORT QGLViewerWidget : public QGLWidget
{

  Q_OBJECT
  
public:
  typedef QGLWidget Super;
   
  // Default constructor.
  QGLViewerWidget( QWidget* _parent=0 );

  // 
  QGLViewerWidget( QGLFormat& _fmt, QWidget* _parent=0 );

  // Destructor.
  virtual ~QGLViewerWidget();
public slots:
  void setBackgroundColor(QColor);
  void rect_select(bool select){
      rect_select_=select;
      rect_selecting_ = false;
  }
signals:
  void have_been_transfomed(arma::fmat R,arma::fvec t,int index);

public:
 void notify_been_transfomed(arma::fmat R,arma::fvec t,int index)
 {
     emit have_been_transfomed(R,t,index);
 }


private:

  virtual void init(void);

public:
  static bool draw_coordinate_;

  /* Sets the center and size of the whole scene. 
     The _center is used as fixpoint for rotations and for adjusting
     the camera/viewer (see view_all()). */
  void set_scene_pos( const OpenMesh::Vec3f& _center, float _radius );  

  /* view the whole scene: the eye point is moved far enough from the
     center so that the whole scene is visible. */
  void view_all();

  /// add draw mode to popup menu, and return the QAction created
  QAction *add_draw_mode(const std::string& _s);

  /// delete draw mode from popup menu
  void del_draw_mode(const std::string& _s);

  // set draw mode to
  void set_draw_mode(const std::string& _s);

  const std::string& current_draw_mode() const 
  { return draw_mode_ ? draw_mode_names_[draw_mode_-1] : nomode_; }

  float radius() const { return radius_; }

  const OpenMesh::Vec3f& center() const { return center_; }

  const GLdouble* modelview_matrix() const  { return modelview_matrix_;  }
  const GLdouble* projection_matrix() const { return projection_matrix_; }

  float fovy() const { return 45.0f; }

  QAction* findAction(const char *name);
  void addAction(QAction* action, const char* name);
  void removeAction(const char* name);
  void removeAction(QAction* action);

protected:

  // draw the scene: will be called by the paintGL() method.
  virtual void draw_scene(const std::string& _draw_mode);

  double performance(void);
  
  void setDefaultMaterial(void);
  void setDefaultLight(void);

public slots:

  // popup menu clicked
  void slotDrawMode(QAction *_mode);  
  void slotSnapshot( void );

  
private: // inherited

  // initialize OpenGL states (triggered by Qt)
  void initializeGL();

  // draw the scene (triggered by Qt)
  void paintGL();

  // handle resize events (triggered by Qt)
  void resizeGL( int w, int h );

protected:
   
  // Qt mouse events
  virtual void mousePressEvent( QMouseEvent* );
  virtual void mouseReleaseEvent( QMouseEvent* );
  virtual void mouseMoveEvent( QMouseEvent* );
  virtual void wheelEvent( QWheelEvent* );
  virtual void keyPressEvent( QKeyEvent* );

  virtual void add_box(void){}
  virtual void add_sub_box(void){}
  virtual void del_box(void){}
  virtual bool mod_box(void){}
  virtual void next_box(void){}

protected:
  void draw_coordinate();
  virtual void processSelections();
  PointSelections selections_;
  float            radius_;
  float            point_size_;
private:
   
  // updates projection matrix
  void update_projection_matrix();

  // translate the scene and update modelview matrix
  void translate(const OpenMesh::Vec3f& _trans);

  // rotate the scene (around its center) and update modelview matrix
  void rotate(const OpenMesh::Vec3f& _axis, float _angle);

  OpenMesh::Vec3f  center_;
	      
  GLdouble    projection_matrix_[16],
              modelview_matrix_[16];


  // popup menu for draw mode selection
  QMenu*               popup_menu_;
  QActionGroup*        draw_modes_group_;
  typedef std::map<QString,QAction*> ActionMap;
  ActionMap            names_to_actions;
  unsigned int              draw_mode_;
  unsigned int              n_draw_modes_;
  std::vector<std::string>  draw_mode_names_;
  static std::string        nomode_;



  // virtual trackball: map 2D screen point to unit sphere
  bool map_to_sphere(const QPoint& _point, OpenMesh::Vec3f& _result);
  bool map_to_far(const QPoint& _point, OpenMesh::Vec3f& _result);
  void map_to_world(const arma::fvec&,arma::fvec&);

  QPoint           last_point_2D_;
  OpenMesh::Vec3f  last_point_3D_;
  bool             last_point_ok_;
  bool             show_coordinates_;

  bool             rect_select_;
  bool             rect_selecting_;
  QPoint           rect_point_2D_;

  float coord_[18];
  GLushort coord_index_[6];
  uint8_t coord_color_[18];
  float blight;
};


//=============================================================================
#endif // QGLVIEWERWIDGET_H
//=============================================================================

