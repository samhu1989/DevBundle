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


#ifndef MESHVIEWERWIDGETT_H
#define MESHVIEWERWIDGETT_H


//== INCLUDES =================================================================

#include <string>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Utils/GenProg.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <OpenMesh/Tools/Utils/StripifierT.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include "QGLViewerWidget.h"
#include "MeshColor.h"
#include <memory>


//== FORWARDS =================================================================

class QImage;


//== CLASS DEFINITION =========================================================

	      
template <typename M>
class MeshViewerWidgetT : public QGLViewerWidget
{
public:

  typedef M                             Mesh;
  typedef OpenMesh::StripifierT<Mesh>   MyStripifier;
public:

  /// default constructor
  MeshViewerWidgetT(QWidget* _parent=0)
    : QGLViewerWidget(_parent),
      f_strips_(false), 
      tex_id_(0),
      tex_mode_(GL_MODULATE),
      mesh_ptr_(new Mesh),
      strips_(new MyStripifier(*mesh_ptr_)),
      custom_color_(new MeshColor<Mesh>(*mesh_ptr_)),
      use_color_(true),
      show_vnormals_(false),
      show_fnormals_(false)
  {
    QAction* a = add_draw_mode("Points");
    slotDrawMode(a);
    add_draw_mode("Hidden-Line");
#if defined(OM_USE_OSG) && OM_USE_OSG
    add_draw_mode("OpenSG Indices");
#endif
  }
  
  /// destructor
  ~MeshViewerWidgetT() {}
  
public:

  /// open mesh
  virtual bool open_mesh(const char* _filename, OpenMesh::IO::Options _opt);
  virtual bool save_mesh(const char* _filename, OpenMesh::IO::Options _opt);
  
  /// load texture
  virtual bool open_texture( const char *_filename );
  bool set_texture( QImage& _texsrc );
 
  void enable_strips();
  void disable_strips();

  virtual void set_center_at_mesh(const Mesh&);
  
  Mesh& mesh() { return *mesh_ptr_; }
  const Mesh& mesh() const { return *mesh_ptr_; }

  std::shared_ptr<Mesh>& mesh_ptr(){ return mesh_ptr_;}
  void reset_ptr(std::shared_ptr<Mesh>& ptr)
  {
      mesh_ptr_ = ptr;
      std::cerr<<"setting up centers:"<<std::endl;
      if(mesh_ptr_&&mesh_ptr_->n_vertices()!=0)set_center_at_mesh(*mesh_ptr_);
      strips_.reset(new MyStripifier(*mesh_ptr_));
      custom_color_.reset(new MeshColor<Mesh>(*mesh_ptr_));
      updateGL();
  }
  
protected:
  
  /// inherited drawing method
  virtual void draw_scene(const std::string& _draw_mode);
  
protected:
  
  /// draw the mesh
  virtual void draw_openmesh(const std::string& _drawmode);


  void glVertex( const typename Mesh::VertexHandle _vh )
  { glVertex3fv( &mesh_ptr_->point( _vh )[0] ); }

  void glVertex( const typename Mesh::Point& _p )
  { glVertex3fv( &_p[0] ); }
  
  void glNormal( const typename Mesh::VertexHandle _vh )
  { glNormal3fv( &mesh_ptr_->normal( _vh )[0] ); }

  void glTexCoord( const typename Mesh::VertexHandle _vh )
  { glTexCoord2fv( &mesh_ptr_->texcoord(_vh)[0] ); }
  
  void glColor( const typename Mesh::VertexHandle _vh )
  { glColor3ubv( &mesh_ptr_->color(_vh)[0] ); }
  
  // face properties

  void glNormal( const typename Mesh::FaceHandle _fh )
  { glNormal3fv( &mesh_ptr_->normal( _fh )[0] ); }

  void glColor( const typename Mesh::FaceHandle _fh )
  { glColor3ubv( &mesh_ptr_->color(_fh)[0] ); }

  void glMaterial( const typename Mesh::FaceHandle _fh, 
		   int _f=GL_FRONT_AND_BACK, int _m=GL_DIFFUSE )
  { 
    OpenMesh::Vec3f c=OpenMesh::color_cast<OpenMesh::Vec3f>(mesh_ptr_->color(_fh));
    OpenMesh::Vec4f m( c[0], c[1], c[2], 1.0f );

    glMaterialfv(_f, _m, &m[0]); 
  }


protected: // Strip support
  
  void compute_strips(void)
  {
    if (f_strips_)
    {
      strips_->clear();
      strips_->stripify();
    }
  }    

protected: // inherited
   
  virtual void keyPressEvent( QKeyEvent* _event);

protected:
   
  bool                   f_strips_; // enable/disable strip usage
  GLuint                 tex_id_;
  GLint                  tex_mode_;
  OpenMesh::IO::Options  opt_; // mesh file contained texcoords?
  
  std::shared_ptr<Mesh>  mesh_ptr_;
  std::shared_ptr<MeshColor<Mesh>> custom_color_;
  std::shared_ptr<MyStripifier>    strips_;
  bool                   use_color_;
  bool                   show_vnormals_;
  bool                   show_fnormals_;
  float                  normal_scale_;
  OpenMesh::FPropHandleT< typename Mesh::Point > fp_normal_base_;
};


//=============================================================================
#include "MeshViewerWidgetT.hpp"
//=============================================================================
#endif // MESHVIEWERWIDGETT_H defined
//=============================================================================

