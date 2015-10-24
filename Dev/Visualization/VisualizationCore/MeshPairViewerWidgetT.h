#ifndef MESHPAIRVIEWERWIDGET_H
#define MESHPAIRVIEWERWIDGET_H
#include <memory>
#include <string>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/Utils/GenProg.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <OpenMesh/Tools/Utils/StripifierT.hh>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include "QGLViewerWidget.h"
template <typename M>
class MeshPairViewerWidgetT : public QGLViewerWidget
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
      strips_(mesh_),
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
  virtual bool open_mesh(const char*,Mesh&,OpenMesh::IO::Options);
  virtual bool save_mesh(const char*,const Mesh&,OpenMesh::IO::Options);

  /// load texture
  virtual bool open_texture( const char *_filename );
  bool set_texture( QImage& _texsrc );

  void enable_strips();
  void disable_strips();

  Mesh& first() { return *mesh_first; }
  const Mesh& first() const { return *mesh_first; }

  Mesh& second() { return *mesh_second; }
  const Mesh& second() const { return *mesh_second; }

  std::shared_ptr<Mesh> first_ptr() { return mesh_first; }
  std::shared_ptr<Mesh> _ptr() { return mesh_second; }

protected:

  /// inherited drawing method
  virtual void draw_scene(const std::string& _draw_mode);

protected:

  /// draw the mesh
  virtual void draw_openmesh(const Mesh& mesh_,const std::string& _drawmode);


  void glVertex( const Mesh& m,  const typename Mesh::VertexHandle _vh)
  { glVertex3fv( &m.point( _vh )[0] ); }

  void glVertex( const Mesh &m, const typename Mesh::Point& _p)
  { glVertex3fv( &m._p[0] ); }

  void glNormal( const Mesh &m,  const typename Mesh::VertexHandle _vh )
  { glNormal3fv( &m.normal( _vh )[0] ); }

  void glTexCoord( const Mesh &m, const typename Mesh::VertexHandle _vh )
  { glTexCoord2fv( &m.texcoord(_vh)[0] ); }

  void glColor( const Mesh &m, const typename Mesh::VertexHandle _vh )
  { glColor3ubv( &m.color(_vh)[0] ); }

  // face properties

  void glNormal( const Mesh &m, const typename Mesh::FaceHandle _fh )
  { glNormal3fv( &m.normal( _fh )[0] ); }

  void glColor( const Mesh &m, const typename Mesh::FaceHandle _fh )
  { glColor3ubv( &m.color(_fh)[0] ); }

  void glMaterial( const Mesh &mesh, const typename Mesh::FaceHandle _fh,
           int _f=GL_FRONT_AND_BACK, int _m=GL_DIFFUSE )
  {
    OpenMesh::Vec3f c=OpenMesh::color_cast<OpenMesh::Vec3f>(mesh.color(_fh));
    OpenMesh::Vec4f m( c[0], c[1], c[2], 1.0f );

    glMaterialfv(_f, _m, &m[0]);
  }


protected: // Strip support

  void compute_strips(void)
  {
    if (f_strips_)
    {
      strips_first.clear();
      strips_first.stripify();
      strips_second.clear();
      strips_second.stripify();
    }
  }

protected: // inherited

  virtual void keyPressEvent( QKeyEvent* _event);

protected:

  bool                   f_strips_; // enable/disable strip usage
  GLuint                 tex_id_;
  GLint                  tex_mode_;
  OpenMesh::IO::Options  opt_; // mesh file contained texcoords?

  std::shared_ptr<Mesh>  mesh_first;
  std::shared_ptr<Mesh>  mesh_second;
  MyStripifier           strips_first;
  MyStripifier           strips_second;
  bool                   use_color_;
  bool                   show_vnormals_;
  bool                   show_fnormals_;
  float                  normal_scale_;
  OpenMesh::FPropHandleT< typename Mesh::Point > fp_normal_base_;
};


//=============================================================================
#if defined(OM_INCLUDE_TEMPLATES) && !defined(OPENMESHAPPS_MESHVIEWERWIDGET_CC)
#  define OPENMESH_MESHVIEWERWIDGET_TEMPLATES
#  include "MeshPairViewerWidgetT.hpp"
#endif
//=============================================================================
#endif // OPENMESHAPPS_MESHVIEWERWIDGETT_HH defined
//=============================================================================


#endif // MESHPAIRVIEWERWIDGET_H
