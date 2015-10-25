#ifndef MESHPAIRVIEWERWIDGET_H
#define MESHPAIRVIEWERWIDGET_H
#include "MeshType.h"
//== INCLUDES =================================================================
#include <iostream>
#include <visualizationcore_global.h>
#include <QWidget>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
#include <OpenMesh/Tools/Utils/getopt.h>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include "MeshType.h"
#include "MeshPairViewerWidgetT.h"


//== CLASS DEFINITION =========================================================
class VISUALIZATIONCORESHARED_EXPORT MeshPairViewerWidget : public MeshPairViewerWidgetT<DefaultMesh>
{
    Q_OBJECT
public:
    using MeshPairViewerWidgetT<DefaultMesh>::Mesh;
    /// default constructor
    MeshPairViewerWidget(QWidget* parent=0) : MeshPairViewerWidgetT<DefaultMesh>(parent)
    {}
    OpenMesh::IO::Options& options() { return _options; }
    const OpenMesh::IO::Options& options() const { return _options; }
    void setOptions(const OpenMesh::IO::Options& opts) { _options = opts; }

    void open_mesh_gui(QString fname,MeshBundle<Mesh>&bundle)
    {
        OpenMesh::Utils::Timer t;
        t.start();
        if ( fname.isEmpty() || !open_mesh( fname.toStdString().c_str(),bundle.mesh_,bundle.strips_,_options) )
        {
            QString msg = "Cannot read mesh from file:\n '";
            msg += fname;
            msg += "'";
            QMessageBox::critical( NULL, windowTitle(), msg);
        }
        t.stop();
        std::cout << "Loaded mesh in ~" << t.as_string() << std::endl;
    }

    void open_texture_gui(QString fname)
    {
        if ( fname.isEmpty() || !open_texture( fname.toStdString().c_str() ) )
        {
            QString msg = "Cannot load texture image from file:\n '";
            msg += fname;
            msg += "'\n\nPossible reasons:\n";
            msg += "- Mesh file didn't provide texture coordinates\n";
            msg += "- Texture file does not exist\n";
            msg += "- Texture file is not accessible.\n";
            QMessageBox::warning( NULL, windowTitle(), msg );
        }
    }

public slots:
    void query_open_source_file() {
        QString fileName = QFileDialog::getOpenFileName(this,
            tr("Open mesh file"),
            tr("../../Dev_Data/"),
            tr("OBJ Files (*.obj);;"
            "OFF Files (*.off);;"
            "STL Files (*.stl);;"
            "PLY Files (*.ply);;"
            "All Files (*)"));
        if (!fileName.isEmpty())
            open_mesh_gui(fileName,*first_);
    }
    void query_open_target_file() {
        QString fileName = QFileDialog::getOpenFileName(this,
            tr("Open mesh file"),
            tr("../../Dev_Data/"),
            tr("OBJ Files (*.obj);;"
            "OFF Files (*.off);;"
            "STL Files (*.stl);;"
            "PLY Files (*.ply);;"
            "All Files (*)"));
        if (!fileName.isEmpty())
            open_mesh_gui(fileName,*second_);
    }
    void query_open_texture_file() {
        QString fileName = QFileDialog::getOpenFileName(this,
            tr("Open texture file"),
            tr("../../Dev_Data/"),
            tr("PNG Files (*.png);;"
            "BMP Files (*.bmp);;"
            "GIF Files (*.gif);;"
            "JPEG Files (*.jpg);;"
            "TIFF Files (*.tif);;"
            "All Files (*)"));
        if (!fileName.isEmpty())
            open_texture_gui(fileName);
    }
private:
    OpenMesh::IO::Options _options;
};

#endif // MESHPAIRVIEWERWIDGET

