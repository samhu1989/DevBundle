#ifndef PATCHPAIRVIEW_H
#define PATCHPAIRVIEW_H

#include <QFrame>
#include "MeshViewerWidget.h"
namespace Ui {
class PatchPairView;
}

class PatchPairView : public QFrame
{
    Q_OBJECT

public:
    explicit PatchPairView(QWidget *parent = 0);
    ~PatchPairView();
    MeshViewerWidget* oview(){return object_view_;}
    MeshViewerWidget* pview(){return patch_view_;}
private:
    Ui::PatchPairView *ui;
    MeshViewerWidget* patch_view_;
    MeshViewerWidget* object_view_;
};

#endif // PATCHPAIRVIEW_H
