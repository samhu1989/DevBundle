#ifndef MESHLABELVIEWERWIDGET_H
#define MESHLABELVIEWERWIDGET_H

#include <QWidget>
#include "MeshLabelViewerWidgetPrivate.h"
namespace Ui {
class MeshLabelViewerWidget;
}

class MeshLabelViewerWidget : public QWidget
{
    Q_OBJECT

public:
    explicit MeshLabelViewerWidget(QWidget *parent = 0);
    ~MeshLabelViewerWidget();

private:
    Ui::MeshLabelViewerWidget *ui;
    MeshLabelViewerWidgetPrivate* view;
};

#endif // MESHLABELVIEWERWIDGET_H
