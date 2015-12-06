#ifndef OBJECTVIEW_H
#define OBJECTVIEW_H
#include <QFrame>
#include "objectmodel.h"
#include "MeshListViewerWidget.h"
namespace Ui {
class ObjectView;
}

class ObjectView : public QFrame
{
    Q_OBJECT
public:
    explicit ObjectView(
            std::vector<ObjModel::Ptr>& objects,
            QWidget *parent = 0
            );
    ~ObjectView();
protected slots:
    void view_color_weight(bool);
    void view_normal_weight(bool);
    void view_spatial_weight(bool);
    void view_original_color(bool);
private:
    std::vector<ObjModel::Ptr>& objects_;
    MeshListViewerWidget* widget_;
    Ui::ObjectView *ui;
};

#endif // OBJECTVIEW_H
