#ifndef GLOBALALIGN_H
#define GLOBALALIGN_H
#include "common.h"
#include <QFrame>
#include <QTimer>
#include "MeshListViewerWidget.h"
namespace Ui {
class GlobalAlign;
}

class GlobalAlign : public QFrame
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> IMeshList;
    explicit GlobalAlign(IMeshList& input,QWidget *parent = 0);
    ~GlobalAlign();
public:
    bool configure(Config::Ptr);
signals:
    void message(QString,int);
protected slots:
    void loadDownSampled(void);
    void start_AlignEachOther(void);
    void finish_AlignEachOther(void);
    void alignUptoZ(void);
private:
    MeshListViewerWidget* geo_view_;
    Ui::GlobalAlign *ui;
    IMeshList& inputs_;
    Config::Ptr config_;
    QTimer gl_timer;
    QThread* geo_thread_;
};

#endif // GLOBALALIGN_H
