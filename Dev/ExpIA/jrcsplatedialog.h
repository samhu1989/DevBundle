#ifndef JRCSPLATEDIALOG_H
#define JRCSPLATEDIALOG_H

#include <QDialog>
#include "jrcsplatedialog.h"
#include "jrcsprimitive.h"
#include "MeshPairViewerWidget.h"
#include <QTimer>
namespace Ui {
class JRCSPlateDialog;
}

class JRCSPlateDialog : public QDialog
{
    Q_OBJECT

public:
    explicit JRCSPlateDialog(QWidget *parent = 0);
    ~JRCSPlateDialog();
    void init_plate();
    void init_points();
protected slots:
    void start_transform();
    void transform();
    void start_fit();
    void fit();
private:
    Ui::JRCSPlateDialog *ui;
    JRCS::Plate* plate;
    MeshPairViewerWidget* geo_view_;
    arma::fmat dR_;
    arma::fmat translate_;
    float dt_;
    float time_;
    QTimer* timer_;
};

#endif // JRCSPLATEDIALOG_H
