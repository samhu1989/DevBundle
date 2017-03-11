#ifndef JRCSPLATEDIALOG_H
#define JRCSPLATEDIALOG_H

#include <QDialog>
#include "jrcsprimitive.h"
#include "jrcscube.h"
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
    void init_for_plate();
    void init_for_cube();
protected:
    void init_plate();
    void init_cube();
    void init_points_for_plate();
    void init_points_for_cube();
    void sample_points_for_cube(JRCS::Cube *cube);
protected slots:
    void start_transform_plate();
    void transform_plate();
    void start_fit_plate();
    void fit_plate();

    void start_transform_cube();
    void transform_cube();
private:
    Ui::JRCSPlateDialog *ui;
    JRCS::Plate* plate;
    JRCS::Cube* cube;
    JRCS::Cube* cube2;
    MeshPairViewerWidget* geo_view_;
    arma::fmat dR_;
    arma::fmat translate_;
    float dt_;
    float time_;
    QTimer* timer_;
    static const int pN_ = 3000;
};

#endif // JRCSPLATEDIALOG_H
