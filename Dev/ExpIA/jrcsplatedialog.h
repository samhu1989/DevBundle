#ifndef JRCSPLATEDIALOG_H
#define JRCSPLATEDIALOG_H

#include <QDialog>
#include "jrcsplatedialog.h"
#include "jrcsprimitive.h"
#include "MeshPairViewerWidget.h"
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
private:
    Ui::JRCSPlateDialog *ui;
    JRCS::Plate* plate;
    MeshPairViewerWidget* geo_view_;
};

#endif // JRCSPLATEDIALOG_H
