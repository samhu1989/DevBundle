#ifndef FEATUREVIEW_H
#define FEATUREVIEW_H

#include <QWidget>
#include "featureviewerwidget.h"
#include <QWheelEvent>
namespace Ui {
class featureview;
}

class featureview : public QWidget
{
    Q_OBJECT
public:
    explicit featureview(
            std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs,
            std::vector<arma::uvec>& labels,
            arma::mat& base,
            arma::mat& center,
            QWidget *parent = 0
            );
    ~featureview();
    bool configure(Config::Ptr);
    void init(void);
signals:
    void message(QString,int);
protected:
    void extract_patch_features();
    void set_patch_features();
private:
    Ui::featureview *ui;
    FeatureViewerWidget* viewer_;
    Config::Ptr config_;
    MeshBundle<DefaultMesh>::PtrList& inputs_;
    std::vector<arma::uvec>& labels_;
    arma::mat& feature_base_;
    arma::mat& feature_center_;
    std::vector<arma::urowvec> input_patch_label_value_;
    std::vector<arma::mat> patch_features_;
};

#endif // FEATUREVIEW_H
