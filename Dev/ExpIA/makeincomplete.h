#ifndef MAKEINCOMPLETE_H
#define MAKEINCOMPLETE_H
#include "common.h"
#include <QThread>
#include <QString>
#include <QList>
class MakeIncomplete : public QThread
{
    Q_OBJECT
public:
    explicit MakeIncomplete(
            MeshBundle<DefaultMesh>::PtrList inputs,
            QList<float> ratio,
            QString path,
            QObject *parent = 0
            );
signals:
public slots:
protected slots:
    void run(void);
protected:
    void saveByRatio(MeshBundle<DefaultMesh>::Ptr in);
private:
    MeshBundle<DefaultMesh>::PtrList inputs_;
    QList<float> ratio_;
    QString path_;
};

#endif // MAKEINCOMPLETE_H
