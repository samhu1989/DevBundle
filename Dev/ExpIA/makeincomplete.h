#ifndef MAKEINCOMPLETE_H
#define MAKEINCOMPLETE_H

#include <QThread>

class MakeIncomplete : public QThread
{
    Q_OBJECT
public:
    explicit MakeIncomplete(QObject *parent = 0);
signals:
public slots:
};

#endif // MAKEINCOMPLETE_H
