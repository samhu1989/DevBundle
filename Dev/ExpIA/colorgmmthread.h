#ifndef COLORGMMTHREAD_H
#define COLORGMMTHREAD_H
#include <QThread>
class ColorGMMThread:public QThread
{
    Q_OBJECT
public:
    explicit ColorGMMThread(QObject* parent=0);
    bool init();
protected:
    void run();
private:
};

#endif // COLORGMMTHREAD_H
