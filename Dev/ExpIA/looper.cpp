#include "looper.h"
#include <QDebug>
bool Looper::configure(Config::Ptr config)
{
    config_= config;
    return true;
}

void Looper::loop()
{
    qDebug("Looping");
    emit end();
}

