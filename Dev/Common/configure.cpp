#include "configure.h"
#include <QString>
#include <QTextStream>
#include <QFile>
#include <iostream>
Config::Config(const std::string& path)
{
    reload(path);
}

void Config::reload()
{
    QFile inFile(QString::fromStdString(_SourcePath));
    inFile.open(inFile.ReadOnly);
    QTextStream stream(&inFile);
    QString str0;
    QString str1;
    while(!stream.atEnd())
    {
        stream >> str0;
        if(str0.startsWith("#"))
        {
            continue;
        }
        stream >> str1;
        _Config.insert(str0.toStdString(),str1.toStdString());
    }
}

void Config::reload(const std::string& path)
{
    _SourcePath = path;
    reload();
}

void Config::add(const std::string& key,const std::string& value)
{
    _Config.insert(key,value);
}

bool Config::has(const std::string& key)
{
    if(_Config.contains(key))
    {
        return true;
    }else{
        std::cerr<<"Missing Config: "<<key;
    }
    return false;
}

std::string Config::getString(const std::string& key)
{
    return _Config.value(key);
}

int Config::getInt(const std::string& key)
{
    QString str = QString::fromStdString(_Config.value(key));
    bool ok;
    int result = str.toInt(&ok);
    if( !ok )std::cerr<<"Can't convert "+key+" to Int";
    return result;
}

float Config::getFloat(const std::string& key)
{
    QString str = QString::fromStdString(_Config.value(key));
    bool ok;
    float result = str.toFloat(&ok);
    if( !ok )std::cerr<<"Can't convert "+key+" to Float";
    return result;
}

double  Config::getDouble(const std::string& key)
{
    QString str = QString::fromStdString(_Config.value(key));
    bool ok;
    double result = str.toDouble(&ok);
    if( !ok )std::cerr<<"Can't convert "+key+" to Double";
    return result;
}

Config::~Config()
{

}

