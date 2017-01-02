#ifndef LAYER_H
#define LAYER_H
#include <memory>
#include "common.h"
class Layer
{
public:
    typedef std::shared_ptr<Layer> Ptr;
    explicit Layer();
    virtual bool config(Config::Ptr,size_t& n)=0;
    virtual bool init(
            std::vector<std::shared_ptr<arma::vec>>& in,
            std::vector<std::shared_ptr<arma::vec>>& out
            )=0;
    virtual double forward(
            std::vector<std::shared_ptr<arma::vec>>& in,
            std::vector<std::shared_ptr<arma::vec>>& out
            )=0;
    virtual double backward(
            std::vector<std::shared_ptr<arma::vec>>& out,
            std::vector<std::shared_ptr<arma::vec>>& dout,
            std::vector<std::shared_ptr<arma::vec>>& bp,
            std::vector<std::shared_ptr<arma::vec>>& in,
            std::vector<std::shared_ptr<arma::vec>>& din
            )=0;
protected:
    const std::string type_;
    const std::string name_;
};

#endif // LAYER_H
