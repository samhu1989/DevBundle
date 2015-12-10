#ifndef QGLPOINTSELECTION_H
#define QGLPOINTSELECTION_H
#include <armadillo>
#include <memory>
#include <vector>
class PointSelectionBase
{
public:
    typedef std::shared_ptr<PointSelectionBase> Ptr;
    typedef enum{
        Unknown,
        RayPoint,
        RayPoints,
        BoxPoints,
    }Type;
    PointSelectionBase():type_(Unknown){}
    PointSelectionBase(Type t):type_(t){}
    virtual Type type(){return type_;}
    virtual void debugSelection(){std::cerr<<"Debug Invalid Selection"<<std::endl;}
private:
    Type type_;
};

class PointSelections:public std::vector<PointSelectionBase::Ptr>
{
public:
    template<typename Mesh>
    inline void selectAll(Mesh&,arma::uvec&);
    template<typename Mesh>
    inline void selectAt(size_t,Mesh&,arma::uvec&);
    template<typename Mesh>
    inline void selectAt(PointSelections::iterator,Mesh&,arma::uvec&);
    void debugSelections();
};

class RayPointSelection:public PointSelectionBase
{
public:
    typedef std::shared_ptr<RayPointSelection> Ptr;
    RayPointSelection(arma::fvec& from,arma::fvec toward):
        PointSelectionBase(PointSelectionBase::RayPoint),from_(from),toward_(toward)
    {
    }
    template<typename Mesh>
    inline void select(Mesh&,arma::uvec&);
    virtual void debugSelection();
private:
    arma::fvec from_;
    arma::fvec toward_;
};
#include "qglpointselection.hpp"
#endif // QGLPOINTSELECTION_H
