#ifndef EXTRACTPATCHFEATURE_H
#define EXTRACTPATCHFEATURE_H
#include "common.h"
#include "colorhistogram.h"
#include "blockbasedfeature.h"
void extract_patch_feature(DefaultMesh&,arma::vec&,Config::Ptr);
void extract_patch_expand(DefaultMesh&i,arma::uvec&indices,DefaultMesh&o,int k);
void extract_patch_expand(DefaultMesh&i,arma::uvec&indices,DefaultMesh&o,float r);
#endif // EXTRACTPATCHFEATURE_H
