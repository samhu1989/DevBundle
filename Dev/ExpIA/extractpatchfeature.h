#ifndef EXTRACTPATCHFEATURE_H
#define EXTRACTPATCHFEATURE_H
#include "common.h"
#include "colorhistogram.h"
#include "blockbasedfeature.h"
void extract_patch_feature(DefaultMesh&,arma::vec&,Config::Ptr);
#endif // EXTRACTPATCHFEATURE_H
