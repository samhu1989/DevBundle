#include "unifylabelcolorsizethread.h"

void UnifyLabelColorSizeThread::run()
{

}

bool UnifyLabelColorSizeThread::configure(Config::Ptr config)
{
    config_ = config;
    if(!config_->has("NormalEstimation_k"))return false;
    if(!config_->has("RegionGrow_k"))return false;
    if(!config_->has("RegionGrow_cluster_min_num"))return false;
    if(!config_->has("RegionGrow_cluster_max_num"))return false;
    if(!config_->has("RegionGrow_max_norm_angle"))return false;
    return true;
}

