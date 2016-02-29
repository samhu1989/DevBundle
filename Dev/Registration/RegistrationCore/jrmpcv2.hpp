#include "jrmpcv2.h"
namespace Registration {
template<typename M>
bool JRMPCV2<M>::configure(Config::Ptr& config,InfoPtr& info)
{
    info = std::make_shared<Info>();
    if( !config || 0==config.use_count() )
    {
        std::cerr<<"no configure"<<std::endl;
        return false;
    }else{
        if(config->has("Align_Max_Iter")){
            info->max_iter = config->getInt("Align_Max_Iter");
        }
        if(config->has("Align_Eps"))
        {
            info->eps = config->getFloat("Align_Eps");
        }
        return true;
    }
    return false;
}
}

