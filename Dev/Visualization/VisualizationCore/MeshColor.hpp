#include "MeshColor.h"
template<typename M>
MeshColor<M>::MeshColor(const Mesh& m):Ref_(m)
{

}
template<typename M>
ColorArray::RGBArray MeshColor<M>::vertex_colors_array(void)
{
    return v_colors;
}

template<typename M>
void* MeshColor<M>::vertex_colors(void)
{
    if(!v_colors.data_)
    {
        v_colors.reset(Ref_.n_vertices(),250,218,141);
    }else if(Ref_.n_vertices()!=v_colors.size_)
    {
        v_colors.reset(Ref_.n_vertices(),250,218,141);
    }
    return (void*)v_colors.data_;
}

