#ifndef KDTREE
#define KDTREE
template<typename M>
struct MeshKDTreeInterface{
    MeshKDTreeInterface(const M&mesh):mesh_(mesh){}
    inline size_t kdtree_get_point_count() const { return mesh_.n_vertices(); }
    inline float kdtree_distance(const float *p1, const size_t idx_p2,size_t /*size*/) const
    {
        const float* p2 = (const float*)(mesh_.points());
        const float d0=p1[0]-p2[3*idx_p2];
        const float d1=p1[1]-p2[3*idx_p2+1];
        const float d2=p1[2]-p2[3*idx_p2+2];
        return d0*d0+d1*d1+d2*d2;
    }
    inline float kdtree_get_pt(const size_t idx, int dim) const
    {
        const float* p2 = (const float*)(mesh_.points());
        return p2[3*idx+dim];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&)const{return false;}
    const M& mesh_;
};
#endif // KDTREE

