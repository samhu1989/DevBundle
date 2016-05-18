#include "hierarchicalization.h"
#include "pcaplaneequ.h"
Hierarchicalization::Hierarchicalization()
{

}


bool Hierarchicalization::configure(Config::Ptr config)
{
    if(config->has("JRCSInit_neighbor_radius"))
    {
        neighbor_radius_ = config->getFloat("JRCSInit_neighbor_radius");
    }else return false;
    if(config->has("JRCSInit_angle_tight"))
    {
        anglethres_tight_ = config->getFloat("JRCSInit_angle_tight");
    }else return false;
    if(config->has("JRCSInit_angle_relax"))
    {
        anglethres_relax_ = config->getFloat("JRCSInit_angle_relax");
    }else return false;
    if(config->has("JRCSInit_new_normal"))
    {
        if(config->getInt("JRCSInit_new_normal"))
        {
            force_new_normal_ = true;
        }else force_new_normal_ = false;
    }else force_new_normal_ = false;
    if(config->has("JRCSInit_point2plane"))
    {
        point2plane_th_ = config->getFloat("JRCSInit_point2plane");
    }else point2plane_th_ = 0.05;
    return true;
}

void Hierarchicalization::compute(DefaultMesh& m)
{
    std::cerr<<"reseting"<<std::endl;
    reset(m);
    std::cerr<<"computing neighbor"<<std::endl;
    calneighbor(m);
    std::cerr<<"building"<<std::endl;
    build(m);
}

void Hierarchicalization::getObjectLabel(arma::uvec&result)
{
    arma::ivec lbl = label_;
    lbl(lbl>=0) += 1;
    lbl(lbl<0).fill(0);
    std::vector<IdNode>::iterator iter;
    arma::uword id;
    id = 0;
    for(iter=idtree_.begin();iter!=idtree_.end();++iter)
    {
        if(!iter->is_obj_)
        {
            lbl(lbl==(id+1)).fill(0);
        }
        ++id;
    }
    result = arma::conv_to<arma::uvec>::from(lbl);
}

void Hierarchicalization::getPlaneLabel(arma::uvec&result)
{
    arma::ivec lbl = label_;
    lbl(lbl>=0) += 1;
    lbl(lbl<0).fill(0);
    std::vector<IdNode>::iterator iter;
    arma::uword id;
    id = 0;
    for(iter=idtree_.begin();iter!=idtree_.end();++iter)
    {
        if(iter->is_obj_)
        {
            lbl(lbl==(id+1)).fill(0);
        }
        ++id;
    }
    result = arma::conv_to<arma::uvec>::from(lbl);
}

void Hierarchicalization::calneighbor(DefaultMesh& mesh)
{
    MeshKDTreeInterface<DefaultMesh> tree_mesh_(mesh);
    KDTree tree(3,tree_mesh_,nanoflann::KDTreeSingleIndexAdaptorParams(3));
    tree.buildIndex();
    arma::uword index;
    float* pptr = (float*)mesh.points();
    std::vector<std::pair<arma::uword,float>> search_result;
    std::vector<std::pair<arma::uword,float>>::iterator riter;
    for(index=0;index<mesh.n_vertices();++index)
    {
        search_result.clear();
        tree.radiusSearch(pptr,neighbor_radius_,search_result,nanoflann::SearchParams(3));
        nei[index].index.clear();
        for(riter=search_result.begin();riter!=search_result.end();++riter)
        {
            nei[index].index.push_back(riter->first);
        }
        pptr += 3;
    }
    if(mesh.has_vertex_normals()&&!force_new_normal_)
    {
        float* nptr = (float*)mesh.vertex_normals();
        arma::fmat normals(nptr,3,mesh.n_vertices(),false,true);
        for(index=0;index<mesh.n_vertices();++index)
        {
            nei[index].normal = normals.col(index);
        }
    }else{
        arma::fvec target, temp;
        std::vector<arma::fvec> list;
        long s=0;
        pcaplaneequ plane;
        arma::fmat cloud((float*)mesh.points(),3,mesh.n_vertices(),false,true);
        for (long i=0; i<mesh.n_vertices(); i++)
            if (nei[i].index.size()>=4)
            {
                list.clear();
                target = cloud.col(i);
                for (long j=0; j<nei[i].index.size(); j++)
                {
                    temp = cloud.col(nei[i].index[j]);
                    list.push_back(temp);
                }
                for (long j=0; j<nei[i].index.size()-1; j++)
                {
                    s = j;
                    for (long k=j+1; k<nei[i].index.size(); k++)
                        if ( arma::norm( target - list[k] ) < arma::norm( target - list[s]) )
                            s=k;
                    if (s!=j)
                    {
                        temp = list[j];
                        list[j] = list[s];
                        list[s] = temp;
                    }
                }
                plane.clear();
                plane.push_point(target);
                for (long j=0; j<nei[i].index.size(); j++)
                    if (j<4)
                        plane.push_point(list[j]);
                    else
                        break;
                nei[i].normal = plane.getnormal();
            }
            else
                label_[i] = -1;
    }
}

void Hierarchicalization::reset(const DefaultMesh& mesh)
{
    label_ = arma::ivec(mesh.n_vertices(),arma::fill::zeros);
    nei.clear();
    nei.resize(mesh.n_vertices());
    //set root
    IdNode tempid;
    tempid.id_ = 0;
    tempid.is_obj_ = true;
    tempid.child_plane_.clear();
    tempid.child_obj_.clear();
    tempid.father_id_ = 0;
    idtree_.clear();
    idtree_.push_back(tempid);

}

void Hierarchicalization::getObjectBox(DefaultMesh& mesh)
{
    std::vector<IdNode>::iterator iter;
    arma::fmat boxes;
    for(iter=idtree_.begin()+1;iter!=idtree_.end();++iter)
    {
        if(iter->is_obj_)
        {
            buildBB(mesh);
            if(boxes.empty())boxes =  iter->bbox_.boxmat;
            else boxes = arma::join_rows(boxes,iter->bbox_.boxmat);
        }
    }
    arma::fmat meshboxmat((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    meshboxmat = boxes;
}

void Hierarchicalization::build(DefaultMesh& mesh)
{
    arma::fmat v((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    calsize(v,0);

    //width first search
    uint32_t listf, listr;
    uint32_t i,j;
    listf = 0;
    listr = 1;
    iden = 0;

    std::cerr<<"start build"<<std::endl;

    while (listf<listr)
    {
        //region grow idtree[listf], if idtree[listf] is obj
        float volumes;
        volumes = idtree_[listf].size_x*idtree_[listf].size_y*idtree_[listf].size_z;
        if ((idtree_[listf].is_obj_)&&(volumes>0.001))
        {
            regiongrow(v,idtree_[listf].id_);
        }
        listr = idtree_.size();
        listf++;
    }
}

void Hierarchicalization::calsize(const arma::fmat& cloud, const arma::uword cloudid)
{
    float x_min, x_max, y_min, y_max, z_min, z_max, temp_x, temp_y, temp_z;
    int k;
    for (int i=0; i<cloud.n_cols; i++)
        if (label_[i]==cloudid)
        {
            x_min = x_max = cloud(0,i);
            y_min = y_max = cloud(1,i);
            z_min = z_max = cloud(2,i);
            k = i+1;
            break;
        }
    for (int i=k; i<cloud.n_cols; i++)
        if (label_[i]==cloudid)
        {
            temp_x = cloud(0,i);
            temp_y = cloud(1,i);
            temp_z = cloud(2,i);
            if (temp_x<x_min) x_min = temp_x;
            if (temp_x>x_max) x_max = temp_x;
            if (temp_y<y_min) y_min = temp_y;
            if (temp_y>y_max) y_max = temp_y;
            if (temp_z<z_min) z_min = temp_z;
            if (temp_z>z_max) z_max = temp_z;
        }
    idtree_[cloudid].size_x = x_max-x_min;
    idtree_[cloudid].size_y = y_max-y_min;
    idtree_[cloudid].size_z = z_max-z_min;
    return;
}

float Hierarchicalization::calarea(float size_x, float size_y, float size_z)
{
    if (size_x<size_y)
        if(size_x<size_z)
            return size_y*size_z;
        else
            return size_y*size_x;
    else
        if (size_y<size_z)
            return size_x*size_z;
        else
            return size_x*size_y;
}

float Hierarchicalization::angle(const arma::fvec& v0 ,const arma::fvec& v1)
{
    float angle = std::fabs(arma::dot(v0,v1));
    angle /= (arma::norm(v0)*arma::norm(v1));
    return std::acos(angle)*180.0/M_PI;
}

BBox Hierarchicalization::calbbox(arma::fmat& pts)
{
    BBox result;
    arma::fmat box;
    get3DMBB(pts,2,box);
    result.boxmat = box;
    result.center = arma::mean(box,1);
    result.width = arma::norm(box.col(0)-box.col(1));
    result.height = arma::norm(box.col(0)-box.col(3));
    result.depth = arma::norm(box.col(0)-box.col(4));
    return result;
}


void Hierarchicalization::regiongrow(const arma::fmat& cloud,const arma::uword targetid)
{
    arma::uword num = cloud.n_cols;
    std::vector<arma::uword> list;
    std::vector<bool> t;	//if the point has been found
    arma::fvec temp_point;
    arma::fvec plane_normal;
    arma::fvec new_point;
    pcaplaneequ plane;
    arma::uword f, r, i, j;

    float x_min, x_max, y_min, y_max, z_min, z_max, size_x, size_y, size_z, area, plane_area;
    area = calarea(idtree_[targetid].size_x, idtree_[targetid].size_y, idtree_[targetid].size_z);

    t.clear();
    for (j=0; j<num; j++) t.push_back(true);

    for (i=0; i<num; i++)
        if ((label_[i]==targetid)&&(t[i]))
        {
            x_min = x_max = cloud(0,i);
            y_min = y_max = cloud(1,i);
            z_min = z_max = cloud(2,i);
            list.clear();
            plane.clear();
            f = 0;
            r = 1;
            list.push_back(i);
            t[i] = false;

            for (j=0; j<nei[i].index.size(); j++)
                if ((label_[nei[i].index[j]]==targetid)&&(t[nei[i].index[j]]))
                {
                    list.push_back(nei[i].index[j]);
                    temp_point = cloud.col(nei[i].index[j]);

                    if (temp_point(0)<x_min) x_min = temp_point(0);
                    if (temp_point(0)>x_max) x_max = temp_point(0);
                    if (temp_point(1)<y_min) y_min = temp_point(1);
                    if (temp_point(1)>y_max) y_max = temp_point(1);
                    if (temp_point(2)<z_min) z_min = temp_point(2);
                    if (temp_point(2)>z_max) z_max = temp_point(2);

                    plane.push_point(temp_point);
                    t[nei[i].index[j]] = false;
                    r++;
                }

            while (f<r)
            {
                for (j=0; j<nei[list[f]].index.size(); j++)
                {
                    temp_point = cloud.col(nei[list[f]].index[j]);

                    if ((label_[nei[list[f]].index[j]]==targetid)&&	//这个点没被找过
                            (t[nei[list[f]].index[j]])&&	// 这个点不在队列内
                            ( ((plane.getsize()<=6)&&(angle(plane.getnormal(), nei[nei[list[f]].index[j]].normal)<anglethres_tight_)) // 这个点的法向和平面法向相差anglethres度
                              ||((plane.getsize()>6)&&(plane.dist(temp_point)<point2plane_th_)&&(angle(plane.getnormal(), nei[nei[list[f]].index[j]].normal)<anglethres_relax_))))
                    {
                        if (temp_point(0)<x_min) x_min = temp_point(0);
                        if (temp_point(0)>x_max) x_max = temp_point(0);
                        if (temp_point(1)<y_min) y_min = temp_point(1);
                        if (temp_point(1)>y_max) y_max = temp_point(1);
                        if (temp_point(2)<z_min) z_min = temp_point(2);
                        if (temp_point(2)>z_max) z_max = temp_point(2);
                        list.push_back(nei[list[f]].index[j]);
                        plane.push_point(temp_point);
                        t[nei[list[f]].index[j]] = false;
                        r++;
                    }
                }
                f++;
            }

            //decide if it is a plane by size
            size_x = x_max-x_min;
            size_y = y_max-y_min;
            size_z = z_max-z_min;
            plane_area = calarea(size_x, size_y, size_z);
            if (plane_area>(area/7))
            {
                std::cerr<<"plane found"<<std::endl;
                //printf("list size = %d , plane_area = %5.3f , total_area = %5.3f \n", list.size(), plane_area, area);
                iden++;
                for (j=0; j<list.size(); j++)
                    label_[list[j]] = iden;
                IdNode tempnode;
                tempnode.id_ = iden;
                tempnode.size_x = size_x;
                tempnode.size_y = size_y;
                tempnode.size_z = size_z;
                tempnode.is_obj_ = false;
                tempnode.child_obj_.clear();
                tempnode.child_plane_.clear();
                tempnode.father_id_ = targetid;
                idtree_.push_back(tempnode);
                idtree_[targetid].child_plane_.push_back(iden);
            }
        }

    t.clear();
    for (j=0; j<num; j++) t.push_back(true);
    for (i=0;i<num;i++)
        if ((label_[i]==targetid)&&t[i])
        {
            x_min = x_max = cloud(0,i);
            y_min = y_max = cloud(1,i);
            z_min = z_max = cloud(2,i);
            list.clear();
            f = 0;
            r = 1;
            list.push_back(i);
            while (f<r)
            {
                for (j=0; j<nei[list[f]].index.size(); j++)
                    if ((label_[nei[list[f]].index[j]]==targetid)&&
                            (t[nei[list[f]].index[j]]))
                    {
                        temp_point = cloud.col(nei[list[f]].index[j]);

                        if (temp_point(0)<x_min) x_min = temp_point(0);
                        if (temp_point(0)>x_max) x_max = temp_point(0);
                        if (temp_point(1)<y_min) y_min = temp_point(1);
                        if (temp_point(1)>y_max) y_max = temp_point(1);
                        if (temp_point(2)<z_min) z_min = temp_point(2);
                        if (temp_point(2)>z_max) z_max = temp_point(2);

                        list.push_back(nei[list[f]].index[j]);
                        t[nei[list[f]].index[j]]=false;
                        r++;

                    }
                f++;
            }
            size_x = x_max-x_min;
            size_y = y_max-y_min;
            size_z = z_max-z_min;

            //printf("volcompare = %5.4lf\n", volcompare);

            arma::fmat pts; //find minimal bounding box
            arma::uvec pts_index(list);
            pts = cloud.cols(pts_index);
            BBox temp_box;
            temp_box = calbbox(pts);

            float volcompare;
            volcompare = temp_box.width * temp_box.height * temp_box.depth/(idtree_[targetid].size_x*idtree_[targetid].size_y*idtree_[targetid].size_z);
            if ((list.size()>100)&&(volcompare<0.4))
            {
                std::cerr<<"obj found"<<std::endl;
                iden++;
                for (j=0; j<list.size(); j++)
                    label_[list[j]] = iden;
                IdNode tempnode;
                tempnode.id_ = iden;
                tempnode.size_x = size_x;
                tempnode.size_y = size_y;
                tempnode.size_z = size_z;
                tempnode.is_obj_ = true;
                tempnode.child_obj_.clear();
                tempnode.child_plane_.clear();
                tempnode.father_id_ = targetid;
                assert(temp_box.boxmat.n_cols==8);
                tempnode.bbox_ = temp_box;
                idtree_.push_back(tempnode);
                idtree_[targetid].child_obj_.push_back(tempnode.id_);
            }
        }
    std::cerr<<"End region grow"<<std::endl;
}
