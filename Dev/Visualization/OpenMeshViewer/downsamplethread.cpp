#include "downsamplethread.h"
#include "filter.h"
void DownSampleThread::run()
{
    Filter::OctreeGrid<DefaultMesh> filter;
    filter.set_seed_resolution(0.1);
    filter.extract(m_);
}
