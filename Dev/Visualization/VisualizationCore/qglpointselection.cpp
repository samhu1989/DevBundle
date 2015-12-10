#include "qglpointselection.h"
#include <QGL>
void PointSelections::debugSelections()
{
    for(iterator iter=begin();iter!=end();++iter)
    {
        PointSelectionBase::Ptr ptr = *iter;
        if(ptr&&0!=ptr.use_count())
        {
            ptr->debugSelection();
        }
    }
}

void RayPointSelection::debugSelection()
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glColor3f(0.1f, 0.9f, 0.9f); // greenish
    glVertex3fv(from_.memptr());
    glVertex3fv(toward_.memptr());
    glEnd();
}
