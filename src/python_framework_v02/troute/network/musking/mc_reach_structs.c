#include <stdlib.h>
#include "../reach_structs.h"

void init_mc_reach(_Reach* reach, int num_segments)
{
  reach->reach.mc_reach.num_segments = num_segments;
  reach->reach.mc_reach._segments = (_MC_Segment*) malloc(sizeof(_MC_Segment)*(num_segments));
}

void free_mc_reach(_Reach* reach)
{
  if( reach != NULL && reach->reach.mc_reach._segments != NULL)
    free(reach->reach.mc_reach._segments);
}

void set_mc_segment(_Reach* reach, int index, long id,
    float dt, float dx, float bw, float tw, float twcc,
    float n, float ncc, float cs, float s0,
    float qdp, float velp, float depthp)
{
  if(index > -1 && index < reach->reach.mc_reach.num_segments)
  {
    _MC_Segment segment;
    //TODO what about segment id???
    segment.id = id;
    segment.dt = dt;
    segment.dx = dx;
    segment.bw = bw;
    segment.tw = tw;
    segment.twcc = twcc;
    segment.n = n;
    segment.ncc = ncc;
    segment.cs = cs;
    segment.s0 = s0;
    segment.qdp = qdp;
    segment.velp = velp;
    segment.depthp = depthp;
    reach->reach.mc_reach._segments[index] = segment;
  }
  //FIXME else what?
}

_MC_Segment get_mc_segment(_Reach* reach, int index)
{
  _MC_Segment seg;
  if( reach != NULL && index > -1 && index < reach->reach.mc_reach.num_segments){
    if(reach->reach.mc_reach._segments != NULL )
      seg = reach->reach.mc_reach._segments[index];
  }
  return seg;
}
