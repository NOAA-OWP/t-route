#include <stdlib.h>
#include "reach_structs.h"

void init_reach(_Reach* reach, long* upstream_ids, int num_upstream_ids, int reach_type)
{
  int i;
  if( reach != NULL){
    reach->type = reach_type;
    reach->_num_upstream_ids = num_upstream_ids;
    if(num_upstream_ids > 0 && upstream_ids != NULL){
      reach->_upstream_ids = (long*) malloc(sizeof(long)*(num_upstream_ids));
      for(i = 0; i < num_upstream_ids; i++){
        reach->_upstream_ids[i] = upstream_ids[i];
      }
    }
    else {
      reach->_upstream_ids = (long*) malloc(sizeof(long));
      reach->_upstream_ids[0] = -1;
    }
  }
}

void free_reach(_Reach* reach)
{
  if(reach != NULL && reach->_upstream_ids != NULL)
    free(reach->_upstream_ids);
}

long* get_upstream_ids(_Reach* reach)
{
  return reach->_upstream_ids;
}
