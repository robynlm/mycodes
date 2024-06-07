/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ProperTime_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ProperTime";
  CCTK_RegisterBanner(banner);
  return 0;
}
