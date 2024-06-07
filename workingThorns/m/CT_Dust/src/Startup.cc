/*  File produced by Kranc */

#include "cctk.h"

extern "C" int CT_Dust_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "CT_Dust";
  CCTK_RegisterBanner(banner);
  return 0;
}
