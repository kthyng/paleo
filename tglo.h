/*
** svn $Id: upwelling.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2014 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for TGLO sim.
**
** Application flag:   TGLO
** Input script:       ocean_tglo.in
*/

#define ROMS_MODEL

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define VISC_GRID


#define TS_MPDATA
#define TS_DIF2
#define MIX_GEO_TS
#define SALINITY
#define DIFF_GRID

#define MASKING
#define CURVGRID
#define SOLVE3D
#define SPLINES
#define DJ_GRADPS

#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
