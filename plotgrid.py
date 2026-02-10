"""
$Id. Plotgrid mimic plotgrid.ncl script for checking domain created by WPS namelist setting.

 Program ini membutuhkan paket f90nml

namelist /share/ wrf_core, max_dom, start_date, end_date, &
                     start_year, end_year, start_month, end_month, &
                     start_day, end_day, start_hour, end_hour, &
                     interval_seconds, io_form_geogrid, opt_output_from_geogrid_path, &
                     debug_level, active_grid
   namelist /geogrid/ parent_id, parent_grid_ratio, &
                      i_parent_start, j_parent_start, s_we, e_we, s_sn, e_sn, &
                      map_proj, ref_x, ref_y, ref_lat, ref_lon, &
                      truelat1, truelat2, stand_lon, dx, dy, pole_lat, pole_lon, &
                      geog_data_res, geog_data_path, opt_geogrid_tbl_path

   !Initialize array

   ref_x = -999.0
   ref_y = -999.0
   ref_lat = -999.0
   ref_lon = -999.0
   do i=1,NVAR
    do j=1,MAX_DOMAINS
      plotvar(i,j) = -999.0
    end do
   end do

   ! Read parameters from Fortran namelist
   do funit=10,100
      inquire(unit=funit, opened=is_used)
      if (.not. is_used) exit
   end do
   open(funit,file=fname,status='old',form='formatted')
   read(funit,share)
   read(funit,geogrid)
   close(funit)

   !Assign integers to map projections

   if (index(map_proj, 'lambert') /= 0) then
      mproj_int = 1
   else if (index(map_proj, 'mercator') /= 0) then
      mproj_int = 2
   else if (index(map_proj, 'polar') /= 0) then
      mproj_int = 3
   else if (index(map_proj, 'lat-lon') /= 0) then
      mproj_int = 4
   end if

   !Put all the variables into an array
   plotvar(1,1) = max_dom
   plotvar(2,1) = dx
   plotvar(3,1) = dy
   plotvar(4,1) = ref_lat
   plotvar(5,1) = ref_lon
   plotvar(6,1) = ref_x
   plotvar(7,1) = ref_y
   plotvar(8,1) = truelat1
   plotvar(9,1) = truelat2
   plotvar(10,1) = stand_lon
   plotvar(11,1) = mproj_int
   plotvar(12,1) = pole_lat
   plotvar(13,1) = pole_lon
   do j=1,int(max_dom)
     plotvar(14,j) = parent_id(j)
     plotvar(15,j) = parent_grid_ratio(j)
     plotvar(16,j) = i_parent_start(j)
     plotvar(17,j) = j_parent_start(j)
     plotvar(18,j) = e_we(j)
     plotvar(19,j) = e_sn(j)
   end do

Author: jpatiani@itb.ac.id
Date: 20220418
"""

import f90nml
from wrf import xy_to_ll_proj
import cartopy.crs as crs
import cartopy.feature as cfeat
import matplotlib.pyplot as p
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import sys

p.rcParams.update({'font.size':16})
filename    = sys.argv[1] #"namelist.wps.all_options"
nml         = f90nml.read(filename)

var         = ('max_dom','dx','dy','ref_lat','ref_lon','ref_x','ref_y',
        'truelat1','truelat2','stand_lon','mproj_int','pole_lat','pole_lon',
        'parent_id','parent_grid_ratio','i_parent_start','j_parent_start',
        'e_we','e_sn')

pvar        = {}
for i in var:
    if i in nml['share']:
        pvar[i] = nml['share'][i]
    else:
        if i in nml['geogrid']:
            pvar[i] = nml['geogrid'][i]
        else:
            pvar[i] = -999.0

## Change map projection
map_proj    = {"lambert": 1, "polar": 2, "mercator": 3, "lat-lon": 6}
pvar["mproj_int"] = map_proj[nml["geogrid"]["map_proj"]]


## Kebutuhan res
## DX, DY, LATINC, LONINC, MAP_PROJ, TRUELAT1, TRUELAT2, STAND_LON, REF_LAT
## REF_LON, REF_X (KNOWNI), REF_Y (KNOWNJ), POLE_LAT, POLE_LON, 
pvar["latinc"]  = 0.0
pvar["loninc"]  = 0.0
e_we            = float(pvar["e_we"][0]) if type(pvar["e_we"]) is list else float(pvar["e_we"])
e_sn            = float(pvar["e_sn"][0]) if type(pvar["e_sn"]) is list else float(pvar["e_sn"])
pvar["known_x"] = pvar["ref_x"] if pvar["ref_x"] != -999.0 else e_we/2.
pvar["known_y"] = pvar["ref_y"] if pvar["ref_y"] != -999.0 else e_sn/2.
pvar["pole_lat"] = pvar["pole_lat"] if pvar["pole_lat"] != -999.0 else 90.0
pvar["pole_lon"] = pvar["pole_lon"] if pvar["pole_lon"] != -999.0 else 0.0

# Deal with global wrf domain
if pvar["dx"] < 1e-10 and pvar["dy"] < 1e-10:
    dx = 360./(e_we - 1)
    dy = 180./(e_sn - 1)
    pvar["ref_lat"] = 0.0
    pvar["ref_lon"] = 180.0

grid_to_plot= 1
adjust_grid = 0.5 if grid_to_plot == 1 else 0.0

xx          = 1.0 - adjust_grid
yy          = 1.0 - adjust_grid

loc         = xy_to_ll_proj(xx,yy,map_proj=pvar["mproj_int"],
        truelat1=pvar["truelat1"],truelat2=pvar["truelat2"],
        stand_lon=pvar["stand_lon"],ref_lat=pvar["ref_lat"],
        ref_lon=pvar["ref_lon"],known_x=pvar["known_x"],
        known_y=pvar["known_y"],pole_lat=pvar["pole_lat"],
        pole_lon=pvar["pole_lon"],dx=pvar["dx"],dy=pvar["dy"],
        latinc=pvar["latinc"],loninc=pvar["loninc"])
start_lat   = loc[0].values
start_lon   = loc[1].values

xx          = (e_we-1) + adjust_grid
yy          = (e_sn-1) + adjust_grid

loc         = xy_to_ll_proj(xx,yy,map_proj=pvar["mproj_int"],
        truelat1=pvar["truelat1"],truelat2=pvar["truelat2"],
        stand_lon=pvar["stand_lon"],ref_lat=pvar["ref_lat"],
        ref_lon=pvar["ref_lon"],known_x=pvar["known_x"],
        known_y=pvar["known_y"],pole_lat=pvar["pole_lat"],
        pole_lon=pvar["pole_lon"],dx=pvar["dx"],dy=pvar["dy"],
        latinc=pvar["latinc"],loninc=pvar["loninc"])
end_lat     = loc[0].values
end_lon     = loc[1].values

## Plot
if pvar["mproj_int"] == 1:
    globe = crs.Globe(ellipse='sphere', semimajor_axis=6370000, semiminor_axis=6370000)
    proj = crs.LambertConformal(central_longitude=pvar["stand_lon"],
            standard_parallels=(pvar["truelat1"], pvar["truelat2"]),
            globe=None)
    pgrid= True
    # mpLambertParallel1F(truelat1) mpLambertParallel2F(truelat2) 
    # mpLambertMeridianF(stand_lon)
elif pvar["mproj_int"] == 2:
    proj = crs.Stereographic(central_longitude=pvar["stand_lon"],
            central_latitude=pvar["ref_lat"])
    pgrid= False
    # mpCenterLatF(cenlat)/mpCenterLatF(ref_lat)
    # mpCenterLonF(stand_lon)
elif pvar["mproj_int"] == 3:
    proj = crs.Mercator(central_longitude=pvar["stand_lon"])
    pgrid= False
    # mpCenterLatF(0.0) mpCenterLonF(stand_lon)
else:
    pgrid= False
    if "pole_lat" in pvar and "pole_lon" in pvar and "stand_lon" in pvar:
        if pvar["pole_lon"] == 0 and pvar["pole_lat"] == 90:
            # not rotated
            proj = crs.PlateCarree(central_longitude=180-pvar["stand_lon"])
            # mpCenterLatF(0.0) mpCenterLonF(180 - stand_lon)
        else:
            # rotated
            southern = False
            if pvar["pole_lon"] == 0.0:
                southern = True
            elif pvar["pole_lon"] != 180:
                southern = True

            if not southern:
                proj = crs.PlateCarree(central_longitude=-pvar["stand_lon"])
            else:
                proj = crs.PlateCarree(central_longitude=180-pvar["stand_lon"])
    else:
        proj = crs.PlateCarree(central_longitude=180)
                        

fig     = p.figure()
ax      = fig.add_subplot(1,1,1, projection=proj)
ax.add_feature(cfeat.LAND)
ax.add_feature(cfeat.OCEAN)
ax.add_feature(cfeat.LAKES)
ax.add_feature(cfeat.COASTLINE)
ax.add_feature(cfeat.STATES)
ax.add_feature(cfeat.BORDERS)
if pvar["mproj_int"] != 6:
    ax.set_extent((start_lon,end_lon,start_lat,end_lat))

# Tick label major
gl      = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True,
        linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlines = pgrid
gl.ylines = pgrid
#gl.xlocator = mticker.FixedLocator()
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlable_style = {'size': 26, 'color': 'gray'}

# Draw text Dom01
ax.text(start_lon,end_lat,"d01",va="top",fontsize=16,
        color="black",transform=crs.PlateCarree())

# Dom Line
domLineColors   = ("white","red","green","blue")

if pvar["max_dom"] > 1:
    numLineColors   = 0
    if type(domLineColors) is tuple:
        numLineColors = len(domLineColors)

    for idom in range(1,pvar["max_dom"]):
        if numLineColors > 0:
            if numLineColors >= idom:
                gsLineColor = domLineColors[idom-1]
                txFontColor = domLineColors[idom-1]
            else:
                gsLineColor = domLineColors[numLineColors-1]
                txFontColor = domLineColors[numLineColors-1]

        if pvar["parent_id"][idom] == 1:
            goffset = 0.5*(1-(1.0/pvar["parent_grid_ratio"][idom]))
            i_start = pvar["i_parent_start"][idom]-goffset
            j_start = pvar["j_parent_start"][idom]-goffset
            # end point
            # change to mass point
            i_end   = (pvar["e_we"][idom]-2)/(1.0*pvar["parent_grid_ratio"][idom]) + i_start
            j_end   = (pvar["e_sn"][idom]-2)/(1.0*pvar["parent_grid_ratio"][idom]) + j_start

            if grid_to_plot == 0:
                adjust_grid = 0.0
            elif grid_to_plot == 1:
                adjust_grid = 0.5/(1.0*pvar["parent_grid_ratio"][idom])
            else:
                print(f"Error: Invalid value for grid_to_plot = {grid_to_plot}")
                adjust_grid = 0.0

        if pvar["parent_id"][idom] == 2:
            nd      = pvar["parent_id"][idom]
            i_points= ((pvar["e_we"][idom]-2)/(1.0*pvar["parent_grid_ratio"][idom]))
            j_points= ((pvar["e_sn"][idom]-2)/(1.0*pvar["parent_grid_ratio"][idom]))
            goffset = 0.5*(1-(1.0/(1.0*pvar["parent_grid_ratio"][idom])))
            ai_start= pvar["i_parent_start"][idom]*1.0-goffset
            aj_start= pvar["j_parent_start"][idom]*1.0-goffset

            if grid_to_plot == 0:
                adjust_grid = 0.0
            elif grid_to_plot == 1:
                adjust_grid = 0.5/(1.0*pvar["parent_grid_ratio"][idom])
            else:
                print(f"Error: Invalid value for grid_to_plot = {grid_to_plot}")
                adjust_grid = 0.0

            while nd > 1:
                goffset = 0.5*(1-(1.0/(1.0*pvar["parent_grid_ratio"][nd-1])))
                ai_start= (ai_start-1)/(1.0*pvar["parent_grid_ratio"][nd-1]) + pvar["i_parent_start"][nd-1]-goffset
                aj_start= (aj_start-1)/(1.0*pvar["parent_grid_ratio"][nd-1]) + pvar["j_parent_start"][nd-1]-goffset
                i_points= (i_points/(1.0*pvar["parent_grid_ratio"][nd-1]))
                j_points= (j_points/(1.0*pvar["parent_grid_ratio"][nd-1]))
                if grid_to_plot == 0:
                    adjust_grid = 0.0
                elif grid_to_plot == 1:
                    adjust_grid = adjust_grid/(1.0*pvar["parent_grid_ratio"][nd-1])
                else:
                    print(f"Error: Invalid value for grid_to_plot = {grid_to_plot}")
                    adjust_grid = 0.0
                nd = pvar["parent_id"][nd-1]

            i_start = ai_start
            j_start = aj_start
            i_end   = i_points + i_start
            j_end   = j_points + j_start

        # get the four corners
        xx  = i_start - adjust_grid
        yy  = j_start - adjust_grid
        loc = xy_to_ll_proj(xx,yy,map_proj=pvar["mproj_int"],
                truelat1=pvar["truelat1"],truelat2=pvar["truelat2"],
                stand_lon=pvar["stand_lon"],ref_lat=pvar["ref_lat"],
                ref_lon=pvar["ref_lon"],known_x=pvar["known_x"],
                known_y=pvar["known_y"],pole_lat=pvar["pole_lat"],
                pole_lon=pvar["pole_lon"],dx=pvar["dx"],dy=pvar["dy"],
                latinc=pvar["latinc"],loninc=pvar["loninc"])
        lat_SW = loc[0].values
        lon_SW = loc[1].values

        xx  = i_end + adjust_grid
        yy  = j_start - adjust_grid
        loc = xy_to_ll_proj(xx,yy,map_proj=pvar["mproj_int"],
                truelat1=pvar["truelat1"],truelat2=pvar["truelat2"],
                stand_lon=pvar["stand_lon"],ref_lat=pvar["ref_lat"],
                ref_lon=pvar["ref_lon"],known_x=pvar["known_x"],
                known_y=pvar["known_y"],pole_lat=pvar["pole_lat"],
                pole_lon=pvar["pole_lon"],dx=pvar["dx"],dy=pvar["dy"],
                latinc=pvar["latinc"],loninc=pvar["loninc"])
        lat_SE = loc[0].values
        lon_SE = loc[1].values

        xx  = i_start - adjust_grid
        yy  = j_end + adjust_grid
        loc = xy_to_ll_proj(xx,yy,map_proj=pvar["mproj_int"],
                truelat1=pvar["truelat1"],truelat2=pvar["truelat2"],
                stand_lon=pvar["stand_lon"],ref_lat=pvar["ref_lat"],
                ref_lon=pvar["ref_lon"],known_x=pvar["known_x"],
                known_y=pvar["known_y"],pole_lat=pvar["pole_lat"],
                pole_lon=pvar["pole_lon"],dx=pvar["dx"],dy=pvar["dy"],
                latinc=pvar["latinc"],loninc=pvar["loninc"])
        lat_NW = loc[0].values
        lon_NW = loc[1].values

        xx  = i_end + adjust_grid
        yy  = j_end + adjust_grid
        loc = xy_to_ll_proj(xx,yy,map_proj=pvar["mproj_int"],
                truelat1=pvar["truelat1"],truelat2=pvar["truelat2"],
                stand_lon=pvar["stand_lon"],ref_lat=pvar["ref_lat"],
                ref_lon=pvar["ref_lon"],known_x=pvar["known_x"],
                known_y=pvar["known_y"],pole_lat=pvar["pole_lat"],
                pole_lon=pvar["pole_lon"],dx=pvar["dx"],dy=pvar["dy"],
                latinc=pvar["latinc"],loninc=pvar["loninc"])
        lat_NE = loc[0].values
        lon_NE = loc[1].values

        xbox    = [lon_SW, lon_SE, lon_NE, lon_NW, lon_SW]
        ybox    = [lat_SW, lat_SE, lat_NE, lat_NW, lat_SW]
        ax.plot(xbox, ybox, color=gsLineColor, transform=crs.Geodetic())

        idd     = idom + 1
        dom_txt = f"d0{idd}"
        ax.text(lon_NW, lat_NW, dom_txt, va="top", fontsize=16, 
                color=txFontColor, transform=crs.PlateCarree())

p.show()
