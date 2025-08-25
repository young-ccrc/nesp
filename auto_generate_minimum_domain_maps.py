import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle
import os

# === Setup ===
file_path = '/g/data/w28/yk8692/nesp/testing_script'
out_figure_path = '/g/data/w28/yk8692/nesp/figure/minimum_domains'
os.makedirs(out_figure_path, exist_ok=True)

# === Sydney ===
city_name = 'Sydney'
center_lat, center_lon = -33.86, 150.85
dx_km = 1
nx = ny = 400
deg_per_km = 1 / 111
half_width = dx_km * nx / 2 * deg_per_km
half_height = dx_km * ny / 2 * deg_per_km
west, east = center_lon - half_width, center_lon + half_width
south, north = center_lat - half_height, center_lat + half_height

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)
sydney_shp = gccsa[gccsa["GCC_NAME21"].str.contains(city_name, case=False)]

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([145, 155, -38, -30], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')

for geom in sydney_shp.geometry:
    ax.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=2, label=f'Greater {city_name}')

ax.add_patch(Rectangle((west, south), east - west, north - south,
                 edgecolor='red', facecolor='none', lw=2,
                 transform=ccrs.PlateCarree(), label='Minimal Domain'))

ax.set_title(f"Minimum Domain Justification for Extreme Rainfall Simulations\nGreater {city_name} Region", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_sydney.png", dpi=300, bbox_inches='tight')
plt.close()

# === Melbourne ===
city_name = 'Melbourne'
center_lat, center_lon = -37.81, 144.96
dx_km = 1
nx = ny = 400
deg_per_km = 1 / 111
half_width = dx_km * nx / 2 * deg_per_km
half_height = dx_km * ny / 2 * deg_per_km
west, east = center_lon - half_width, center_lon + half_width
south, north = center_lat - half_height, center_lat + half_height

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)
melbourne_shp = gccsa[gccsa["GCC_NAME21"].str.contains(city_name, case=False)]

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([140, 150, -42, -32], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')

for geom in melbourne_shp.geometry:
    ax.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=2, label=f'Greater {city_name}')

ax.add_patch(Rectangle((west, south), east - west, north - south,
                 edgecolor='red', facecolor='none', lw=2,
                 transform=ccrs.PlateCarree(), label='Minimal Domain'))

ax.set_title(f"Minimum Domain Justification for Extreme Rainfall Simulations\nGreater {city_name} Region", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_melbourne.png", dpi=300, bbox_inches='tight')
plt.close()

# === Brisbane ===
city_name = 'Brisbane'
center_lat, center_lon = -27.47, 153.03
dx_km = 1
nx = ny = 400
deg_per_km = 1 / 111
half_width = dx_km * nx / 2 * deg_per_km
half_height = dx_km * ny / 2 * deg_per_km
west, east = center_lon - half_width, center_lon + half_width
south, north = center_lat - half_height, center_lat + half_height

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)
brisbane_shp = gccsa[gccsa["GCC_NAME21"].str.contains(city_name, case=False)]

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([147, 158, -32, -24], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')

for geom in brisbane_shp.geometry:
    ax.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=2, label=f'Greater {city_name}')

ax.add_patch(Rectangle((west, south), east - west, north - south,
                 edgecolor='red', facecolor='none', lw=2,
                 transform=ccrs.PlateCarree(), label='Minimal Domain'))

ax.set_title(f"Minimum Domain Justification for Extreme Rainfall Simulations\nGreater {city_name} Region", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_brisbane.png", dpi=300, bbox_inches='tight')
plt.close()

# === Perth ===
city_name = 'Perth'
center_lat, center_lon = -31.95, 115.86
dx_km = 1
nx = ny = 400
deg_per_km = 1 / 111
half_width = dx_km * nx / 2 * deg_per_km
half_height = dx_km * ny / 2 * deg_per_km
west, east = center_lon - half_width, center_lon + half_width
south, north = center_lat - half_height, center_lat + half_height

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)
perth_shp = gccsa[gccsa["GCC_NAME21"].str.contains(city_name, case=False)]

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([110, 122, -36, -25], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')

for geom in perth_shp.geometry:
    ax.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=2, label=f'Greater {city_name}')

ax.add_patch(Rectangle((west, south), east - west, north - south,
                 edgecolor='red', facecolor='none', lw=2,
                 transform=ccrs.PlateCarree(), label='Minimal Domain'))

ax.set_title(f"Minimum Domain Justification for Extreme Rainfall Simulations\nGreater {city_name} Region", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_perth.png", dpi=300, bbox_inches='tight')
plt.close()

# === Adelaide ===
city_name = 'Adelaide'
center_lat, center_lon = -34.93, 138.6
dx_km = 1
nx = ny = 400
deg_per_km = 1 / 111
half_width = dx_km * nx / 2 * deg_per_km
half_height = dx_km * ny / 2 * deg_per_km
west, east = center_lon - half_width, center_lon + half_width
south, north = center_lat - half_height, center_lat + half_height

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)
adelaide_shp = gccsa[gccsa["GCC_NAME21"].str.contains(city_name, case=False)]

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([132, 145, -40, -28], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')

for geom in adelaide_shp.geometry:
    ax.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=2, label=f'Greater {city_name}')

ax.add_patch(Rectangle((west, south), east - west, north - south,
                 edgecolor='red', facecolor='none', lw=2,
                 transform=ccrs.PlateCarree(), label='Minimal Domain'))

ax.set_title(f"Minimum Domain Justification for Extreme Rainfall Simulations\nGreater {city_name} Region", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_adelaide.png", dpi=300, bbox_inches='tight')
plt.close()

# === Hobart ===
city_name = 'Hobart'
center_lat, center_lon = -42.88, 147.33
dx_km = 1
nx = ny = 400
deg_per_km = 1 / 111
half_width = dx_km * nx / 2 * deg_per_km
half_height = dx_km * ny / 2 * deg_per_km
west, east = center_lon - half_width, center_lon + half_width
south, north = center_lat - half_height, center_lat + half_height

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)
hobart_shp = gccsa[gccsa["GCC_NAME21"].str.contains(city_name, case=False)]

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([140, 150, -47, -39], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')

for geom in hobart_shp.geometry:
    ax.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=2, label=f'Greater {city_name}')

ax.add_patch(Rectangle((west, south), east - west, north - south,
                 edgecolor='red', facecolor='none', lw=2,
                 transform=ccrs.PlateCarree(), label='Minimal Domain'))

ax.set_title(f"Minimum Domain Justification for Extreme Rainfall Simulations\nGreater {city_name} Region", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_hobart.png", dpi=300, bbox_inches='tight')
plt.close()

# === Darwin ===
city_name = 'Darwin'
center_lat, center_lon = -12.46, 130.84
dx_km = 1
nx = ny = 400
deg_per_km = 1 / 111
half_width = dx_km * nx / 2 * deg_per_km
half_height = dx_km * ny / 2 * deg_per_km
west, east = center_lon - half_width, center_lon + half_width
south, north = center_lat - half_height, center_lat + half_height

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)
darwin_shp = gccsa[gccsa["GCC_NAME21"].str.contains(city_name, case=False)]

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([124, 138, -18, -8], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')

for geom in darwin_shp.geometry:
    ax.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=2, label=f'Greater {city_name}')

ax.add_patch(Rectangle((west, south), east - west, north - south,
                 edgecolor='red', facecolor='none', lw=2,
                 transform=ccrs.PlateCarree(), label='Minimal Domain'))

ax.set_title(f"Minimum Domain Justification for Extreme Rainfall Simulations\nGreater {city_name} Region", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_darwin.png", dpi=300, bbox_inches='tight')
plt.close()

# === Canberra ===
city_name = 'Canberra'
center_lat, center_lon = -35.28, 149.13
dx_km = 1
nx = ny = 400
deg_per_km = 1 / 111
half_width = dx_km * nx / 2 * deg_per_km
half_height = dx_km * ny / 2 * deg_per_km
west, east = center_lon - half_width, center_lon + half_width
south, north = center_lat - half_height, center_lat + half_height

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)
canberra_shp = gccsa[gccsa["GCC_NAME21"].str.contains(city_name, case=False)]

fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([142, 155, -40, -30], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')

for geom in canberra_shp.geometry:
    ax.add_geometries([geom], crs=ccrs.PlateCarree(), edgecolor='blue', facecolor='none', linewidth=2, label=f'Greater {city_name}')

ax.add_patch(Rectangle((west, south), east - west, north - south,
                 edgecolor='red', facecolor='none', lw=2,
                 transform=ccrs.PlateCarree(), label='Minimal Domain'))

ax.set_title(f"Minimum Domain Justification for Extreme Rainfall Simulations\nGreater {city_name} Region", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_canberra.png", dpi=300, bbox_inches='tight')
plt.close()

