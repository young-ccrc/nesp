import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import box
import os

file_path = '/g/data/w28/yk8692/nesp/testing_script'
out_figure_path = '/g/data/w28/yk8692/nesp/figure/minimum_domains'
os.makedirs(out_figure_path, exist_ok=True)

gccsa_shp = f"{file_path}/GCCSA_2021_AUST_SHP_GDA2020/GCCSA_2021_AUST_GDA2020.shp"
gccsa = gpd.read_file(gccsa_shp)

# === Sydney ===
sydney = gccsa[gccsa["GCC_NAME21"].str.contains("Sydney", case=False)].to_crs(epsg=4326)
buffer_deg = 0.9
minx, miny, maxx, maxy = sydney.total_bounds
west, east = minx - buffer_deg - 0.5, maxx + buffer_deg
south, north = miny - buffer_deg, maxy + buffer_deg
domain_box = box(west, south, east, north)
domain_gdf = gpd.GeoDataFrame(geometry=[domain_box], crs='EPSG:4326')

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([145, 155, -38, -30], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')
sydney.boundary.plot(ax=ax, edgecolor='blue', linewidth=2, transform=ccrs.PlateCarree(), label='Greater Sydney')
domain_gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, transform=ccrs.PlateCarree(), label='Minimum Domain')
ax.set_title("Minimum Domain Justification for Extreme Rainfall Simulations\nGreater Sydney", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_sydney_simplified.png", dpi=300, bbox_inches='tight')
plt.close()

# === Melbourne ===
melbourne = gccsa[gccsa["GCC_NAME21"].str.contains("Melbourne", case=False)].to_crs(epsg=4326)
buffer_deg = 0.9
minx, miny, maxx, maxy = melbourne.total_bounds
west, east = minx - buffer_deg - 0.5, maxx + buffer_deg
south, north = miny - buffer_deg, maxy + buffer_deg
domain_box = box(west, south, east, north)
domain_gdf = gpd.GeoDataFrame(geometry=[domain_box], crs='EPSG:4326')

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([140, 150, -44, -35], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')
melbourne.boundary.plot(ax=ax, edgecolor='blue', linewidth=2, transform=ccrs.PlateCarree(), label='Greater Melbourne')
domain_gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, transform=ccrs.PlateCarree(), label='Minimum Domain')
ax.set_title("Minimum Domain Justification for Extreme Rainfall Simulations\nGreater Melbourne", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_melbourne_simplified.png", dpi=300, bbox_inches='tight')
plt.close()

# === Brisbane ===
brisbane = gccsa[gccsa["GCC_NAME21"].str.contains("Brisbane", case=False)].to_crs(epsg=4326)
buffer_deg = 0.9
minx, miny, maxx, maxy = brisbane.total_bounds
west, east = minx - buffer_deg - 0.5, maxx + buffer_deg
south, north = miny - buffer_deg, maxy + buffer_deg
domain_box = box(west, south, east, north)
domain_gdf = gpd.GeoDataFrame(geometry=[domain_box], crs='EPSG:4326')

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([147, 158, -32, -24], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')
brisbane.boundary.plot(ax=ax, edgecolor='blue', linewidth=2, transform=ccrs.PlateCarree(), label='Greater Brisbane')
domain_gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, transform=ccrs.PlateCarree(), label='Minimum Domain')
ax.set_title("Minimum Domain Justification for Extreme Rainfall Simulations\nGreater Brisbane", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_brisbane_simplified.png", dpi=300, bbox_inches='tight')
plt.close()

# === Perth ===
perth = gccsa[gccsa["GCC_NAME21"].str.contains("Perth", case=False)].to_crs(epsg=4326)
buffer_deg = 0.9
minx, miny, maxx, maxy = perth.total_bounds
west, east = minx - buffer_deg - 0.5, maxx + buffer_deg
south, north = miny - buffer_deg, maxy + buffer_deg
domain_box = box(west, south, east, north)
domain_gdf = gpd.GeoDataFrame(geometry=[domain_box], crs='EPSG:4326')

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([110, 122, -36, -25], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')
perth.boundary.plot(ax=ax, edgecolor='blue', linewidth=2, transform=ccrs.PlateCarree(), label='Greater Perth')
domain_gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, transform=ccrs.PlateCarree(), label='Minimum Domain')
ax.set_title("Minimum Domain Justification for Extreme Rainfall Simulations\nGreater Perth", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_perth_simplified.png", dpi=300, bbox_inches='tight')
plt.close()

# === Adelaide ===
adelaide = gccsa[gccsa["GCC_NAME21"].str.contains("Adelaide", case=False)].to_crs(epsg=4326)
buffer_deg = 0.9
minx, miny, maxx, maxy = adelaide.total_bounds
west, east = minx - buffer_deg - 0.5, maxx + buffer_deg
south, north = miny - buffer_deg, maxy + buffer_deg
domain_box = box(west, south, east, north)
domain_gdf = gpd.GeoDataFrame(geometry=[domain_box], crs='EPSG:4326')

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([132, 145, -40, -28], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')
adelaide.boundary.plot(ax=ax, edgecolor='blue', linewidth=2, transform=ccrs.PlateCarree(), label='Greater Adelaide')
domain_gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, transform=ccrs.PlateCarree(), label='Minimum Domain')
ax.set_title("Minimum Domain Justification for Extreme Rainfall Simulations\nGreater Adelaide", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_adelaide_simplified.png", dpi=300, bbox_inches='tight')
plt.close()

# === Hobart ===
hobart = gccsa[gccsa["GCC_NAME21"].str.contains("Hobart", case=False)].to_crs(epsg=4326)
buffer_deg = 0.9
minx, miny, maxx, maxy = hobart.total_bounds
west, east = minx - buffer_deg - 0.5, maxx + buffer_deg
south, north = miny - buffer_deg, maxy + buffer_deg
domain_box = box(west, south, east, north)
domain_gdf = gpd.GeoDataFrame(geometry=[domain_box], crs='EPSG:4326')

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([140, 150, -47, -39], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')
hobart.boundary.plot(ax=ax, edgecolor='blue', linewidth=2, transform=ccrs.PlateCarree(), label='Greater Hobart')
domain_gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, transform=ccrs.PlateCarree(), label='Minimum Domain')
ax.set_title("Minimum Domain Justification for Extreme Rainfall Simulations\nGreater Hobart", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_hobart_simplified.png", dpi=300, bbox_inches='tight')
plt.close()

# === Darwin ===
darwin = gccsa[gccsa["GCC_NAME21"].str.contains("Darwin", case=False)].to_crs(epsg=4326)
buffer_deg = 0.9
minx, miny, maxx, maxy = darwin.total_bounds
west, east = minx - buffer_deg - 0.5, maxx + buffer_deg
south, north = miny - buffer_deg, maxy + buffer_deg
domain_box = box(west, south, east, north)
domain_gdf = gpd.GeoDataFrame(geometry=[domain_box], crs='EPSG:4326')

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([124, 138, -18, -8], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')
darwin.boundary.plot(ax=ax, edgecolor='blue', linewidth=2, transform=ccrs.PlateCarree(), label='Greater Darwin')
domain_gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, transform=ccrs.PlateCarree(), label='Minimum Domain')
ax.set_title("Minimum Domain Justification for Extreme Rainfall Simulations\nGreater Darwin", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_darwin_simplified.png", dpi=300, bbox_inches='tight')
plt.close()

# === Canberra ===
canberra = gccsa[gccsa["GCC_NAME21"].str.contains("Canberra", case=False)].to_crs(epsg=4326)
buffer_deg = 0.9
minx, miny, maxx, maxy = canberra.total_bounds
west, east = minx - buffer_deg - 0.5, maxx + buffer_deg
south, north = miny - buffer_deg, maxy + buffer_deg
domain_box = box(west, south, east, north)
domain_gdf = gpd.GeoDataFrame(geometry=[domain_box], crs='EPSG:4326')

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([142, 155, -40, -30], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linewidth=1)
ax.add_feature(cfeature.STATES, linestyle='--')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgrey')
canberra.boundary.plot(ax=ax, edgecolor='blue', linewidth=2, transform=ccrs.PlateCarree(), label='Greater Canberra')
domain_gdf.boundary.plot(ax=ax, edgecolor='red', linewidth=2, transform=ccrs.PlateCarree(), label='Minimum Domain')
ax.set_title("Minimum Domain Justification for Extreme Rainfall Simulations\nGreater Canberra", fontsize=14)
ax.legend(loc='lower left')
plt.tight_layout()
plt.savefig(f"{out_figure_path}/minimum_domain_canberra_simplified.png", dpi=300, bbox_inches='tight')
plt.close()

