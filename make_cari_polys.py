""" Pull the beach and rocky polygons from the CARI database """
from shared import *

if __name__ == '__main__':
    # slivers removed, no CARI ID
    srcd = Path('/scratch/tws_grace/data/Sea_Level/SFEI/Coastal_Habitats_v2.gdb')
    dstd = srcd.parent / 'CARI_polygons'
    gdf = gpd.read_file(srcd)

    habits_beach = ['Marine Beach Natural Intertidal Non-vegetated']
    habits_rocky = ['Marine Rocky Shore Natural Intertidal', 'Marine Rocky Shore Unnatural Intertidal']
    for i, habits in enumerate([habits_beach, habits_rocky]):
        k = 'Beach' if 'Beach' in habits[0] else 'Rocky' 
        
        gdf_habit = gdf[gdf['clicklabel'].isin(habits)].copy()
        gdf_habit['CARI_id'] = gdf_habit.index
        dst_habit = dstd / f'Cari_{k}.GeoJSON'
        gdf_habit.to_file(dst_habit)
        print (f'Wrote: {dst_habit}')
