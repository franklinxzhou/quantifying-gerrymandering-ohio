#!/usr/bin/env python3
"""Create geometry-based node/MCT overlap for the full-weight variant.

Includes ZERO-population Census blocks and positive-area polygon intersections.
This is separate from the POP20 atoms. It is not possible to infer empty-land
intersection from a positive-population-only crosswalk.

Uses the SAME block-to-graph-node resolver as build_mct_atoms.py; preserve the
original graph and GPKG input. Outputs a list of geometrically nonnested places
and an area-overlap incidence TSV. No source file is modified.
"""
import argparse
import csv
import json
from pathlib import Path
import geopandas as gpd
import pandas as pd
from build_mct_atoms import graph_node_lookup, resolve_block_node, selected_geographies


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--graph',type=Path,required=True)
    ap.add_argument('--blocks',type=Path,required=True)
    ap.add_argument('--places',type=Path,required=True)
    ap.add_argument('--townships',type=Path,required=True)
    ap.add_argument('--units',type=Path,required=True)
    ap.add_argument('--cross-ids-out',type=Path,required=True)
    ap.add_argument('--incidence-out',type=Path,required=True)
    ap.add_argument('--min-area-m2',type=float,default=0.01,
                    help='Require intersection area greater than this threshold; 0.01 m2 excludes mere boundary touches')
    ap.add_argument('--chunk-size',type=int,default=3000)
    a=ap.parse_args()
    if a.min_area_m2<0 or a.chunk_size<1:raise ValueError('Invalid geometry threshold/chunk size')
    for target in (a.cross_ids_out,a.incidence_out):
        if target.exists():raise FileExistsError(f'Refusing overwrite: {target}')
    units=list(csv.DictReader(a.units.open(newline='')))
    if not units or 'geometric_status' not in units[0]:
        raise ValueError('Expected population-aware units.csv with geometric_status')
    cross_ids={r['geoid'] for r in units if r['geography']=='place' and r['geometric_status']=='cross'}
    geo_nested={r['geoid'] for r in units if r['geography']=='place' and r['geometric_status']=='nested'}
    valid_towns={r['geoid'] for r in units if r['geography']=='township'}
    if len(cross_ids)!=437 or len(geo_nested)!=488 or len(valid_towns)!=1314:
        raise ValueError(f'Unexpected geography counts: cross={len(cross_ids)}, nested={len(geo_nested)}, towns={len(valid_towns)}')
    towns=selected_geographies(a.townships,kind='township')
    places=selected_geographies(a.places,kind='place')
    towns=towns[towns.GEOID.isin(valid_towns)].copy()
    places=places[places.GEOID.isin(cross_ids)].copy()
    assert len(towns)==len(valid_towns) and len(places)==len(cross_ids)
    blocks=gpd.read_file(a.blocks)
    required={'NAME','COUNTY','ENACTED_UPPER','ENACTED_LOWER','geometry'}
    if not required.issubset(blocks): raise ValueError(f'Missing block columns: {required-set(blocks)}')
    if blocks.crs is None:raise ValueError('Blocks have no CRS')
    nodes, graph_pop=graph_node_lookup(json.loads(a.graph.read_text())['nodes'])
    node_names=[]
    for row in blocks.itertuples(index=False):
        node_names.append(resolve_block_node(row,nodes))
    blocks=blocks[['geometry']].copy()
    blocks['node_name']=node_names
    blocks=blocks.to_crs('EPSG:26917')
    targets=[]
    for kind,frame in (('township',towns),('place',places)):
        tmp=frame.to_crs('EPSG:26917')[['GEOID','geometry']].copy()
        tmp['geography']=kind
        targets.append(tmp)
    targets=gpd.GeoDataFrame(pd.concat(targets,ignore_index=True),crs='EPSG:26917')
    if not targets.geometry.is_valid.all():raise ValueError('Invalid target geometries')
    incidence=set()
    for start in range(0,len(blocks),a.chunk_size):
        part=blocks.iloc[start:start+a.chunk_size]
        # bbox/geometry candidate search, then eliminate boundary-only touches.
        hits=gpd.sjoin(part,targets,how='inner',predicate='intersects')
        for idx,h in hits.iterrows():
            actual=part.loc[idx].geometry.intersection(targets.geometry.iloc[int(h['index_right'])]).area
            if actual>a.min_area_m2:
                incidence.add((h['node_name'],h['geography'],h['GEOID']))
        if start==0 or start//a.chunk_size%20==0:
            print(f'Blocks scanned {min(len(blocks),start+a.chunk_size):,}/{len(blocks):,}; incidences={len(incidence):,}',flush=True)
    a.cross_ids_out.parent.mkdir(parents=True,exist_ok=True)
    a.incidence_out.parent.mkdir(parents=True,exist_ok=True)
    a.cross_ids_out.write_text('\n'.join(sorted(cross_ids))+'\n')
    with a.incidence_out.open('x',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['node_name','geography','geoid'])
        w.writerows(sorted(incidence))
    print('Geometric cross places:',len(cross_ids))
    print('Geometric node–unit incidences:',len(incidence))
    print('Output:',a.incidence_out)

if __name__=='__main__':main()
