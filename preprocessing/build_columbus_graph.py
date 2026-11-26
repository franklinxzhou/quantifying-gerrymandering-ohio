#!/usr/bin/env python3
"""Build non-overwriting Columbus-only Ohio precinct and five-city sampling graphs.

Source graphs must be verified, GEOID-based PLACE_COUNTY graph variants.
Only PLACE and PLACE_COUNTY values are filtered; topology, populations,
county labels, votes, geometry metadata, and protected-city contractions remain
byte-equivalent after JSON decoding. No geopandas or Census preprocessing.

The five-city output is for the sampler. Keep the ORIGINAL unfiltered
OH_pct20.json as the assignment graph for B7 and posthoc analyses.
"""
from __future__ import annotations
import argparse
import json
import os
from collections import Counter
from pathlib import Path
import tempfile

COLUMBUS = "3918000"
CITY_NAMES = {"Akron", "Cincinnati", "Cleveland", "Dayton", "Toledo"}
EXPECTED_UNITS = {"3918000|FRANKLIN": 880830,
                  "3918000|DELAWARE": 14570,
                  "3918000|FAIRFIELD": 10849}
EXPECTED_CITY_POP = 906249
EXPECTED_OH_POP = 11799448


def check(ok: bool, message: str):
    if not ok:
        raise ValueError(message)


def read(path: Path) -> dict:
    check(path.is_file(), f"Missing source graph: {path}")
    obj = json.loads(path.read_text())
    check(isinstance(obj, dict) and isinstance(obj.get('nodes'), list)
          and 'adjacency' in obj and not obj.get('directed',False)
          and not obj.get('multigraph',False), f"Unexpected graph JSON schema: {path}")
    return obj


def strict_places(node: dict, index: int) -> tuple[dict, dict]:
    place, county = node.get('PLACE'), node.get('PLACE_COUNTY')
    check(isinstance(place, dict) and isinstance(county, dict),
          f"Node {index} missing PLACE / PLACE_COUNTY dictionaries")
    recon = Counter()
    for key, pop in place.items():
        check(type(key) is str and len(key)==7 and key.startswith('39') and key.isdigit()
              and type(pop) is int and pop>0,
              f"Node {index}: invalid PLACE entry {key!r}={pop!r}")
    for key,pop in county.items():
        parts=key.split('|') if isinstance(key,str) else []
        check(len(parts)==2 and len(parts[0])==7 and parts[0].startswith('39')
              and parts[0].isdigit() and parts[1] and parts[1]==parts[1].upper()
              and type(pop) is int and pop>0,
              f"Node {index}: invalid PLACE_COUNTY entry {key!r}={pop!r}")
        recon[parts[0]]+=pop
    check(dict(recon)==place,
          f"Node {index}: PLACE and PLACE_COUNTY disagree")
    check(type(node.get('POP20')) is int and sum(place.values())<=node['POP20'],
          f"Node {index}: bad population")
    return place,county


def subset(graph: dict, expected_nodes: int, expected_protected: bool):
    check(len(graph['nodes'])==expected_nodes,
          f"Expected {expected_nodes} nodes; got {len(graph['nodes'])}")
    new = dict(graph)
    updated=[]
    input_place=Counter(); input_pc=Counter()
    output_place=Counter(); output_pc=Counter()
    unique_names=set(); unique_ids=set()
    protected=[]
    for idx,node in enumerate(graph['nodes']):
        name=node.get('NAME'); nid=node.get('id')
        check(name not in unique_names and nid not in unique_ids,
              f"Duplicate node NAME/id {name!r}/{nid!r}")
        unique_names.add(name); unique_ids.add(nid)
        origp, origc=strict_places(node,idx)
        input_place.update(origp); input_pc.update(origc)
        p={COLUMBUS:origp[COLUMBUS]} if COLUMBUS in origp else {}
        pc={k:v for k,v in origc.items() if k.startswith(COLUMBUS+'|')}
        check(sum(p.values())==sum(pc.values()),
              f"Node {idx}: filtering introduced PLACE / PLACE_COUNTY mismatch")
        updated_node=dict(node)
        updated_node['PLACE']=p
        updated_node['PLACE_COUNTY']=pc
        updated.append(updated_node)
        output_place.update(p); output_pc.update(pc)
        if 'PROTECTED_CITY' in node:
            protected.append(node['PROTECTED_CITY'])
    check(set(protected)==(CITY_NAMES if expected_protected else set())
          and len(protected)==(5 if expected_protected else 0),
          f"Unexpected protected-city set: {protected}")
    check(set(output_place)=={COLUMBUS},
          f"Retained non-Columbus place: {output_place}")
    check(dict(output_pc)==EXPECTED_UNITS,
          f"Columbus county totals disagree: {dict(output_pc)} vs {EXPECTED_UNITS}")
    check(output_place[COLUMBUS]==EXPECTED_CITY_POP,
          f"Columbus population mismatch: {output_place[COLUMBUS]}")
    check(sum(node['POP20'] for node in graph['nodes'])==EXPECTED_OH_POP,
          'Input graph population changed or wrong source')
    check(set(graph)==set(new), 'Changed graph keys')
    new['nodes']=updated
    # Strong invariant: only PLACE and PLACE_COUNTY may change, including all
    # original adjacency / edge metadata objects and the five-city contraction.
    for key in graph:
        if key!='nodes':
            check(graph[key]==new[key], f"Changed graph-level key {key}")
    for a,b in zip(graph['nodes'],new['nodes']):
        check(a.keys()==b.keys() and all(a[k]==b[k]
             for k in a if k not in ('PLACE','PLACE_COUNTY')),
             f"Changed non-place field at node {a.get('NAME')}")
    return new, {'nodes':len(updated),'protected_cities':sorted(protected),
                 'total_pop20':EXPECTED_OH_POP,
                 'source_place_population':sum(input_place.values()),
                 'source_place_county_units':len(input_pc),
                 'selected_place_geoid':COLUMBUS,
                 'selected_place_population':output_place[COLUMBUS],
                 'selected_columbus_county_population':dict(sorted(output_pc.items()))}


def write_atomic(path: Path, obj: dict):
    check(not path.exists(),f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp=None
    try:
        with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=path.parent,
            prefix='.'+path.name+'.',suffix='.tmp',delete=False) as f:
            tmp=Path(f.name)
            json.dump(obj,f,ensure_ascii=False,separators=(',',':'),allow_nan=False)
            f.write('\n');f.flush();os.fsync(f.fileno())
        os.replace(tmp,path)
    finally:
        if tmp: tmp.unlink(missing_ok=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--precinct',type=Path,default=Path('data/graph/OH_pct20_place_county.json'))
    ap.add_argument('--five-city',type=Path,default=Path('data/graph/OH_pct20_preserved_five_cities_place_county.json'))
    ap.add_argument('--output-dir',type=Path,default=Path('data/graph'))
    args=ap.parse_args()
    inputs=[(args.precinct,'OH_pct20_columbus_only.json',7928,False),
            (args.five_city,'OH_pct20_preserved_five_cities_columbus_only.json',6975,True)]
    outputs=[args.output_dir / name for _,name,_,_ in inputs]
    report_path=args.output_dir/'OH_columbus_only_graph_audit.json'
    for p in [*outputs,report_path]:check(not p.exists(),f"Output exists, refusing overwrite: {p}")
    bundles=[];reports=[]
    for src,name,n,protected in inputs:
        transformed,report=subset(read(src),n,protected)
        bundles.append((args.output_dir/name,transformed))
        reports.append({'input':str(src),'output':name, **report})
    check(reports[0]['selected_columbus_county_population']==reports[1]['selected_columbus_county_population'],
          'Original and contracted Columbus populations disagree')
    for path,graph in bundles:
        write_atomic(path,graph)
    report={'mode':'Columbus-only sampler: county-specific PLACE_COUNTY edges and county-specific PLACE_COUNTY energy',
            'precinct_assignment_graph_for_B7':'data/graph/OH_pct20.json (unmodified)',
            'graph_variants':reports,
            'assertions':['Original adjacency, edge metadata and non-place node attributes unchanged',
              'Five protected cities remain hard-contracted',
              'Only Columbus PLACE and three Columbus|COUNTY units retained',
              'Both graph variants have identical Columbus and statewide populations']}
    write_atomic(report_path,report)
    print(json.dumps({'result':'PASS','output_dir':str(args.output_dir),
          'selected_columbus_total':EXPECTED_CITY_POP,'county_pieces':EXPECTED_UNITS,
          'source_counts':[row['nodes'] for row in reports],
          'audit':str(report_path)},indent=2))

if __name__=='__main__':main()
