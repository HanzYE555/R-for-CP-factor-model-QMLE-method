"""Prepare an ignored reproduction workspace from tracked, compressed data.
Never fits models or runs optimizers. Original audited scripts use their original paths.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=ROOT/'.reproduce')
args=parser.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
manifest=json.loads((ROOT/'DATA_MANIFEST.json').read_text())
for item in manifest['datasets']:
    stage=item.get('stage_path')
    if not stage:continue
    target=out/stage;target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()==item['sha256']:continue
    with target.open('wb') as dest:
        if item['storage']=='plain':
            with (ROOT/item['files'][0]).open('rb') as src:shutil.copyfileobj(src,dest)
        elif item['storage']=='gzip':
            with gzip.open(ROOT/item['files'][0],'rb') as src:shutil.copyfileobj(src,dest)
        else:
            for i,part in enumerate(item['parts']):
                with gzip.open(ROOT/part['path'],'rb') as src:
                    header=src.readline()
                    if i==0:dest.write(header)
                    shutil.copyfileobj(src,dest)
    assert hashlib.sha256(target.read_bytes()).hexdigest()==item['sha256'],stage
(out/'outputs/Id_Prob_QMLE_CP_Model_JASA_S23/additions/lambda').mkdir(parents=True,exist_ok=True)
print('Prepared:',out)
print('No model fitting or geometric optimization has been run.')
