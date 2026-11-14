"""Retain real pre-Dotplot render records under an exact coordinate-only delta."""
import collections,hashlib,json,re,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def sha(value):return hashlib.sha256(value).hexdigest()
def blob(name):return subprocess.check_output(['git','show','3a56988:'+name],cwd=ROOT)
def dotplot_record():
    path=ROOT/'tools/site-audit/dotplot-rendering.json'
    if not path.exists():return None
    record=json.loads(path.read_text())
    assert record['status']=='complete' and record['other_expressions_and_signatures_identical']
    assert record['changed_function']=='ggchord_prepare_dotplot'
    for name,digest in record['source_sha256'].items():assert sha((ROOT/name).read_bytes())==digest,'Dotplot render source moved: '+name
    assert record['engine_sha256']=={str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in (ROOT/'R').glob('*.R')}
    counts=collections.Counter()
    for name,digest in record['raw_source_verification']['files'].items():
        data=(ROOT/name).read_bytes();assert sha(data)==digest
        for line in data.decode().splitlines():
            if line and not line.startswith('#'):
                fields=line.split('\t');counts[fields[0],fields[1]]+=1
    assert sum(counts.values())==32 and len(counts)==6
    assert not any(a==b or (b,a) in counts for a,b in counts)

    assert all(record['checks'][key] for key in ('alignment_and_axis_fields_identical','reversed_axis_retains_32','inverse_mapping','explicit_position_scale_wins'))
    return record

def validate_home_sources(record):
    dot=dotplot_record()
    for name,digest in record['source_sha256'].items():
        current=(ROOT/name).read_bytes()
        if sha(current)==digest:continue
        if dot and name=='R/coord-dotplot.R':
            assert sha(blob(name))==digest
            continue
        if dot and name=='tools/site-flagships.R':
            previous=blob(name)
            assert sha(previous)==digest
            chunks=lambda value:dict(re.findall(r'^# ---- ([a-z]+) ----\n(.*?)(?=^# ---- [a-z]+ ----)',value.decode(),re.M|re.S))
            a=chunks(previous);b=chunks(current)
            assert {k for k in a.keys()|b.keys() if a.get(k)!=b.get(k)}=={'flagshipdot'}
            continue
        if name=='website/manual-figures.json':
            metadata=record['post_render_metadata'];previous=blob(name)
            assert sha(previous)==metadata['after_sha256']
            old=json.loads(previous);new=json.loads(current)
            if dot:
                assert [r for r in old if r['id']!='flagdotplot']==[r for r in new if r['id']!='flagdotplot']
                before=next(r for r in old if r['id']=='flagdotplot');after=next(r for r in new if r['id']=='flagdotplot')
                assert {k:v for k,v in before.items() if k not in {'source','width','height','dpi'}}=={k:v for k,v in after.items() if k not in {'source','width','height','dpi'}}
                item=dot['figures']['flagdotplot']
                assert all(after[k]==item[k] for k in ('width','height','dpi'))
            else:assert json.loads(current)==old
            next(r for r in old if r['id']=='flaglinear').update(metadata['previous_fields'])
            assert sha((json.dumps(old,ensure_ascii=False,indent=2)+'\n').encode())==digest
            continue
        raise AssertionError('Rendering source changed outside reviewed boundary: '+name)
    return dot
