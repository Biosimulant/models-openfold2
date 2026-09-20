"""Synthetic contract fixtures only; actual NVIDIA calls are recorded separately."""
import copy,importlib.util,io,json
from pathlib import Path
import httpx,pytest,yaml
from Bio.PDB import MMCIFIO
from biosim import BioModule,BioWorld,ExecutionContext,ExecutionPolicy,SignalSpec,make_signal
from biosim.pack import validate_lab_source
LAB=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('openfold2_model',LAB/'owned/models/main/model.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
CTX=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN,run_start=0,run_end=1)
RID='12345678-1234-5678-1234-567812345678'
def pdb():
    lines=[];serial=0
    for number,resname in [(1,'ALA'),(2,'CYS')]:
        for offset,name in enumerate(['N','CA','C','O']):
            serial+=1;lines.append(f'ATOM  {serial:5d} {name:^4} {resname} A{number:4d}    {float(serial):8.3f}{float(offset):8.3f}{0.0:8.3f}  1.00 75.00')
    return '\n'.join(lines)+'\n'
@pytest.fixture
def values():
    text='>query\nAC\n>hit\nAc-\n'
    return {'protein_sequence':'AC','protein_msa':{'query':{'sequence':'AC','sequence_sha256':m.sha(b'AC')},'alignments':[{'database':'Uniref30_2302','format':'a3m','text':text,'artifact':{'sha256':m.sha(text.encode()),'bytes':len(text.encode())}}]},'templates':{'items':[]},'run_options':{'parameters':{'selected_models':[1,2]}},'query_metadata':{'query_id':'synthetic-query','provenance':{'purpose':'contract-only'}}}
def prepared(v):return m.prepare_request(v['protein_sequence'],v['protein_msa'],v['templates'],v['run_options'],v['query_metadata'])
def response():return {'input_id':'synthetic-query','of2_nim_handled_error_message':'no-handled-error','metrics':{},'structures_in_ranked_order':[{'structure':pdb(),'format':'pdb','relaxed':False,'rank_by_confidence':0,'confidence':80.0,'model_param_set':2,'plddt':[80,80],'predicted_aligned_error':[[0,4],[8,1]],'ptm_score':0.5,'max_predicted_aligned_error':31.75,'aligned_confidence_probs':None},{'structure':pdb(),'format':'pdb','relaxed':False,'rank_by_confidence':1,'confidence':70.0,'model_param_set':1}]}
def parse(v,d):
    request,identity=prepared(v);return m.parse_response(m.encode(d),request,identity)
def signals(v):return {k:make_signal(spec=m.NvidiaOpenFold2().inputs()[k],value=x,source='fixture',name=k,emitted_at=0) for k,x in v.items()}
def mock(monkeypatch,handler):
    real=httpx.Client;monkeypatch.setattr(m.httpx,'Client',lambda **kw:real(transport=httpx.MockTransport(handler),**kw));monkeypatch.setenv('NVIDIA_API_KEY','unit-test-key');monkeypatch.setattr(m.time,'sleep',lambda _:None)
def test_manifest_runtime_lab_parity():
    model=m.NvidiaOpenFold2();manifest=yaml.safe_load((LAB/'owned/models/main/model.yaml').read_text());lab=yaml.safe_load((LAB/'lab.yaml').read_text())
    for direction,ports in [('inputs',model.inputs()),('outputs',model.outputs())]:
        entries={p['name']:p for p in manifest['io'][direction]};assert set(entries)==set(ports)
        assert {p['maps_to'] for p in lab['io'][direction]}=={'main.'+k for k in ports}
        for name,port in ports.items():
            for key in ['signal_type','schema','dtype','required','default','emitted_unit']:assert entries[name].get(key)==getattr(port,key)
            assert tuple(entries[name].get('accepted_units',()))==(port.accepted_units or ())
    result=validate_lab_source(LAB);assert result.valid,result.errors
    assert len(lab['models'])==1 and model.execution_policy==ExecutionPolicy.ONCE_BEFORE_RUN
@pytest.mark.parametrize('query',['','ac','AX','A C','AC*','A'*1001,None])
def test_invalid_sequence(values,query):
    values['protein_sequence']=query
    with pytest.raises(ValueError):prepared(values)
@pytest.mark.parametrize('key,value',[('selected_models',[]),('selected_models',[0]),('selected_models',[6]),('selected_models',[True]),('selected_models',[1,1]),('relax_prediction',True),('relax_prediction',None),('use_templates',1),('use_templates',True),('random_seed',1)])
def test_invalid_options(values,key,value):
    values['run_options']['parameters'][key]=value
    with pytest.raises(ValueError):prepared(values)
@pytest.mark.parametrize('mutation',['hash','query','format','columns','paired','duplicate'])
def test_bad_msa(values,mutation):
    entry=values['protein_msa']['alignments'][0]
    if mutation=='hash':entry['artifact']['sha256']='0'*64
    if mutation=='query':values['protein_msa']['query']['sequence']='AA'
    if mutation=='format':entry['format']='fasta'
    if mutation=='columns':entry['text']='>q\nAC\n>h\nACD\n'
    if mutation=='paired':entry['chain_id']='A'
    if mutation=='duplicate':values['protein_msa']['alignments'].append(copy.deepcopy(entry))
    with pytest.raises(ValueError):prepared(values)
def test_explicit_sequence_only(values):
    values['protein_msa']={'query':{},'alignments':[]};request,identity=prepared(values)
    assert 'alignments' not in request and identity['msa_mode']=='sequence-only'
def test_all_ranked_models_and_score_availability(values):
    out=parse(values,response());structures=out['predicted_structure']['structures'];scores=out['confidence_scores']['records']
    assert [s['model_param_set'] for s in structures]==[2,1]
    assert len(structures[0]['residue_mapping'])==2
    assert scores[0]['native']['predicted_aligned_error']==[[0,4],[8,1]]
    assert 'plddt' not in scores[1]['native'] and scores[0]['native']['aligned_confidence_probs'] is None
    assert structures[0]['residue_mapping'][0]['atoms'][0]['native_b_factor']==75
@pytest.mark.parametrize('mutation',['handled_error','no_status','wrong_query_id','missing_model','duplicate_model','wrong_rank','relaxed','rank_order','missing_backbone','wrong_residue','nonfinite_coord','two_chains'])
def test_invalid_prediction(values,mutation):
    d=response();s=d['structures_in_ranked_order'][0]
    if mutation=='handled_error':d['of2_nim_handled_error_message']='GPU failed'
    if mutation=='no_status':del d['of2_nim_handled_error_message']
    if mutation=='wrong_query_id':d['input_id']='someone-else'
    if mutation=='missing_model':d['structures_in_ranked_order'].pop()
    if mutation=='duplicate_model':s['model_param_set']=1
    if mutation=='wrong_rank':s['rank_by_confidence']=1
    if mutation=='relaxed':s['relaxed']=True
    if mutation=='rank_order':s['confidence']=60
    if mutation=='missing_backbone':s['structure']='\n'.join(s['structure'].splitlines()[1:])
    if mutation=='wrong_residue':s['structure']=s['structure'].replace('CYS','ALA')
    if mutation=='nonfinite_coord':s['structure']=s['structure'][:30]+'     nan'+s['structure'][38:]
    if mutation=='two_chains':s['structure']=s['structure'].replace('CYS A','CYS B')
    with pytest.raises(ValueError):parse(values,d)
@pytest.mark.parametrize('key,value',[('plddt',[80]),('plddt',[80,101]),('predicted_aligned_error',[[0,1]]),('predicted_aligned_error',[[0,-1],[2,0]]),('ptm_score',1.2),('confidence',float('nan')),('max_predicted_aligned_error',[]),('aligned_confidence_probs',[[[]]])])
def test_invalid_confidence(values,key,value):
    d=response();d['structures_in_ranked_order'][0][key]=value
    with pytest.raises(ValueError):parse(values,d)
def test_mmcif_structure_support(values):
    structure=m.PDBParser(PERMISSIVE=False,QUIET=True).get_structure('test',io.StringIO(pdb()));writer=MMCIFIO();writer.set_structure(structure);out=io.StringIO();writer.save(out)
    d=response();d['structures_in_ranked_order'][0].update({'structure':out.getvalue(),'format':'mmcif'})
    assert parse(values,d)['predicted_structure']['structures'][0]['format']=='mmcif'
def template_item():
    structure=m.PDBParser(PERMISSIVE=False,QUIET=True).get_structure('test',io.StringIO(pdb()));writer=MMCIFIO();writer.set_structure(structure);out=io.StringIO();writer.save(out)
    text=out.getvalue()+'\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n1 1 ALA\n1 2 CYS\n#\n'
    return {'structure':text,'format':'mmcif','sha256':m.sha(text.encode()),'name':'synthetic','source':'synthetic-test-only','provenance':{'purpose':'contract test; no scientific template claim'}}
def test_explicit_template_preserved(values):
    item=template_item();values['templates']['items']=[item];values['run_options']['parameters']['use_templates']=True
    request,identity=prepared(values)
    assert request['explicit_templates'][0]['structure']==item['structure']
    assert request['explicit_templates'][0]['rank']==-1
    assert identity['templates'][0]['native_selected_label_chain']=='A'
@pytest.mark.parametrize('mutation',['hash','compressed','missing_fields','disabled','chain'])
def test_invalid_template(values,mutation):
    item=template_item();values['templates']['items']=[item];values['run_options']['parameters']['use_templates']=True
    if mutation=='hash':item['sha256']='0'*64
    if mutation=='compressed':item['format']='mmcif.gz'
    if mutation=='missing_fields':item['structure']='data_empty\n';item['sha256']=m.sha(item['structure'].encode())
    if mutation=='disabled':values['run_options']['parameters']['use_templates']=False
    if mutation=='chain':item['structure']=item['structure'].replace(' A ',' B ');item['sha256']=m.sha(item['structure'].encode())
    with pytest.raises(ValueError):prepared(values)
@pytest.mark.parametrize('status',[401,422,429,500,503])
def test_native_failure_no_retries(values,monkeypatch,status):
    calls=[]
    def handler(req):calls.append(req.method);return httpx.Response(status,json={'error':'synthetic provider failure'},headers={'nvcf-reqid':RID})
    mock(monkeypatch,handler);model=m.NvidiaOpenFold2()
    with pytest.raises(RuntimeError):model.execute(signals(values),context=CTX)
    assert calls==['POST'] and model.snapshot()['last_http_status']==status
    assert model.snapshot()['provider_request_id']==RID and 'unit-test-key' not in json.dumps(model.snapshot())
    assert not model.get_outputs()
def test_async_poll_identity(values,monkeypatch):
    calls=[]
    def handler(req):
        calls.append(req.method)
        if len(calls)==1:return httpx.Response(202,json={},headers={'nvcf-reqid':RID})
        assert str(req.url)==m.POLL+RID
        return httpx.Response(200,json=response(),headers={'nvcf-reqid':RID})
    mock(monkeypatch,handler);model=m.NvidiaOpenFold2();out=model.execute(signals(values),context=CTX)
    assert calls==['POST','GET'] and out['run_provenance']['provider']['request_id']==RID
def test_exact_resume(values,monkeypatch):
    calls=[]
    def handler(req):calls.append(req.method);return httpx.Response(200,json=response(),headers={'nvcf-reqid':RID})
    mock(monkeypatch,handler);request,_=prepared(values);model=m.NvidiaOpenFold2(resume_request_id=RID,resume_request_sha256=m.sha(m.encode(request)))
    model.execute(signals(values),context=CTX);assert calls==['GET']
def test_wrong_resume_hash(values,monkeypatch):
    mock(monkeypatch,lambda req:pytest.fail('Must not send request'))
    with pytest.raises(ValueError):m.NvidiaOpenFold2(resume_request_id=RID,resume_request_sha256='0'*64).execute(signals(values),context=CTX)
def test_missing_credential(values,monkeypatch):
    monkeypatch.delenv('NVIDIA_API_KEY',raising=False)
    with pytest.raises(ValueError,match='NVIDIA_API_KEY'):m.NvidiaOpenFold2().execute(signals(values),context=CTX)
def test_response_bound(values,monkeypatch):
    mock(monkeypatch,lambda req:httpx.Response(200,content=b' '*1100))
    with pytest.raises(RuntimeError):m.NvidiaOpenFold2(max_response_bytes=1024).execute(signals(values),context=CTX)
def test_handled_error_does_not_emit_predictions(values,monkeypatch):
    data=response();data['of2_nim_handled_error_message']='synthetic failed prediction'
    mock(monkeypatch,lambda req:httpx.Response(200,json=data));model=m.NvidiaOpenFold2()
    with pytest.raises(RuntimeError):model.execute(signals(values),context=CTX)
    assert not model.get_outputs()
def test_bioworld_once_all_terminal_ports(values,monkeypatch):
    calls=[]
    def handler(req):calls.append(req.method);return httpx.Response(200,json=response(),headers={'nvcf-reqid':RID})
    mock(monkeypatch,handler)
    class Source(BioModule):
        execution_policy=ExecutionPolicy.ONCE_BEFORE_RUN
        def outputs(self):return {'protein_sequence':SignalSpec.scalar(dtype='str',emitted_unit='1'),'protein_msa':SignalSpec.record(schema=m.MSA,emitted_unit='1'),'templates':SignalSpec.record(schema=m.TEMPLATES,emitted_unit='1'),'run_options':SignalSpec.record(schema=m.OPTIONS,emitted_unit='1'),'query_metadata':SignalSpec.record(schema=m.META,emitted_unit='1')}
        def execute(self,inputs,*,context):return values
    model=m.NvidiaOpenFold2();world=BioWorld(communication_step=1);world.add_biomodule('source',Source());world.add_biomodule('main',model)
    for name in values:world.connect('source.'+name,'main.'+name)
    world.run(duration=1)
    assert calls==['POST'] and set(model.get_outputs())==set(model.outputs())
    assert all(signal.emitted_at==0 for signal in model.get_outputs().values())
