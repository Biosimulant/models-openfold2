"""NVIDIA OpenFold2 adapter; hosted runtime/scientific qualification pending."""
import hashlib,io,json,math,os,re,time,uuid
from datetime import datetime,timezone
import httpx
from Bio.PDB import PDBParser,MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Data.PDBData import protein_letters_3to1
from biosim import BioModule,ExecutionPolicy,SignalSpec
ENDPOINT='https://health.api.nvidia.com/v1/biology/openfold/openfold2/predict-structure-from-msa-and-template'
POLL='https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/'
AA='ARNDCQEGHILKMFPSTWYV'
MSA={'query':'json','alignments':'json'}
TEMPLATES={'items':'json'}
OPTIONS={'parameters':'json'}
META={'query_id':'str','provenance':'json'}
STRUCTURES={'query':'json','structures':'json'}
SCORES={'records':'json','semantics':'json'}
METRICS={'native':'json','semantics':'json'}
PROVENANCE={'request':'json','provider':'json','inputs':'json','backend_identity':'json','limitations':'json'}
def require(ok,message):
    if not ok:raise ValueError(message)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def encode(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def exact(value,keys,label):require(isinstance(value,dict) and set(value)==set(keys),label+' schema mismatch')
def strict_json(raw):
    def pairs(items):
        result={}
        for key,value in items:
            require(key not in result,'Duplicate native JSON key');result[key]=value
        return result
    data=json.loads(raw,object_pairs_hook=pairs);encode(data);return data
def finite(value):return type(value) in (int,float) and math.isfinite(value)

def check_a3m(text,sequence,max_rows=20000):
    require(isinstance(text,str) and 0<len(text.encode())<=5_000_000,'Invalid alignment size')
    rows=[];parts=[];header=None
    for line in text.splitlines():
        if not line:continue
        if line.startswith('>'):
            if header is not None:
                require(bool(parts),'Empty MSA row');rows.append(''.join(parts));require(len(rows)<max_rows,'Too many MSA rows')
            header=line[1:];parts=[]
            require(header.strip() and not any(ord(c)<32 and c!='\t' for c in header),'Invalid MSA header')
            require(not header.startswith(('ss_pred','ss_conf','ss_dssp','sa_')),'Annotation tracks need an explicit adapter')
        else:require(header is not None and not any(c.isspace() for c in line),'Malformed MSA row');parts.append(line)
    require(header is not None and parts,'Empty MSA');rows.append(''.join(parts));require(len(rows)<=max_rows,'Too many MSA rows')
    require(rows[0]==sequence,'MSA query differs from protein input')
    alphabet=set(AA+'XBZJUO'+(AA+'XBZJUO').lower()+'-.')
    require(all(set(s)<=alphabet and len(''.join(c for c in s if not c.islower() and c!='.'))==len(sequence) for s in rows),'Invalid A3M match columns/alphabet')
    return len(rows)

def prepare_request(sequence,msa,templates,options,metadata):
    require(isinstance(sequence,str) and 1<=len(sequence)<=1000 and re.fullmatch('['+AA+']+',sequence),'Expected 1–1000 uppercase standard protein residues')
    exact(metadata,META,'Query metadata');require(isinstance(metadata['query_id'],str) and len(metadata['query_id'])<=128 and not any(ord(c)<32 for c in metadata['query_id']),'Invalid query ID')
    require(isinstance(metadata['provenance'],dict),'Query provenance must be an object');encode(metadata)
    query={'id':metadata['query_id'] or 'sha256:'+sha(sequence.encode()),'sequence':sequence,'sha256':sha(sequence.encode()),'length':len(sequence),'provenance':metadata['provenance']}
    exact(options,OPTIONS,'Run options');params=options['parameters'];require(isinstance(params,dict) and set(params)<= {'selected_models','use_templates','relax_prediction'},'Unsupported option')
    settings={'selected_models':[1,2,3,4,5],'use_templates':False,'relax_prediction':False,**params}
    selected=settings['selected_models'];require(isinstance(selected,list) and 1<=len(selected)<=5 and all(type(v) is int and 1<=v<=5 for v in selected) and len(set(selected))==len(selected),'Invalid selected parameter sets')
    require(settings['relax_prediction'] is False,'Current NIM does not support relaxation');require(type(settings['use_templates']) is bool,'Invalid template toggle')
    request={'sequence':sequence,'input_id':query['id'],**settings}
    exact(msa,MSA,'MSA collection');entries=msa['alignments'];require(isinstance(entries,list) and len(entries)<=10 and isinstance(msa['query'],dict),'Invalid MSA collection')
    alignment_info=[]
    if entries:
        require(msa['query'].get('sequence')==sequence and msa['query'].get('sequence_sha256')==sha(sequence.encode()),'MSA query identity/hash mismatch')
        require(msa['query'].get('mode')!='paired','Paired MSA requires an explicit monomer chain selector')
        request['alignments']={}
        for entry in entries:
            require(isinstance(entry,dict) and not entry.get('chain_id') and entry.get('format')=='a3m','Only monomer A3M entries supported')
            db=entry.get('database');require(isinstance(db,str) and 0<len(db)<=128 and db not in request['alignments'],'Invalid/duplicate database key')
            text=entry.get('text');count=check_a3m(text,sequence);artifact=entry.get('artifact')
            require(isinstance(artifact,dict) and artifact.get('sha256')==sha(text.encode()) and artifact.get('bytes')==len(text.encode()),'MSA artifact hash/size mismatch')
            request['alignments'][db]={'a3m':{'alignment':text,'format':'a3m'}}
            alignment_info.append({'database':db,'rows':count,'sha256':sha(text.encode()),'bytes':len(text.encode())})
    else:require(msa['query']=={},'Empty MSA collection cannot carry a silently ignored query')
    template_info=prepare_templates(templates,settings['use_templates'],request)
    require(len(encode(request))<=10_000_000,'Request exceeds client byte bound')
    return request,{'query':query,'alignments':alignment_info,'templates':template_info,'msa_mode':'provided' if entries else 'sequence-only'}

def prepare_templates(envelope,enabled,request):
    exact(envelope,TEMPLATES,'Templates');items=envelope['items'];require(isinstance(items,list) and len(items)<=10,'Invalid template count')
    require(enabled==bool(items),'Template enablement and provided items must agree')
    result=[];native=[]
    for item in items:
        exact(item,('structure','format','sha256','name','source','provenance'),'Template item')
        text=item['structure'];require(item['format']=='mmcif' and isinstance(text,str) and 0<len(text.encode())<=2_000_000,'Uncompressed mmCIF text required within client bounds')
        require(item['sha256']==sha(text.encode()),'Template checksum mismatch')
        require(isinstance(item['name'],str) and len(item['name'])<=128 and isinstance(item['source'],str) and len(item['source'])<=2048 and isinstance(item['provenance'],dict),'Invalid template provenance')
        d=MMCIF2Dict(io.StringIO(text));required=('_atom_site.label_asym_id','_atom_site.label_seq_id','_atom_site.label_comp_id','_atom_site.label_atom_id','_atom_site.Cartn_x','_atom_site.Cartn_y','_atom_site.Cartn_z','_entity_poly_seq.num','_entity_poly_seq.mon_id')
        require(all(k in d for k in required),'Incomplete mmCIF identity/coordinate fields')
        n=len(d['_atom_site.label_asym_id']);require(n>0 and all(len(d[k])==n for k in required[:7]),'Inconsistent mmCIF atom columns')
        require('A' in d['_atom_site.label_asym_id'],'Native template processing requires label_asym_id A')
        require(all(math.isfinite(float(v)) for k in required[4:7] for v in d[k]),'Nonfinite template coordinate')
        native.append({'structure':text,'format':'mmcif','name':item['name'],'source':item['source'],'rank':-1})
        result.append({'name':item['name'],'source':item['source'],'sha256':item['sha256'],'bytes':len(text.encode()),'label_chains':sorted(set(d['_atom_site.label_asym_id'])),'native_selected_label_chain':'A','provenance':item['provenance'],'query_to_template_mapping':'Native internal alignment; not exposed or independently established by wrapper.'})
    if native:request['explicit_templates']=native
    return result

def parse_structure(text,fmt,sequence):
    require(isinstance(text,str) and bool(text) and len(text.encode())<=5_000_000,'Invalid structure artifact size')
    if fmt=='pdb':structure=PDBParser(PERMISSIVE=False,QUIET=True).get_structure('prediction',io.StringIO(text))
    else:
        require(fmt=='mmcif','Unsupported native structure format');structure=MMCIFParser(QUIET=True).get_structure('prediction',io.StringIO(text))
    models=list(structure);require(len(models)==1,'Expected one model per ranked structure artifact')
    chains=list(models[0]);require(len(chains)==1,'Monomer output must contain one identified chain')
    residues=list(chains[0]);require(len(residues)==len(sequence),'Output residue count differs from query')
    mapped=[];observed=[]
    for position,residue in enumerate(residues,1):
        require(residue.id[0]==' ' and not residue.is_disordered(),'Unexpected hetero/disordered residue')
        aa=protein_letters_3to1.get(residue.resname);require(aa==sequence[position-1],'Output residue sequence differs from query');observed.append(aa)
        require({'N','CA','C','O'}<=set(residue.child_dict),'Incomplete expected backbone')
        atoms=[]
        for atom in residue:
            require(not atom.is_disordered(),'Unresolved alternate atom in prediction')
            xyz=[float(v) for v in atom.coord];require(all(math.isfinite(v) for v in xyz),'Nonfinite predicted coordinate')
            bf=atom.get_bfactor();occ=atom.get_occupancy();require(finite(bf) and finite(occ) and 0<=occ<=1,'Invalid occupancy/B-factor')
            atoms.append({'name':atom.id,'serial':atom.serial_number,'element':atom.element,'xyz_angstrom':xyz,'occupancy':occ,'native_b_factor':bf})
        mapped.append({'query_position':position,'amino_acid':aa,'chain_id':chains[0].id,'author_residue_number':residue.id[1],'insertion_code':residue.id[2].strip(),'atoms':atoms})
    require(''.join(observed)==sequence,'Output/query identity mismatch')
    return mapped

def check_confidence(item,n):
    confidence=item.get('confidence');require(confidence is None or finite(confidence),'Nonfinite ranking confidence')
    scores={key:item[key] for key in ('confidence','plddt','predicted_aligned_error','max_predicted_aligned_error','aligned_confidence_probs','ptm_score') if key in item}
    plddt=item.get('plddt')
    if plddt is not None:require(isinstance(plddt,list) and len(plddt)==n and all(finite(v) and 0<=v<=100 for v in plddt),'Invalid pLDDT dimension/domain')
    pae=item.get('predicted_aligned_error')
    if pae is not None:require(isinstance(pae,list) and len(pae)==n and all(isinstance(row,list) and len(row)==n and all(finite(v) and v>=0 for v in row) for row in pae),'Invalid PAE dimensions/domain')
    ptm=item.get('ptm_score')
    if ptm is not None:require(finite(ptm) and 0<=ptm<=1,'Invalid predicted TM score')
    upper=item.get('max_predicted_aligned_error')
    if upper is not None:
        vals=upper if isinstance(upper,list) else [upper];require(bool(vals) and all(finite(v) and v>=0 for v in vals),'Invalid PAE bin upper-bound field')
    probs=item.get('aligned_confidence_probs')
    if probs is not None:
        require(isinstance(probs,list) and len(probs)==n and all(isinstance(row,list) and len(row)==n and all(isinstance(bins,list) and len(bins)==64 and all(finite(v) and 0<=v<=1 for v in bins) for bins in row) for row in probs),'Invalid PAE-bin probability tensor')
        scores['probability_sum_diagnostics']={'minimum':min(sum(bins) for row in probs for bins in row),'maximum':max(sum(bins) for row in probs for bins in row),'acceptance':'Normalization tolerance not scientifically reviewed; range/shape checks only.'}
    return scores

def parse_response(raw,request,identity):
    data=strict_json(raw);require(isinstance(data,dict) and set(data)<= {'structures_in_ranked_order','input_id','metrics','of2_nim_handled_error_message'},'Unexpected native response schema')
    require(data.get('of2_nim_handled_error_message')=='no-handled-error','Native handled-error or missing success status')
    require(data.get('input_id')==request['input_id'],'Returned input identity mismatch')
    items=data.get('structures_in_ranked_order');require(isinstance(items,list) and len(items)==len(request['selected_models']) and bool(items),'Missing requested structure predictions')
    required={'structure','format','relaxed','rank_by_confidence','confidence','model_param_set'}
    optional={'plddt','predicted_aligned_error','max_predicted_aligned_error','aligned_confidence_probs','ptm_score'}
    structures=[];scores=[];params=[];rank_values=[]
    for rank,item in enumerate(items):
        require(isinstance(item,dict) and required<=set(item) and set(item)<=required|optional,'Unexpected native structure schema')
        require(type(item['rank_by_confidence']) is int and item['rank_by_confidence']==rank,'Native rank is not contiguous zero-based order')
        require(item['relaxed'] is False,'Unexpected relaxed structure')
        model=item['model_param_set'];require(type(model) is int and model in request['selected_models'],'Unresolved or unexpected native parameter-set identity');params.append(model)
        mapping=parse_structure(item['structure'],item['format'],request['sequence']);native_scores=check_confidence(item,len(request['sequence']))
        if item['confidence'] is not None:rank_values.append(item['confidence'])
        body=item['structure'].encode();structures.append({'rank':rank,'model_param_set':model,'format':item['format'],'structure':item['structure'],'relaxed':False,'artifact':{'handling':'inline UTF-8 content','sha256':sha(body),'bytes':len(body)},'residue_mapping':mapping,'coordinate_unit':'angstrom'})
        scores.append({'rank':rank,'model_param_set':model,'available_fields':list(native_scores),'native':native_scores})
    require(set(params)==set(request['selected_models']) and len(set(params))==len(params),'Requested model identities incomplete/duplicated')
    require(rank_values==sorted(rank_values,reverse=True),'Ranking confidence contradicts returned order')
    metrics=data.get('metrics');require(metrics is None or isinstance(metrics,dict),'Invalid native metrics')
    return {'predicted_structure':{'query':identity['query'],'structures':structures},'confidence_scores':{'records':scores,'semantics':{'confidence':'Native ranking confidence; exact formula/scale unverified. Not measured accuracy.','plddt':'Predicted local distance difference test, 0–100 score, query residue order.','predicted_aligned_error':'Expected aligned error in angstrom, query-by-query matrix; not necessarily symmetric.','ptm_score':'Predicted TM score, dimensionless; not measured structural similarity.','native_b_factor':'Preserved coordinate-file field; no automatic confidence or experimental uncertainty conversion.','availability':'Absent and null fields preserved; hosted release identity unresolved.'}},'inference_metrics':{'native':metrics,'semantics':{'relaxation_time_in_seconds':'seconds if present; current requested relaxation false','other_fields':'Units/meaning unverified unless explicitly documented','temporal':'Finite prediction, no biological time evolution'}},'native_response':raw.decode('utf-8')}

class NvidiaOpenFold2(BioModule):
    execution_policy=ExecutionPolicy.ONCE_BEFORE_RUN
    def __init__(self,observation_timeout_s=600.0,max_response_bytes=20000000,resume_request_id=None,resume_request_sha256=None):
        require(finite(observation_timeout_s) and 1<=observation_timeout_s<=1200,'Invalid observation bound')
        require(type(max_response_bytes) is int and 1024<=max_response_bytes<=50000000,'Invalid response byte bound')
        require((resume_request_id is None)==(resume_request_sha256 is None),'Resume requires UUID and original request SHA-256')
        if resume_request_sha256 is not None:require(isinstance(resume_request_sha256,str) and re.fullmatch('[0-9a-f]{64}',resume_request_sha256),'Invalid resume digest')
        self.timeout=float(observation_timeout_s);self.max_bytes=max_response_bytes
        self.resume_id=str(uuid.UUID(resume_request_id)) if resume_request_id else None;self.resume_sha=resume_request_sha256;self.state={'status':'not_started'}
    def inputs(self):
        return {'protein_sequence':SignalSpec.scalar(dtype='str',accepted_units=('1',),required=True),'protein_msa':SignalSpec.record(schema=MSA,accepted_units=('1',),required=False,default={'query':{},'alignments':[]}),'templates':SignalSpec.record(schema=TEMPLATES,accepted_units=('1',),required=False,default={'items':[]}),'run_options':SignalSpec.record(schema=OPTIONS,accepted_units=('1',),required=False,default={'parameters':{}}),'query_metadata':SignalSpec.record(schema=META,accepted_units=('1',),required=False,default={'query_id':'','provenance':{}})}
    def outputs(self):
        return {'predicted_structure':SignalSpec.record(schema=STRUCTURES,emitted_unit='1'),'confidence_scores':SignalSpec.record(schema=SCORES,emitted_unit='1'),'inference_metrics':SignalSpec.record(schema=METRICS,emitted_unit='1'),'run_provenance':SignalSpec.record(schema=PROVENANCE,emitted_unit='1'),'native_response':SignalSpec.scalar(dtype='str',emitted_unit='1')}
    def snapshot(self):return dict(self.state)
    def execute(self,inputs,*,context):
        require('protein_sequence' in inputs,'Missing protein query')
        values={k:v.value for k,v in inputs.items()}
        values.setdefault('protein_msa',{'query':{},'alignments':[]});values.setdefault('templates',{'items':[]});values.setdefault('run_options',{'parameters':{}});values.setdefault('query_metadata',{'query_id':'','provenance':{}})
        request,identity=prepare_request(values['protein_sequence'],values['protein_msa'],values['templates'],values['run_options'],values['query_metadata'])
        body=encode(request);digest=sha(body);require(not self.resume_id or self.resume_sha==digest,'Resume input digest mismatch')
        key=os.environ.get('NVIDIA_API_KEY');require(isinstance(key,str) and bool(key.strip()),'Secure NVIDIA_API_KEY not provisioned on worker')
        started=time.monotonic();deadline=started+self.timeout
        self.state={'status':'observing','request_sha256':digest,'provider_request_id':self.resume_id}
        method,url=('GET',POLL+self.resume_id) if self.resume_id else ('POST',ENDPOINT)
        try:
            with httpx.Client(follow_redirects=False,trust_env=False) as client:
                while True:
                    remaining=deadline-time.monotonic();require(remaining>0,'Observation deadline reached')
                    headers={'Authorization':'Bearer '+key,'Content-Type':'application/json','NVCF-POLL-SECONDS':str(min(30,max(1,int(remaining)-1)))}
                    with client.stream(method,url,headers=headers,content=body if method=='POST' else None,timeout=min(60,remaining)) as response:
                        rid=response.headers.get('nvcf-reqid') or response.headers.get('nvcf-request-id')
                        if rid:
                            rid=str(uuid.UUID(rid));require(self.state['provider_request_id'] in (None,rid),'Provider request identity changed');self.state['provider_request_id']=rid
                        safe_headers={k:v for k,v in response.headers.items() if k.lower() in {'date','content-type','nvcf-reqid','nvcf-request-id','retry-after'}}
                        raw=b''
                        for chunk in response.iter_bytes():
                            require(len(raw)+len(chunk)<=self.max_bytes and time.monotonic()<deadline,'Response observation limit exceeded');raw+=chunk
                        status=response.status_code
                    self.state.update({'last_http_status':status,'response_bytes':len(raw),'response_sha256':sha(raw),'safe_headers':safe_headers})
                    if status==202:
                        require(bool(self.state['provider_request_id']),'Accepted response has no resumable UUID')
                        method,url='GET',POLL+self.state['provider_request_id'];time.sleep(min(1,max(0,deadline-time.monotonic())));continue
                    require(status==200,'NVIDIA HTTP failure');require(key.encode() not in raw,'Unexpected credential in native output');break
            self.state['phase']='output_contract';outputs=parse_response(raw,request,identity)
        except Exception:
            self.state['status']='observation_or_contract_failed'
            raise RuntimeError('OpenFold2 observation/contract failed; inspect safe snapshot and preserve the existing request ID. No automatic resubmission.') from None
        self.state['status']='completed'
        outputs['run_provenance']={'request':{'body':request,'sha256':digest},'provider':{'endpoint':ENDPOINT,'request_id':self.state['provider_request_id'],'safe_headers':safe_headers,'response_sha256':sha(raw),'response_bytes':len(raw),'observed_at':datetime.now(timezone.utc).isoformat(),'wall_seconds':time.monotonic()-started,'resumed':bool(self.resume_id)},'inputs':{'values':values,'checked_identity':identity},'backend_identity':{'documented_hosted_schema':'2.1.0','self_hosted_target':'NIM 2.6.0','deployed_version':None,'checkpoint_digests':None,'container_digest':None,'provider_gpu':None},'limitations':['Exact hosted identity, managed runtime and independent scientific acceptance remain unresolved.','No measured accuracy or biological/clinical claim follows from confidence or template-conditioned reconstruction.','Native template-to-query mapping is internal; compressed template conventions need additional integration evidence.']}
        return outputs
